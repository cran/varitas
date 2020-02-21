#!/usr/bin/perl -w

### AUTHORS ###################################################################
# Name			Date			Description
# Erle			2017-08-17		Initial development
# Erle			2017-09-08		Refactored to accept more low-level arguments

### DESCRIPTION ###############################################################
#
# Run alignment and associated tools on a single sample.
# 
## QUESTIONS 
# 	- why include paths to tools as variables rather than just module loading?
#	- should panel be a command line argument or included in the config file? 
# 	- should the output BAM file name be a parameter? What about the names of the other output files? 
#
# TO DO: 
#	- intelligent warnings

### LIBRARIES #################################################################
use strict;
use Getopt::Long;
use YAML qw(LoadFile);
use File::Basename;


### COMMAND LINE ARGUMENTS ####################################################
#
# fastq: paths to FastQ files to be aligned
# sample_id: ID of sample to be aligned
# config_file: path to config file 
# output_directory: path to directory where output files should be saved
#	- required unless all FASTQ files are in the same directory (in which case output defaults to this directory too.
# ontarget_bam_filename: filename of main output file – BAM with ontarget reads
# code_directory: directory where code should be saved (defaults to output_directory)
# log_directory: directory where logs should be saved (defaults to output_directory)
# job_dependencies: jobs the script depends on 
# job_group: group job belongs to

GetOptions(
	# required
	"fastq=s{1,}" => \my @fastq_files,
	"sample_id=s" => \ my $sample_id,
	"config_file=s" => \my $config_file,
	# optional
	"output_directory=s" => \my $output_directory, 
	"ontarget_bam_filename=s" => \my $ontarget_bam_filename,
	"code_directory=s" => \my $code_directory,
	"log_directory=s" => \my $log_directory,
	"job_dependencies=s{1,}" => \my @job_dependencies,
	"job_name=s" => \my $job_name,
	"job_group=s" => \my $job_group
	);

die "$0 requires the fastq argument\n" unless @fastq_files;
die "$0 requires the sample_id argument\n" unless $sample_id;
die "$0 requires the config_file argument\n" unless $config_file;
die "$0 can accept at most two FASTQ files\n" unless (scalar @fastq_files < 3);

my @fastq_dirnames = map { dirname($_); } @fastq_files;
 
# if all FASTQs are in the same directory, and no output directory has been specified,
# make output directory same as the input
my %string = map { $_, 1 } @fastq_dirnames;
if( keys %string == 1 && not $output_directory ) { 
	$output_directory = $fastq_dirnames[0];
} 

die "$0 requires the output_directory argument" unless $output_directory;

if( not $code_directory ) { 
	$code_directory = $output_directory;
}

if( not $log_directory ) {
	$log_directory = $output_directory;
}

## NAME OF OUTPUT FILE
if( not $ontarget_bam_filename ) {
	$ontarget_bam_filename = "${sample_id}.sorted.bam.ontarget.bam";
}

# idea: pass "frequently changing parameters" as command line
# include config file for more complicated 
my $config = LoadFile($config_file);


my $reference_build = $config->{reference_build};
my $reference_genome = $config->{reference_genome}->{ $reference_build };
my $germline_vcf = $config->{ germline_vcf };

my $target_panel = $config->{target_panel};

# cluster settings
my $num_cpu = $config->{num_cpu};
my $project_code = $config->{project_code};
my $cluster_scheduler = $config->{cluster_scheduler};
my $cluster_settings = "";

# Picard
my $picard_jar = $config->{picard_jar};
my $sequencing_platform = $config->{sequencing_platform};
my $barcode = $config->{barcode};
my $java_memory = $num_cpu*12.8;


# module load version numbers
my $bwa = "bwa/$config->{bwa_version}";
my $bedtools= "bedtools/$config->{bedtools_version}";
my $java = "java/$config->{java_version}";
my $samtools = "samtools/$config->{samtools_version}";
my $picardtools = "picard-tools/$config->{picardtools_version}";
my $fastqc = "fastqc/$config->{fastqc_version}";

my $gatk_jar = $config->{ gatk_jar };

## JOB GROUP
# submitting -g and and empty string does not seem to work, so include whole line as variable
my $job_group_string = '';
if( $job_group && $cluster_scheduler eq "LSF" ) { 
	$job_group_string = "#BSUB -g /${job_group}";
} elsif ( $job_group && $cluster_scheduler eq "PBS") {
  $job_group_string = "#PBS -g /${job_group}";
}

## QSUB/BSUB SETTINGS
if($cluster_scheduler eq "LSF") {
  $cluster_settings = <<END;
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -o ${log_directory}/${job_name}.out
#BSUB -e ${log_directory}/${job_name}.err 
#BSUB -n 5
#BSUB -q normal
#BSUB -P $project_code
#BSUB -W 168:00
${job_group_string}
END
} elsif ($cluster_scheduler eq "PBS") {
  $cluster_settings = <<END;
#!/bin/bash
#PBS -P $project_code
#PBS -N ${job_name}
#PBS -o ${log_directory}/${job_name}.out
#PBS -e ${log_directory}/${job_name}.err 
#PBS -n 5
#PBS -q normal
#PBS -l cput=0:168:00
${job_group_string}
END
} else {
  $cluster_settings = ""
}

### MAIN ######################################################################

if( not $job_name ) { 
	$job_name = "bwa_${sample_id}";
}

my $sample_alignment_script="${code_directory}/${job_name}.sh";

my $fastq_string = join(" ", @fastq_files);


open(OUT,">$sample_alignment_script") or die "Cannot open output script $sample_alignment_script $!\n";
print OUT <<EOF;
${cluster_settings}

module load $bedtools $samtools $bwa $java $picardtools $fastqc;

mkdir -p ${output_directory} ${code_directory} ${log_directory}

# STEP 1: quality control
#	- input: fastq files 
#	- output: html report and zip file for each input fastq file
fastqc $fastq_string \\
--outdir ${output_directory}

# STEP 2: run alignment with bwa mem
#	- input: fastq files and reference genome
#	- output: sam file with alignments
#		* fq.aln.sam
bwa mem  -M \\
-t $num_cpu \\
$reference_genome \\
$fastq_string \\
> ${output_directory}/${sample_id}.fq.aln.sam

# STEP 3: convert to bam
#	- input: sam file
#	- output: bam file
#		* fq.aln.bam
samtools view -bS \\
${output_directory}/${sample_id}.fq.aln.sam \\
> ${output_directory}/${sample_id}.fq.aln.bam

# STEP 4: sorting
# 	- input: bam file with alignments
#	- output: sorted bam file 
# 		* fq.aln.sorted.bam
samtools sort \\
${output_directory}/${sample_id}.fq.aln.bam \\
-o ${output_directory}/${sample_id}.fq.aln.sorted.bam

# STEP 5: indexing (speeds up downstream steps, and is required for variant calling with VarDict)
#	- input: sorted bam file with alignments
#	- output: .bai file
#		* SAMPLE_ID.sort.bam.bai
samtools index \\
${output_directory}/${sample_id}.fq.aln.sorted.bam

# STEP 6: replace read groups with a single read group
#	- input: sorted bam
# 	- output: bam file with all reads assigned to the same read group
#		* SAMPLE_ID.sorted.bam
java -Xmx${java_memory}g -jar $picard_jar \\
AddOrReplaceReadGroups \\
INPUT=${output_directory}/${sample_id}.fq.aln.sorted.bam \\
OUTPUT=${output_directory}/${sample_id}.sorted.bam \\
RGID=$sample_id \\
RGLB=$sample_id \\
RGPL=$sequencing_platform \\
RGPU=$barcode \\
RGSM=$sample_id \\
CREATE_INDEX=true

# STEP 6.5: GATK BQSR
# - input: single read group bam file
# - output: recalibrated bam file
#		* SAMPLE_ID.recal.bam

# BaseRecalibrator

java -Xmx12g -jar $gatk_jar BaseRecalibrator \\
-I ${output_directory}/${sample_id}.sorted.bam \\
-R $reference_genome \\
--known-sites ${germline_vcf} \\
-O ${output_directory}/${sample_id}.recalibration.table
  
# ApplyBQSR

java -Xmx12g -jar $gatk_jar ApplyBQSR \\
-R $reference_genome \\
-I ${output_directory}/${sample_id}.sorted.bam \\
--bqsr-recal-file ${output_directory}/${sample_id}.recalibration.table \\
-O ${output_directory}/${sample_id}.recal.bam

# STEP 7: restrict to reads aligned to target regions
#	- input: single read group bam file
#	- output: bam containing only ontarget reads
#		* SAMPLE_ID.sorted.bam.ontarget.bam
bedtools intersect \\
-abam ${output_directory}/${sample_id}.recal.bam \\
-b $target_panel \\
-u -wa \\
> ${output_directory}/${ontarget_bam_filename}

samtools index \\
${output_directory}/${ontarget_bam_filename}

# STEP 8: flagstat (read mapping statistics)
#	- input: bam file before and after filtering off-target reads
#	- output: one text file containing read mapping statistics before and after filtering off-target reads
#		* SAMPLE_ID.sorted.bam.stats
echo 'Flagstat before filtering off-target reads' \\
> ${output_directory}/${sample_id}.sorted.bam.stats \\
&& \\
samtools flagstat \\
${output_directory}/${sample_id}.sorted.bam \\
>> ${output_directory}/${sample_id}.sorted.bam.stats \\

echo 'Flagstat after filtering off-target reads' \\
>> \\
${output_directory}/${sample_id}.sorted.bam.stats \\
&& \\
samtools flagstat ${output_directory}/${ontarget_bam_filename} \\
>> ${output_directory}/${sample_id}.sorted.bam.stats

# STEP 9: indexing of ontarget reads
#	- input: bam file of ontarget reads
#	- output: .bai file
#		* SAMPLE_ID.sorted.bam.ontarget.bam.bai
samtools index \\
${output_directory}/${ontarget_bam_filename}

# STEP 10: bedtools basic coverage report
#	- input: bam file of ontarget reads
# 	- output: bed file of features in target panel and basic coverage statistics
#		* SAMPLE_ID.sorted.bam.ontarget.txt
#		* SAMPLE_ID.sorted.bam.ontarget.sort.txt (sorted)
#	- TO DO: fix sorting, currently treating chromosome as character
bedtools coverage \\
-b ${output_directory}/${ontarget_bam_filename} \\
-a $target_panel \\
> ${output_directory}/${sample_id}.sorted.bam.ontarget.txt \\
&& \\
sort -k1,1 -k2,2 -k3,3 ${output_directory}/$sample_id.sorted.bam.ontarget.txt \\
> ${output_directory}/$sample_id.sorted.bam.ontarget.sort.txt 

# STEP 11: bedtools "histogram of coverage" reports for features in target panel
# 	- input: bam file of ontarget reads
#	- output: full coverage report and coverage report subsetted to only include summary for all features
#		* coverage.report (full)
#		* all.coverage.report (subsetted version)
bedtools coverage -hist \\
-a $target_panel \\
-b ${output_directory}/${ontarget_bam_filename} \\
> ${output_directory}/${sample_id}.coverage.report

awk '\$1=="all"' \\
${output_directory}/${sample_id}.coverage.report \\
> ${output_directory}/${sample_id}.all.coverage.report

# STEP 12: calculate insert size metrics
#	- input: bam file with reads assigned to a single read group (including off-target reads!)
#	- output: frequency table of insert sizes
#		* insert_size.metrics.txt
#	- NOTE: removed TMP_DIR here, might have to re-add it if that causes problems
java -Xmx${java_memory}g -jar $picard_jar \\
CollectInsertSizeMetrics  \\
INPUT=${output_directory}/${sample_id}.sorted.bam \\
OUTPUT=${output_directory}/${sample_id}.insert_size.metrics.txt \\
H=${output_directory}/${sample_id}.insert_size.metrics.pdf \\
REFERENCE_SEQUENCE=$reference_genome \\
MAX_RECORDS_IN_RAM=5000000 \\
VALIDATION_STRINGENCY=LENIENT

EOF
close(OUT);
if ($cluster_scheduler eq "LSF") {
  system("bsub < $sample_alignment_script");
} elsif ($cluster_scheduler eq "PBS") {
  system("qsub $sample_alignment_script");
}