#!/usr/bin/perl -w

### AUTHORS ###################################################################
#
# NAME			DATE			DESCRIPTION
# Erle 			2017-09-11		Initial development
#
### DESCRIPTION ###############################################################
#
# Run coverage quality control
#
# TO DO: 
#	- set output_directory automatically if full_bam and ontarget_bam are in the same directory

### LIBRARIES #################################################################
use strict;
use Getopt::Long;
use YAML qw(LoadFile);
use File::Basename;


### COMMAND LINE ARGUMENTS ####################################################
#
# full_bam: Path to BAM file with all aligned reads
# sample_id: ID of sample to be aligned
# config_file: path to config file 
# output_directory: path to directory where output files should be saved
# ontarget_bam_filename: 
# code_directory: directory where code should be saved (defaults to output_directory)
# log_directory: directory where logs should be saved (defaults to output_directory)
# job_dependencies: jobs the script depends on 
# job_group: group job belongs to

GetOptions(
	# required
	"bam_file=s" => \my $bam_file,
	"sample_id=s" => \ my $sample_id,
	"config_file=s" => \my $config_file,
	"output_directory=s" => \my $output_directory, 
	# optional
	"ontarget_bam_filename=s" => \my $ontarget_bam_filename,
	"code_directory=s" => \my $code_directory,
	"log_directory=s" => \my $log_directory,
	"job_dependencies=s{1,}" => \my @job_dependencies,
	"job_name=s" => \my $job_name,
	"job_group=s" => \my $job_group,
	);

die "$0 requires the bam_file argument\n" unless $bam_file;
die "$0 requires the output_directory argument\n" unless $output_directory;
die "$0 requires the sample_id argument\n" unless $sample_id;
die "$0 requires the config_file argument\n" unless $config_file;


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

my $input_directory = dirname($bam_file);

my $reference_build = $config->{reference_build};
my $reference_genome = $config->{reference_genome}->{$reference_build};
my $germline_vcf = $config->{ germline_vcf };

my $target_panel = $config->{target_panel};


# cluster settings
my $num_cpu = $config->{num_cpu};
my $project_code = $config->{project_code};
my $cluster_scheduler = $config->{cluster_scheduler};
my $cluster_settings = "";

# Picard
my $picard_jar = $config->{picard_jar};
my $java_memory = $num_cpu*12.8;

# module load version numbers
my $bwa = "bwa/$config->{bwa_version}";
my $bedtools = "bedtools/$config->{bedtools_version}";
my $java = "java/$config->{java_version}";
my $samtools = "samtools/$config->{samtools_version}";
my $picardtools = "picard-tools/$config->{picardtools_version}";

my $gatk_jar = $config->{ gatk_jar };

if( not $job_name ) { 
	$job_name = "target_qc_${sample_id}";
}

## JOB GROUP
# submitting -g and and empty string does not seem to work, so include whole line as variable
my $job_group_string = '';
if( $job_group && $cluster_scheduler eq "LSF" ) { 
	$job_group_string = "#BSUB -g /${job_group}";
} elsif ( $job_group && $cluster_scheduler eq "PBS" ) {
  # Is there such a thing? Does it matter?
  $job_group_string = '';
}

## QSUB/BSUB SETTINGS
if( $cluster_scheduler eq "LSF" ) {
  $cluster_settings = <<END;
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -o ${log_directory}/${job_name}.out
#BSUB -e ${log_directory}/${job_name}.err 
#BSUB -n 8
#BSUB -q normal
#BSUB -P ${project_code}
#BSUB -W 168:00
${job_group_string}
END
} elsif ( $cluster_scheduler eq "PBS" ) {
  $cluster_settings = <<END;
#!/bin/bash
#PBS -P ${project_code}
#PBS -N ${job_name}
#PBS -o ${log_directory}/${job_name}.out
#PBS -e ${log_directory}/${job_name}.err 
#PBS -n 8
#PBS -q normal
#PBS -l cput=0:168:00
${job_group_string}
END
} else {
  $cluster_settings = ""
}


### MAIN ######################################################################

my $sample_coverage_qc_script="${code_directory}/${job_name}.sh";


open(OUT,">$sample_coverage_qc_script") or die "Cannot open output script $sample_coverage_qc_script $!\n";
print OUT <<EOF;
${cluster_settings}


module load $bedtools $samtools $java $picardtools

mkdir -p ${output_directory} ${code_directory} ${log_directory}

# STEP 1: indexing (speeds up downstream steps, and is required for variant calling with VarDict)
#	- input: sorted bam file with alignments
#	- output: .bai file
#		* BAM_FILE.bai
samtools index ${bam_file}

# STEP 2: replace read groups with a single read group
#	- input: sorted bam
# 	- output: bam file with all reads assigned to the same read group
#		* SAMPLE_ID.sorted.bam
java -Xmx${java_memory}g -jar $picard_jar \\
AddOrReplaceReadGroups \\
INPUT=${bam_file} \\
OUTPUT=${output_directory}/${sample_id}.sorted.bam \\
RGID=$sample_id \\
RGLB=$sample_id \\
RGPL=$config->{sequencing_platform} \\
RGPU=$config->{barcode} \\
RGSM=$sample_id \\
CREATE_INDEX=true

# STEP 2.5: GATK BQSR
# - input: single read group bam file
# - output: recalibrated bam file
#		* SAMPLE_ID.recal.bam

# BaseRecalibrator

java -Xmx12g -jar $gatk_jar \\
BaseRecalibrator \\
-I ${output_directory}/${sample_id}.sorted.bam \\
-R $reference_genome \\
--known-sites ${germline_vcf} \\
-O ${output_directory}/${sample_id}.recalibration.table
  
# ApplyBQSR

java -Xmx12g -jar $gatk_jar \\
ApplyBQSR \\
-R $reference_genome \\
-I ${output_directory}/${sample_id}.sorted.bam \\
--bqsr-recal-file ${output_directory}/${sample_id}.recalibration.table \\
-O ${output_directory}/${sample_id}.recal.bam

# STEP 3: restrict to reads aligned to target regions
#	- input: single read group bam file
#	- output: bam containing only ontarget reads
#		* SAMPLE_ID.sorted.bam.ontarget.bam
bedtools intersect \\
-abam ${output_directory}/${sample_id}.recal.bam \\
-b $target_panel \\
-u -wa \\
> ${output_directory}/${ontarget_bam_filename}

# STEP 4: indexing ontarget reads
#	- input: bam file with ontarget reads
#	- output: index
#		* ONTARGET_BAM_FILENAME.bai
samtools index \\
${output_directory}/${ontarget_bam_filename}

# STEP 4: flagstat (read mapping statistics)
#	- input: bam file before and after filtering off-target reads
#	- output: one text file containing read mapping statistics before and after filtering off-target reads
#		* SAMPLE_ID.sorted.bam.stats
echo 'Flagstat before filtering off-target reads' \\
> ${output_directory}/${sample_id}.sorted.bam.stats \\
&& \\
samtools flagstat \\
${bam_file} \\
>> ${output_directory}/${sample_id}.sorted.bam.stats \\

echo 'Flagstat after filtering off-target reads' \\
>> \\
${output_directory}/${sample_id}.sorted.bam.stats \\
&& \\
samtools flagstat ${output_directory}/${ontarget_bam_filename} \\
>> ${output_directory}/${sample_id}.sorted.bam.stats

# STEP 5: bedtools basic coverage report
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

# STEP 6: bedtools "histogram of coverage" reports for features in target panel
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

# STEP 7: calculate insert size metrics
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
  system("bsub < $sample_coverage_qc_script");
} elsif ($cluster_scheduler eq "PBS") {
  system("qsub $sample_coverage_qc_script");
}

