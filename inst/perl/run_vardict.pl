#!/usr/bin/perl -w

### AUTHORS ###################################################################
# Name			Date
# Erle			2017-08-18

### DESCRIPTION ###############################################################
# Call variants with VarDict. 
# Creates an output directory named vardict under the $sample_id directory
#
# References:
#	https://github.com/AstraZeneca-NGS/VarDict
#	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4914105/
#
# TO DO:
#	- figure out what teststrandbias.R did and add it in R
#	- figure out how to invoke var2vcf_valid.pl script
#	- check output after removing .dict file
#	- figure out what's going on with the VCF formatting errors in step 1
#	- update Annovar protocol (and associated post-processing) to remove icgc21 with GRCh38
#
# QUESTIONS:
#	- why is the minimum AF threshold so low?
# 	- do we need all the ANNOVAR output files? 

### LIBRARIES #################################################################
use strict;
use Getopt::Long;
use YAML qw(LoadFile);
use File::Basename;

### COMMAND LINE VARIABLES ####################################################
#
# tumour_bam: path to tumour BAM file
# sample_id: name of sample (for file naming purposes)
# config_file: path to config file
# paired: flag indicating if analysis should be paired
# normal_bam: path to normal BAM file (required for paired analysis)
# proton: flag indicating that proton sequencing was used
# output_directory: path to directory where results should be saved (defaults to directory with tumour_bam)
# output_filename: Filename of final output file (defaults to ${sample_id}.vcf)
# code_directory: directory where code should be saved (defaults to output_directory)
# log_directory: directory where logs should be saved (defaults to output_directory)
# job_dependencies: jobs the script depends on 
# job_name: name of job (useful for letting R control job dependencies)
# job_group: group job belongs to

GetOptions(
	# Required
	"tumour_bam=s" => \my $tumour_bam,
	"sample_id=s" => \my $sample_id,
	"config_file=s" => \my $config_file,
	# Optional
	"normal_bam=s" => \my $normal_bam,
	"paired" => \(my $paired = 0),
	"proton" => \( my $proton = 0), 
	"output_directory=s" => \my $output_directory, 
	"output_filename=s" => \my $output_filename,
	"code_directory=s" => \my $code_directory,
	"log_directory=s" => \my $log_directory,
	"job_dependencies=s{1,}" => \my @job_dependencies,
	"job_name=s" => \my $job_name,
	"job_group=s" => \my $job_group,
	);

die "$0 requires the sample_id argument\n" unless $sample_id;
die "$0 requires the tumour_bam argument\n" unless $tumour_bam;
die "$0 requires the config_file argument\n" unless $config_file;


if($paired && not $normal_bam) { 
	print "Paired flag is on but no normal sample BAM has been provided.\n";
	exit;
}

# if no output_directory, set to directory containing tumour BAM file
if( not $output_directory ) {
	$output_directory = dirname($tumour_bam);
}

if( not $code_directory ) { 
	$code_directory = $output_directory;
}

if( not $log_directory ) {
	$log_directory = $output_directory;
}

# if no output filename has been supplied, set to SAMPLE_ID.vcf
if( not $output_filename ) {
	$output_filename = "${sample_id}.vcf";
}

# if no job name provided, set to mutect_SAMPLE_ID.vcf
if( not $job_name ) { 
	my $job_name = "vardict_${sample_id}";
}

### MAIN ######################################################################

my $config = LoadFile($config_file);

my $reference_build = $config->{reference_build};
my $reference_genome = $config->{reference_genome}->{ $reference_build };

# Allele frequency in config file is set for filtering
# set a low number here to catch all variants, and then filter them out in post-processing step
my $af_threshold = 0.001;

# get path to corresponding .dict file
(my $reference_dict = $reference_genome) =~ s{\.[^.]+$}{};
$reference_dict = $reference_dict . ".dict";

# cluster settings
my $num_cpu = $config->{num_cpu};
my $project_code = $config->{project_code};
my $cluster_scheduler = $config->{cluster_scheduler};
my $cluster_settings = "";

# header of output VCF – label samples with TUMOUR and NORMAL rather than sample ID
my $sample_name = "TUMOUR";

# software tools
my $annovar_path = $config->{annovar_path};
my $annovar_database = $config->{annovar_database}->{ $reference_build };

my $gatk_jar = $config->{gatk_jar};
my $picard_jar = $config->{picard_jar};

# module load versions
my $java = "java/$config->{java_version}";
my $vardictjava = "vardictjava/$config->{vardictjava_version}";

# set buildver parameter based on config file
my $buildver = "hg19";
my $gencode_version = "V19";

if( "grch38" eq lc $reference_build ) { 
	$buildver = "hg38";
	$gencode_version = "V24";
}

# adjust ANNOVAR and Picard calls based on desired gene annotation method
my $annotation_db = "refGene";

if( "gencode" eq lc ${config}->{annotation_method} ) { 
	$annotation_db = "wgEncodeGencodeBasic" . $gencode_version;
}

my $vardict_path = $config->{vardict_path};

# set read base quality threshold
# PGM tends to underestimate base qualities, so VarDict recommends setting a lower threshold for proton runs
my $min_base_quality = 25;

if( $proton ) { 
	$min_base_quality = 15;
} 

## JOB DEPENDENCIES

my $job_dependency_string = '';
my $orphan_flag = ''; # if job_dependencies is non-empty, set this to -ti to automatically kill orphan jobs

if( scalar(@job_dependencies) > 0 ) {
  if ( $cluster_scheduler eq "LSF" ) {
  	my @parsed_job_dependencies = map { "done(\"$_\")" } @job_dependencies;
  
  	$job_dependency_string = join( " && ", @parsed_job_dependencies );
  	$orphan_flag = '-ti';
  } elsif ( $cluster_scheduler eq "PBS" ) {
    $job_dependency_string = join( ",", @job_dependencies );
  }
}

## JOB GROUP
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
#BSUB -n 5
#BSUB -q normal
#BSUB -P ${project_code}
#BSUB -w '${job_dependency_string}' ${orphan_flag}
#BSUB -W 124:00
${job_group_string}
END
} elsif ( $cluster_scheduler eq "PBS" ) {
  $cluster_settings = <<END;
#!/bin/bash
#PBS -P ${project_code}
#PBS -N ${job_name}
#PBS -o ${log_directory}/${job_name}.out
#PBS -e ${log_directory}/${job_name}.err 
#PBS -n 5
#PBS -q normal
#PBS -hold_jid ${job_dependency_string}
#PBS -l cput=0:124:00
${job_group_string}
END
} else {
  $cluster_settings = ""
}

my $sample_vardict_script = "${code_directory}/${job_name}.sh";


open(OUT,">$sample_vardict_script") or die "Cannot open output script $sample_vardict_script $!\n";

print OUT <<EOF;
${cluster_settings}

module load $java $vardictjava;

mkdir -p ${output_directory} ${log_directory} ${code_directory}

EOF

if( $paired ) { 
	print OUT <<EOF;

# STEP 1: run VarDict
#	- input: bam file (must be indexed)
#	- output: variants in VCF file
#		* SAMPLE_ID.vcf
VarDict -G $reference_genome \\
-f $af_threshold \\
-N \"TUMOUR|NORMAL\" \\
-b \"${tumour_bam}|${normal_bam}\" \\
-z -c 1 -S 2 -th 5 -E 3 -g 4 -q ${min_base_quality} \\
$config->{target_panel} \\
| $vardict_path/testsomatic.R \\
| $vardict_path/var2vcf_paired.pl \\
-N \"TUMOUR|NORMAL\"  -f $af_threshold > ${output_directory}/all_variants.vcf

awk '\$0 !~ /STATUS=Germline/' ${output_directory}/all_variants.vcf > ${output_directory}/${output_filename} 


EOF

} else {
	print OUT <<EOF;

# STEP 1: run VarDict
#	- input: bam file (must be indexed)
#	- output: variants in VCF file
#		* SAMPLE_ID.vcf
VarDict -G $reference_genome \\
-f ${af_threshold} \\
-N \"TUMOUR\" \\
-b ${tumour_bam} \\
-z -c 1 -S 2 -th 5 -E 3 -g 4 -q ${min_base_quality} \\
$config->{target_panel} \\
| $vardict_path/teststrandbias.R \\
| $vardict_path/var2vcf_valid.pl \\
-N \"TUMOUR\" -f $af_threshold > ${output_directory}/${output_filename}

EOF

}

close(OUT);
if ($cluster_scheduler eq "LSF") {
  system("bsub < $sample_vardict_script");
} elsif ($cluster_scheduler eq "PBS") {
  system("qsub $sample_vardict_script");
}

