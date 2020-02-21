#!/usr/bin/perl -w

### AUTHORS ###################################################################
# Name			Date
# Erle			2017-08-29
# Adam      2018-12-05

### DESCRIPTION ###############################################################
#
# Run MuTect on paired tumour normal samples
#
# TO DO:
#	- parameterize whether to use dbSNP and COSMIC files

### LIBRARIES #################################################################
use strict;
use Getopt::Long;
use YAML qw(LoadFile);
use File::Basename;

### COMMAND LINE ARGUMENTS ####################################################
#
# tumour_bam: path to tumour BAM file
# normal_bam: path to normal BAM file
# output_directory: path to directory where results should be saved (defaults to directory with tumour_bam)
# sample_id: name of sample (for file naming purposes)
# config_file: path to config file
# job_dependencies: jobs the script depends on 
# output_directory: path to directory where results should be saved (defaults to directory with tumour_bam)
# output_filename: Filename of final output file (defaults to ${sample_id}.passed.ontarget.vcf)
# code_directory: directory where code should be saved (defaults to output_directory)
# log_directory: directory where logs should be saved (defaults to output_directory)
# job_name: name of job (useful for letting R control job dependencies)
# job_group: group job belongs to

GetOptions(
	# Required
	"tumour_bam=s" => \my $tumour_bam,
	"normal_bam=s" => \my $normal_bam,
	"sample_id=s" => \my $sample_id,
	"config_file=s" => \my $config_file,
	# Optional
	"paired" => \(my $paired = 0), 
	"output_directory=s" => \my $output_directory, 
	"output_filename=s" => \my $output_filename,
	"code_directory=s" => \my $code_directory,
	"log_directory=s" => \my $log_directory,
	"job_dependencies=s{1,}" => \my @job_dependencies,
	"job_name=s" => \my $job_name, 
	"job_group=s" => \my $job_group,
	"PON" => \(my $PON = 0),
	);

die "$0 requires the sample_id argument\n" unless $sample_id;
die "$0 requires the tumour_bam argument\n" unless $tumour_bam;
die "$0 requires the config_file argument\n" unless $config_file;


if($paired && not defined $normal_bam) { 
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
	$output_filename = "${sample_id}.passed.ontarget.vcf";
}

# if no job name provided, set to mutect_SAMPLE_ID.vcf
if( not $job_name ) { 
	my $job_name = "mutect_${sample_id}";
}


### MAIN ######################################################################

# idea: pass "frequently changing parameters" as command line
# include config file for more complicated 
my $config = LoadFile($config_file);
my $reference_build = $config->{ reference_build };
my $reference_genome = $config->{ reference_genome }->{ $reference_build };
my $germline_vcf = $config->{ germline_vcf };
my $target_panel = $config->{ target_panel };
my $gatk_jar = $config->{ gatk_jar };

# cluster settings
my $num_cpu = $config->{num_cpu};
my $project_code = $config->{project_code};
my $cluster_scheduler = $config->{cluster_scheduler};
my $cluster_settings = "";

# module load version numbers
my $java = "java/$config->{java_version}";
my $samtools = "samtools/$config->{samtools_version}";

## JOB DEPENDENCIES
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

my $normal_bam_string = '';
if( $paired ) {
	$normal_bam_string = "-I ${normal_bam}";
}

my $normal_rg_string = '';
if( $paired ) {
	$normal_rg_string = "-normal ${sample_id}-NORMAL";
}

## QSUB/BSUB SETTINGS
if( $cluster_scheduler eq "LSF" ) {
  $cluster_settings = <<END;
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -P ${project_code}
#BSUB -w '${job_dependency_string}' ${orphan_flag}
#BSUB -W 168:00  
#BSUB -q normal  
#BSUB -n 5
#BSUB -R "span[ptile=5]"
#BSUB -o ${log_directory}/${job_name}.out     
#BSUB -e ${log_directory}/${job_name}.err 
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
#PBS -l cput=0:168:00
${job_group_string}
END
} else {
  $cluster_settings = ""
}
 
### BUILD SCRIPT ##############################################################
# treat paired and unpaired mode separately (see if you can merge afterwards)
#
my $sample_mutect_script="${code_directory}/${job_name}.sh";

open(OUT,">$sample_mutect_script") or die "Cannot open output script $!\n";

# HEADER
print OUT <<EOF;
${cluster_settings}

module load $samtools $java;

mkdir -p ${output_directory} ${log_directory} ${code_directory}

# STEP 1: run MuTect
# - input: tumour (and possibly normal) BAM
# - output: variant file
#	* SAMPLE_ID.vcf
java -Xmx64g \\
-Djava.io.tmpdir=${output_directory} -jar $gatk_jar \\
Mutect2 \\
-I ${tumour_bam} ${normal_bam_string} \\
-tumor ${sample_id} ${normal_rg_string} \\
--germline-resource ${germline_vcf} \\
-R $reference_genome \\
-O ${output_directory}/$sample_id.vcf \\

# STEP 2: Filter variants
#  input: variant VCF 
#  output: Variants with filtering flags set
java -Xmx64g -Djava.io.tmpdir=${output_directory} -jar $gatk_jar \\
	FilterMutectCalls \\
  -V ${output_directory}/$sample_id.vcf \\
  -O ${output_directory}/$sample_id.filter.vcf
  
# STEP 3: remove all variants not flagged as "passed"
#  input: variant VCF 
#  output: VCF with variants that passed filter
java -Xmx64g -Djava.io.tmpdir=${output_directory} -jar $gatk_jar \\
	SelectVariants \\
	-R $reference_genome \\
	-V ${output_directory}/$sample_id.filter.vcf \\
	-O ${output_directory}/$sample_id.passed.vcf \\
	--exclude-filtered 

# STEP 4: remove off-target variants
# - input: VCF with passed variants
# - output: ontarget and passed variant VCF
#	* OUTPUT_FILENAME
java -Xmx64g -jar $gatk_jar \\
   SelectVariants \\
   -R $reference_genome \\
   -V ${output_directory}/$sample_id.passed.vcf \\
   -O ${output_directory}/${output_filename} \\
   -L ${target_panel} 

EOF
	close(OUT);
if ($cluster_scheduler eq "LSF") {
  system("bsub < $sample_mutect_script");
} elsif ($cluster_scheduler eq "PBS") {
  system("qsub $sample_mutect_script");
}

