
#!/usr/bin/perl -w

### AUTHORS ###################################################################
# Name			Date
# Erle			2017-08-18

### DESCRIPTION ###############################################################
#
# 
#

### LIBRARIES #################################################################
use strict;
use Getopt::Long;
use YAML qw(LoadFile);
use File::Basename;

### COMMAND LINE VARIABLES ####################################################
#
# variant_file
# 	- path to variant file to be filtered
# config_file
# 	- path to config file
# output_directory
#	- path to directory where results should be saved (defaults to directory with tumour_bam)
# output_filename
#	- Filename of final output file (defaults to ${sample_id}.vcf)
# code_directory 
#	- directory where code should be saved (defaults to output_directory)
# log_directory
#	- directory where logs should be saved (defaults to output_directory)
# job_dependencies 
#	- jobs the script depends on 
# job_group 
#	- group job belongs to

GetOptions(
	# Required
	"variant_file=s" => \my $variant_file,
	"variant_caller=s" => \my $variant_caller,
	"config_file=s" => \my $config_file,
	# Optional
	"output_directory=s" => \my $output_directory, 
	"output_filename=s" => \my $output_filename,  
	"log_directory=s" => \my $log_directory, 
	"code_directory=s" => \my $code_directory,
	"job_dependencies=s{1,}" => \my @job_dependencies,
	"job_group=s" => \my $job_group,
	);

die "$0 requires the variant_file argument\n" unless $variant_file;
die "$0 requires the config_file argument\n" unless $config_file;


# if no output_directory, set to directory containing tumour BAM file
if( not $output_directory ) {
	$output_directory = dirname($variant_file);
}

if (not $output_filename ) { 
	$output_filename = basename($variant_file) . ".filtered.txt";
}

if( not $code_directory ) { 
	$code_directory = $output_directory;
}

if( not $log_directory ) {
	$log_directory = $output_directory;
}

### MAIN ######################################################################

my $config = LoadFile($config_file);


## JOB DEPENDENCIES

my $job_dependency_string = '';
my $orphan_flag = ''; # if job_dependencies is non-empty, set this to -ti to automatically kill orphan jobs

if( scalar(@job_dependencies) > 0 ) {
	my @parsed_job_dependencies = map { "done(\"$_\")" } @job_dependencies;

	$job_dependency_string = join(" && ", @parsed_job_dependencies );
	$orphan_flag = '-ti';
}

## JOB GROUP
my $job_group_string = '';
if( $job_group ) { 
	$job_group_string = "#BSUB -g /${job_group}";
}

my $job_name = "filter_" . basename($variant_file);
my $sample_filtering_script = "${code_directory}/${job_name}.sh";


open(OUT,">$sample_filtering_script") or die "Cannot open output script $sample_filtering_script $!\n";

print OUT <<EOF;

#!/bin/bash
#BSUB -J ${job_name}
#BSUB -o ${log_directory}/${job_name}.out
#BSUB -e ${log_directory}/${job_name}.err 
#BSUB -n 1
#BSUB -q normal
#BSUB -P $config->{project_code}
#BSUB -w '${job_dependency_string}' ${orphan_flag}
#BSUB -W 124:00
${job_group_string}

mkdir -p ${output_directory} ${log_directory} ${code_directory}

R -q -e " \\
library($config->{pkgname}); \\
filter.variant.file( \\
 variant.file = '${variant_file}', \\
 caller = '${variant_caller}', \\
 config.file = '${config_file}', \\
 output.file = '${output_directory}/${output_filename}' \\
 );
"

EOF

close(OUT);
   system("bsub < $sample_filtering_script");
