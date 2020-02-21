
### DESCRIPTION ###############################################################
#
# Submit post-processing job to the cluster with approporiate job dependencies
#
# TO DO:
#	- make perl module for parsing job dependencies
#
#

### LIBRARIES #################################################################
use strict;
use Getopt::Long;
use YAML qw(LoadFile);
use File::Basename;


### COMMAND LINE VARIABLES ####################################################
# variant_specification: path to data frame specifying which variants to process
# config_file: path to config file
# 
# Optional:
# 	output_directory: directory for output files
# 	code_directory: directory for code
# 	log_directory: directory for log files
#  	job_dependencies: names of jobs that need to have finished before post-processing can start
#	job_name: name of job
#	email: email address that should be notified at the end of job
#
GetOptions(
	# Required
	"variant_specification=s" => \my $variant_specification, 
	"config_file=s" => \my $config_file,
	# Optional
	"output_directory=s" => \my $output_directory,
	"code_directory=s" => \my $code_directory,
	"log_directory=s" => \my $log_directory,
	"job_dependencies=s{1,}" => \my @job_dependencies,
	"job_name=s" => \my $job_name,
	"email=s" => \my $email,
	"quiet" => \my $quiet
	);

die "$0 requires the variant_specification argument" unless $variant_specification;
die "$0 requires the config_file argument" unless $config_file;

if( not $job_name ) {
	$job_name = "post_processing";
}

if( not $output_directory ) {
	$output_directory = dirname($variant_specification);
} 

if( not $code_directory ) {
	$code_directory = $output_directory;
}

if( not $log_directory ) {
	$log_directory = $output_directory;
}

# load config file (needed for project code and package name)
my $config = LoadFile($config_file);

# cluster settings
my $num_cpu = $config->{num_cpu};
my $project_code = $config->{project_code};
my $cluster_scheduler = $config->{cluster_scheduler};
my $cluster_settings = "";

### MAIN ######################################################################

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

# email notification
my $email_string = '';
if( $email && $cluster_scheduler eq "LSF") {
	$email_string = "#BSUB -N -u $email";
} elsif ( $email && $cluster_scheduler eq "PBS" ) {
  $email_string = "#PBS -m e -M $email";
}

## QSUB/BSUB SETTINGS
if( $cluster_scheduler eq "LSF" ) {
  $cluster_settings = <<END;
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -o ${log_directory}/${job_name}.out
#BSUB -e ${log_directory}/${job_name}.err 
#BSUB -n 1
#BSUB -q normal
#BSUB -P $config->{project_code}
#BSUB -w '${job_dependency_string}' ${orphan_flag}
#BSUB -W 168:00
${email_string}
END
} elsif ( $cluster_scheduler eq "PBS" ) {
  $cluster_settings = <<END;
#!/bin/bash
#PBS -P $config->{project_code}
#PBS -N ${job_name}
#PBS -o ${log_directory}/${job_name}.out
#PBS -e ${log_directory}/${job_name}.err 
#PBS -n 5
#PBS -q normal
#PBS -hold_jid $job_dependency_string
#PBS -l cput=0:168:00
${email_string}
END
} else {
  $cluster_settings = ""
}

my $post_processing_script = "${code_directory}/${job_name}.sh";

open(OUT,">$post_processing_script") or die "Cannot open output script $post_processing_script $!\n";

print OUT <<EOF;
${cluster_settings}

module load R/$config->{r_version};

mkdir -p ${output_directory} ${log_directory} ${code_directory}

R -q -e " \\
library($config->{pkgname}); \\
varitas:::post.processing( \\
 variant.specification = '${variant_specification}', \\
 project.directory = '${output_directory}', \\
 config.file = '${config_file}', \\
 sleep = TRUE \\
 );
"

EOF

close(OUT);

if( not $quiet ) {
  if ($cluster_scheduler eq "LSF") {
    system("bsub < $post_processing_script");
  } elsif ($cluster_scheduler eq "PBS") {
    system("qsub $post_processing_script");
  }
}