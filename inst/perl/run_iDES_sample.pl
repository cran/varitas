#!/usr/bin/perl -w

### AUTHORS ###################################################################
# Name			Date
# Erle			2017-08-21

### DESCRIPTION ###############################################################
#
# Run iDES on a sample in a project directory
# Creates an output directory named ides under the sample directory.
#
# TO DO: 
#	- add proton option back in
#	- move to shared Perl module library
#
# QUESTIONS:
# 	- why run this on BAM with all reads rather than ontarget ones?
#	- why run this after replacing read groups?
#
# NOTES: 
#	- Ros had originally implemented this as an array job. I opted to do a single 
#	  sample thing to keep with the conventions from VarDict and alignment
#	- iDES step 1 requires the extension sorted.bam for all input files
#	

### LIBRARIES #################################################################
use strict;
use Getopt::Long;
use YAML qw(LoadFile);

### COMMAND LINE ARGUMENTS ####################################################
#
# project_directory: path to directory where all project files should be stored. 
#		- each sample will have its own subdirectory here
#		- will create subdirectories log and code (for log files and bash scripts, respectively)
# sample_id: id of sample to be aligned
# config_file: path to config file
# job_dependencies: prefix of jobs the script depends on (will add _SAMPLE_ID)

GetOptions(
	"project_directory=s" => \ my $project_directory, 
	"sample_id=s" => \ my $sample_id,
	"config_file=s" => \(my $config_file = 'config.yaml'),
	"job_dependencies=s{1,}" => \my @job_dependencies,
	);

die "$0 requires the sample id argument (--sample_id)\n" unless $sample_id;
die "$0 requires the project directory argument (--project_directory)\n" unless $project_directory;

# idea: pass "frequently changing parameters" as command line
# include config file for more complicated 
my $config = LoadFile($config_file);


my $reference_build = $config->{reference_build};
my $reference_genome = $config->{reference_genome}->{ $reference_build };
my $target_panel = $config->{target_panel};

# cluster settings
my $project_code = $config->{project_code};

# software tools
my $annovar_path = $config->{annovar_path};
my $annovar_database = $config->{annovar_database}->{ $reference_build };

my $ides_path = $config->{ides_path};


# set buildver parameter based on config file
my $buildver = "hg19";
my $gencode_version = "V19";

if( "grch38" eq lc $reference_build ) { 
	$buildver = "hg38";
	$gencode_version = "V24";
}

# adjust ANNOVAR calls based on desired gene annotation method
my $annotation_db = "refGene";

if( "gencode" eq lc ${config}->{annotation_method} ) { 
	$annotation_db = "wgEncodeGencodeBasic" . $gencode_version;
}


### MAIN ######################################################################

# parse job dependencies into string that can be used in job preamble
my $job_dependency_string = '';
my $orphan_flag = ''; # if job_dependencies is non-empty, set this to -ti to automatically kill orphan jobs

if( scalar(@job_dependencies) > 0 ) {
	my @parsed_job_dependencies = map { "done(\"$_\_$sample_id\")" } @job_dependencies;

	$job_dependency_string = join(" && ", @parsed_job_dependencies );
	$orphan_flag = '-ti';

}

## QSUB/BSUB SETTINGS
if($cluster_scheduler eq "LSF") {
  $cluster_settings = <<END;
#!/bin/bash
#BSUB -J ${job_name}
#BSUB -o ${log_directory}/${job_name}.out
#BSUB -e ${log_directory}/${job_name}.err 
#BSUB -n 1
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
#PBS -n 1
#PBS -q normal
#PBS -l cput=0:168:00
${job_group_string}
END
} else {
  $cluster_settings = ""
}

# module load version numbers
my $java = "java/$config->{java_version}";
my $samtools = "samtools/$config->{samtools_version}";
my $r = "R/$config->{r_version}";


my $annovar_protocol = " -protocol $annotation_db,cytoBand,cosmic70,1000g2015aug_all,exac03nontcga,clinvar_20170130,nci60,icgc21,dbnsfp30a -operation g,r,f,f,f,f,f,f,f -nastring .";

my $sample_directory = "$project_directory/$sample_id";
my $sample_ides_script = "$project_directory/code/ides_$sample_id.sh";

my $ides_step1 = "$ides_path/ides-bam2freq.pl";

open(OUT,">$sample_ides_script") or die "Cannot open output script $sample_ides_script $!\n";
print OUT <<EOF;
$cluster_settings

# Only need Java due to the RJava issue on Davros
module load $samtools $java $r;

cd $sample_directory
mkdir ides

# workaround until we get shared Perl library
export PERL5LIB=/scratch/DBC/MONCOLOGY/rcutts/apps/iDES_software_NBT

# STEP 1: run iDES pileup command 
#	- input: BAM file
#	- output: file giving strand-specific calls for each base
# 		* SAMPLE_ID.sorted.freq.allreads.Q30.txt
$ides_step1 -a -o ${sample_directory}/ides/ \\
${sample_directory}/${sample_id}.sorted.bam \\
$reference_genome $target_panel

# STEP 2: convert to ANNOVAR-friendly format
# 	- input: strand specific calls
#	- output: table with calls per alt base and allele frequency
#		* SAMPLE_ID.sorted.freq.allreads.Q30.txt.calls.txt
R -q -e "library($config->{package_name}); convert.ides.output('${sample_directory}/ides/${sample_id}.sorted.freq.allreads.Q30.txt');"

# STEP 3: variant annotation
#	- input: 
#	- output: 
#		* SAMPLE_ID.sorted.freq.allreads.Q30.txt.annovar
perl $annovar_path/table_annovar.pl ${sample_directory}/ides/${sample_id}.sorted.freq.allreads.Q30.txt.calls.txt \\
$annovar_protocol $annovar_database \\
-buildver $buildver  \\
-out ${sample_directory}/ides/${sample_id}.sorted.freq.allreads.Q30.txt.calls.txt.annovar

# STEP 4: merge iDES (potential) variant calls with annotations
# 	- input: iDES variant calls and Annovar annotation
#	- output: merged table with potential calls and annotations
#		* SAMPLE_ID.sorted.freq.allreads.Q30.txt.calls.txt.ann.txt
R -q -e " \
library($config->{package_name}); 
merge.ides.annotation('${sample_directory}/ides/${sample_id}.sorted.freq.allreads.Q30.txt.calls.txt'); 
"

EOF
close(OUT);
   system("bsub < $sample_ides_script");

