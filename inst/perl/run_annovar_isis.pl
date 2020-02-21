#!/usr/bin/perl -w

### AUTHORS ###################################################################
#
# NAME			DATE
# Erle			2017-08-31

### DESCRIPTION ###############################################################
# Run Annovar annotation on a VCF file
#
# TO DO:
#	- decide how to parameterize protocol
#	- make sure log_directory and code_directory end in a slash

### LIBRARIES #################################################################
use strict;
use warnings;
use Getopt::Long;
use YAML qw(LoadFile);
use File::Basename;
use Data::Dumper;

### COMMAND LINE VARIABLES ####################################################
#
# project_directory: path to directory where all project files should be stored. 
#		- each sample will have its own subdirectory here
#		- will create subdirectories log and code (for log files and bash scripts, respectively)
# sample_id: id of sample to be aligned
# config_file: path to config file
# af_threshold: minimum allele frequency
# job_dependencies: jobs the script depends on

GetOptions(
	# required
	"vcf_file=s" => \my $vcf_file,
	"config_file=s" => \my $config_file,
	# optional 
	"output_directory=s" => \my $output_directory,
	"output_filename=s" => \my $output_filename,
	"log_directory=s" => \my $log_directory,
	"code_directory=s" => \my $code_directory,
	"job_dependencies=s{1,}" => \my @job_dependencies,
	"job_name=s" => \my $job_name, 
	"job_group=s" => \my $job_group,
	);

die "$0 requires the vcf_file argument\n" unless $vcf_file;
die "$0 requires the config_file argument" unless $config_file;

# TO DO: fix this
# warn if input file does not end in .vcf – likely to be an error
# my ($vcf_extension) = $vcf_file =~ /((\.[^.\s]+)+)$/;
# warn "Input VCF does not have extension .vcf" unless ( ".vcf" eq $vcf_extension );

# get name of VCF file without direcory
my $vcf_name = basename($vcf_file);

# if no output_directory, set to directory containing tumour BAM file
if( not $output_directory ) {
	$output_directory = dirname($vcf_file);
}

if( not $code_directory ) { 
	$code_directory = $output_directory;
}

if( not $log_directory ) {
	$log_directory = $output_directory;
}

# NOTE: default settings of output_filename are sorted out after the config file has been parsed

### MAIN ######################################################################

my $config = LoadFile($config_file);

my $reference_build = $config->{reference_build};
my $reference_genome = $config->{reference_genome}->{ $reference_build };

# get path to corresponding .dict file
(my $reference_dict = $reference_genome) =~ s{\.[^.]+$}{};
$reference_dict = $reference_dict . ".dict";

# software tools
my $annovar_path = $config->{annovar_path};
my $annovar_database = $config->{annovar_database}->{ $reference_build };

my $gatk_jar = $config->{gatk_jar};
my $picard_jar = $config->{picard_jar};

# module load version numbers
my $java = "java/$config->{java_version}";


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

# cluster settings
my $num_cpu = $config->{num_cpu};
my $project_code = $config->{project_code};
my $cluster_scheduler = $config->{cluster_scheduler};
my $cluster_settings = "";

## ANNOVAR SETTINGS

my @all_annovar_fields = @{ $config->{annovar_fields} };
my @annovar_genotype_fields = @{ $config->{annovar_genotype_fields} };

# add annotation database to annovar fields
# database depends on reference build and annotation method (RefSeq vs GENCODE), so it's easier to add it here than to 
# force user to specify it
my @parsed_annovar_fields = map { (my $field = $_) =~ s/ANNOTATION_DB/$annotation_db/g; $field } @all_annovar_fields;

# if QUAL and STATUS not present in VCF header, drop from parsed fields list
my $field;

if( `grep QUAL $vcf_file | wc -l` == 0 ) {
	@parsed_annovar_fields = grep { !/QUAL/ } @parsed_annovar_fields;
}

if( `grep STATUS $vcf_file | wc -l` == 0 ) {
	@parsed_annovar_fields = grep { !/STATUS/ } @parsed_annovar_fields;
}

# Want each variable with a -F or -GF tag before it. 
# 	=> Add first one manually and use join for the rest.
my $annovar_fields_string = "-F " . join( " -F ", @parsed_annovar_fields);
my $annovar_genotype_string = "-GF " . join( " -GF ", @annovar_genotype_fields);

# if no variants survive the filters, GATK's VariantsToTable will fail. We treat this case separately.
#	=> make a string variable containing only header of VariantsToTable output that can be saved to file
my @header_fields = (@parsed_annovar_fields, @annovar_genotype_fields); 
my $header_string = join( "\t", @header_fields ); 

my $annovar_protocol = " -protocol $annotation_db,cytoBand,cosmic70,1000g2015aug_all,exac03nontcga,clinvar_20170130,nci60,icgc21,dbnsfp30a,dbnsfp31a_interpro -operation g,r,f,f,f,f,f,f,f,f -nastring .";

## JOB DEPENDENCIES

my $job_dependency_string = '';
my $orphan_flag = ''; # if job_dependencies is non-empty, set this to -ti to automatically kill orphan jobs

if( scalar(@job_dependencies) > 0 ) {
	my @parsed_job_dependencies = map { "done(\"$_\")" } @job_dependencies;

	if( scalar(@job_dependencies) > 0 ) {
    if ( $cluster_scheduler eq "LSF" ) {
    	my @parsed_job_dependencies = map { "done(\"$_\")" } @job_dependencies;
    
    	$job_dependency_string = join( " && ", @parsed_job_dependencies );
    	$orphan_flag = '-ti';
    } elsif ( $cluster_scheduler eq "PBS" ) {
      $job_dependency_string = join( ",", @job_dependencies );
    }
  }

}

## JOB GROUP
my $job_group_string = '';
if( $job_group ) { 
	$job_group_string = "#BSUB -g /${job_group}";
}

## OUTPUT FILENAME

if( not $output_filename ) {
	$output_filename = "${vcf_name}.annovar.${buildver}_multianno.vcf.txt";
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

## MAKE SCRIPT

if( not $job_name ) { 
	$job_name = "annotate_${vcf_name}";
}


my $annotation_script = "${code_directory}/${job_name}.sh";

open(OUT,">$annotation_script") or die "Cannot open output script $annotation_script $!\n";
print OUT <<EOF;
$cluster_settings

module load $java;

# STEP 1: run ANNOVAR
#	- input: VCF file with variant calls per sample
#	- output: annotated VCF file
#		* SAMPLE_ID.vcf.annovar.multianno.vcf
perl $annovar_path/table_annovar.pl ${vcf_file} \\
$annovar_protocol $annovar_database \\
-buildver $buildver -vcfinput  \\
-out ${output_directory}/${vcf_name}.annovar


# STEP 2: make tab-delimited table with variants
# 	- input: 
#	- output: tab delimited variant file
#		* SAMPLE_ID.vcf.annovar.${buildver}_multianno.vcf.txt

# Note: VariantsToTable fails if the preceeding VCF file is empty.
n_variants=`grep -v ^# ${output_directory}/${vcf_name}.annovar.${buildver}_multianno.vcf | wc -l`

if [ \$n_variants -gt 0 ]; then
	java -jar $gatk_jar \\
	-T VariantsToTable \\
	-V ${output_directory}/${vcf_name}.annovar.${buildver}_multianno.vcf \\
	-R $reference_genome \\
	-o ${output_directory}/${output_filename} \\
	${annovar_fields_string} \\
	-F DP -GF GT -GF VF;
else
	echo $header_string > ${output_directory}/${output_filename};
fi

EOF
close(OUT);
if ($cluster_scheduler eq "LSF") {
  system("bsub < $annotation_script");
} elsif ($cluster_scheduler eq "PBS") {
  system("qsub $annotation_script");
}
