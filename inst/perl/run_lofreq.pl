#!/usr/bin/perl -w

### AUTHORS ###################################################################
# Name			Date
# Adam			2019-02-07

### DESCRIPTION ###############################################################
#
# Run LoFreq on paired tumour normal samples (or non paired)
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
	$output_filename = "${sample_id}.vcf";
}

# if no job name provided, set to varscan_SAMPLE_ID.vcf
if( not $job_name ) {
	my $job_name = "lofreq_${sample_id}";
}

### MAIN ######################################################################

my $config = LoadFile($config_file);
my $reference_build = $config->{ reference_build };
my $reference_genome = $config->{ reference_genome }->{ $reference_build };
my $lofreq = $config->{lofreq_path};
my $dbsnp = $config->{germline_vcf};
my $bed_file = $config->{target_panel};
my $bam_dir = dirname($tumour_bam);
my $platform = $config->{sequencing_platform};

# Module load versions
my $samtools = "samtools/$config->{samtools_version}";

# cluster settings
my $num_cpu = $config->{num_cpu};
my $project_code = $config->{project_code};
my $cluster_scheduler = $config->{cluster_scheduler};
my $cluster_settings = "";

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

## Enable lofreq viterbi only for illumina data
if( $platform eq "illumina" ) {
  if( $paired ) {
    my $viterbi = "${lofreq} viterbi -f ${reference_genome} ${normal_bam} | samtools sort - -@ 5 -o ${sample_id}_normal_ready.bam; ${lofreq} viterbi -f ${reference_genome} ${tumour_bam} | samtools sort - -@ 5 -o ${sample_id}_tumour_ready.bam"
  } else {
    my $viterbi = "${lofreq} viterbi -f ${reference_genome} ${tumour_bam} | samtools sort - -@ 5 -o ${sample_id}_tumour_ready.bam"
  }
} else {
  if( $paired ) {
    my $viterbi = "cp ${normal_bam} ${sample_id}_normal_ready.bam; cp ${tumour_bam} ${sample_id}_tumour_ready.bam"
  } else {
    my $viterbi = "cp ${tumour_bam} ${sample_id}_tumour_ready.bam"
  }
}


#SCRIPT
my $sample_lofreq_script="${code_directory}/${job_name}.sh";

open(OUT,">$sample_lofreq_script") or die "Cannot open output script $!\n";

if( $paired ) {
  print OUT <<EOF;
${cluster_settings}

module load ${samtools}

mkdir -p ${output_directory} ${log_directory} ${code_directory}

cd ${bam_dir}

${viterbi}

${lofreq} indelqual --dindel -f ${reference_genome} --out ${sample_id}_normal_indel.bam ${sample_id}_normal_ready.bam
${lofreq} indelqual --dindel -f ${reference_genome} --out ${sample_id}_tumour_indel.bam ${sample_id}_tumour_ready.bam
samtools index ${sample_id}_normal_indel.bam
samtools index ${sample_id}_tumour_indel.bam

# Add baseline indel qualities to BAM files as they are input
cd ${output_directory}
${lofreq} somatic \\
-n ${bam_dir}/${sample_id}_normal_indel.bam \\
-t ${bam_dir}/${sample_id}_tumour_indel.bam \\
-f ${reference_genome} \\
-d ${dbsnp} \\
--threads ${num_cpu} \\
--call-indels \\
-o ${sample_id}_ \\
-l ${bed_file}

zcat ${sample_id}_somatic_final*.vcf.gz > ${output_filename}

EOF

} else {
  print OUT <<EOF;
${cluster_settings}

module load ${samtools}

mkdir -p ${output_directory} ${log_directory} ${code_directory}

cd ${bam_dir}

${viterbi}

${lofreq} indelqual --dindel -f ${reference_genome} --out ${sample_id}_tumour_indel.bam ${sample_id}_tumour_ready.bam
samtools index ${sample_id}_tumour_indel.bam

cd ${output_directory}
${lofreq} call-parallel \\
--force-overwrite \\
-f ${reference_genome} \\
--pp-threads ${num_cpu} \\
--call-indels \\
-d ${dbsnp} \\
-o ${output_filename}_dup \\
-l ${bed_file} \\
${bam_dir}/${sample_id}_tumour_indel.bam

# Remove duplicate variants (happens when regions overlap in BED file)
awk '!a[\$0]++' ${output_filename}_dup > ${output_filename}

EOF

}
	close(OUT);
if ($cluster_scheduler eq "LSF") {
  system("bsub < $sample_lofreq_script");
} elsif ($cluster_scheduler eq "PBS") {
  system("qsub $sample_lofreq_script");
}
