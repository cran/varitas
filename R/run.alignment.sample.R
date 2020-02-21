#' Run alignment for a single sample
#'
#' @param fastq.files Paths to FASTQ files (one file if single-end reads, two files if paired-end)
#' @inheritParams run.vardict.sample
#' 
run.alignment.sample <- function(
    fastq.files,
    sample.id,
    output.directory = NULL,
    output.filename = NULL,
    code.directory = NULL,
    log.directory = NULL,
    config.file = NULL,
    job.dependencies = NULL,
    job.name = NULL,
    job.group = NULL,
    quiet = FALSE,
    verify.options = !quiet
    ) {
    
    ### INPUT TESTS ###########################################################
    
    # FASTQ files
    assertthat::assert_that(
        length(fastq.files) <= 2,
        msg = 'Cannot accept more than two FASTQ files'
    );
    
    assertthat::assert_that(
        0 != length(fastq.files), 
        msg = 'Must supply at least one FASTQ file'
    );
    
    assertthat::assert_that(
        quiet || !is.null(job.dependencies) || all( file.exists(fastq.files) ), 
        msg = 'No job dependency supplied, yet FASTQ files do not exist'
    );
    
    
    # sample ID
    assertthat::assert_that( 
        1 == length(sample.id), 
        msg = 'sample.id must have length 1'
    );
    
    
    # could possibly allow numeric.. but be strict to start with
    if( is.factor(sample.id) ) sample.id <- as.character(sample.id);
    assertthat::assert_that( is.character(sample.id) );
    
    ### MAIN ##################################################################
    
    if( verify.options ) {
        verify.varitas.options(
            varitas.options = config.file, 
            stages.to.run = 'alignment'
        );
    }
    
    
    # if no config.file has been supplied, save a temporary one and pass to Perl 	
    if( is.null(config.file) ) {
        config.file <- save.config();
    }
    
    script <- system.file('perl', 'run_alignment.pl', package = get.varitas.options('pkgname') );

    
    # Might want to rename the ontarget_bam_filename option...
    command <- make.command.line.call(
        main.command = c('perl', script), 
        options = list(
            'fastq' = fastq.files,
            'sample_id' = sample.id,
            'config_file' = config.file,
            'output_directory' = output.directory,
            'ontarget_bam_filename' = output.filename,
            'code_directory' = code.directory,
            'log_directory' = log.directory,
            'job_dependencies' = job.dependencies,
            'job_name' = job.name,
            'job_group' = job.group
        )
    );
    
    if(quiet) { 
        cat(command, '\n');
    } else { 
        system(command);
    }
    
}