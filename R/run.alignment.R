#' Run alignment
#'
#' @details 
#' Runs alignment (and related processing steps) on each sample.
#'
#' @param fastq.specification 
#'	Data frame detailing FASTQ files to be processed, typically from prepare.fastq.specification 
#' @param output.directory 
#' 	Path to project directory
#' @param paired.end
#'  Logical indicating whether paired-end sequencing was performed
#' @param sample.directories
#'  Logical indicating whether all sample files should be saved to sample-specific subdirectories (will be created)
#' @param output.subdirectory
#'  If further nesting is required, name of subdirectory. If no further nesting, set to FALSE
#' @param job.name.prefix
#'  Prefix for job names on the cluster
#' @param job.group
#'  Group job should be associated with on cluster
#' @param quiet 
#'  Logical indicating whether to print commands to screen rather than submit them
#' @param verify.options
#'  Logical indicating whether to run verify.varitas.options
#' 
#' @examples
#' run.alignment(
#'   fastq.specification = data.frame(
#'     sample.id = c('1', '2'),
#'     reads = c('1-R1.fastq.gz', '2-R1.fastq.gz'),
#'     mates = c('1-R2.fastq.gz', '2-R2.fastq.gz'),
#'     patient.id = c('P1', 'P1'),
#'     tissue = c('tumour', 'normal')
#'   ),
#'   output.directory = '.',
#'   quiet = TRUE,
#'   paired.end = TRUE
#' )
#'
#' @return None
#'
#' @export
run.alignment <- function(
    fastq.specification,
    output.directory,
    paired.end = FALSE,
    sample.directories = TRUE,
    output.subdirectory = FALSE,
    job.name.prefix = NULL,
    job.group = 'alignment',
    quiet = FALSE,
    verify.options = !quiet
    ) { 
    
    # TO DO:
    #	- make job.dependencies argument more robust to misspecification

    # avoid factor problems
    fastq.specification$sample.id <- as.character(fastq.specification$sample.id);
    fastq.specification$reads <- as.character(fastq.specification$reads);

    if(paired.end) fastq.specification$mates <- as.character(fastq.specification$mates);

    ### INPUT TESTS ###########################################################
    
    # Check that output directory exists. 
    # In theory we could create the directory if it does not exist, but it would have to be created recursively. 
    # To be safe, I have opted to throw an error.
    if( !quiet && !dir.exists(output.directory) ) { 
        error.message <- paste('Directory', output.directory, 'does not exist.');
        stop(error.message);
    }
    
    # if not submitting commands, do not check that files exist
    verify.fastq.specification(
        fastq.specification, 
        paired.end = paired.end, 
        files.ready = !quiet
        );	
    
    # play it safe by refusing to proceed if unclear about whether it's paired-end
    if( 'mates' %in% names(fastq.specification) && !paired.end ) {
        stop('fastq.specification contains a column mates, but paired.end = FALSE');
    }

    ### MAIN ##################################################################
    
    if( verify.options ) {
        verify.varitas.options(stages.to.run = 'alignment');
    }
    
    # create directories for log files and code
    # 	- should this be parameterized?
    if( !quiet ) {
        create.directories(
            directory.names = c('log', 'code'), 
            path = output.directory
            );
    }
    
    # determine output path for each sample
    if(sample.directories) { 
        fastq.specification$output.path <- file.path(
            output.directory, 
            fastq.specification$sample.id
        );
    } else { 
        fastq.specification$output.path <- output.directory;
    }
    
    # if further nesting has been requested, add to path
    if( !identical(output.subdirectory, FALSE) && !is.null(output.subdirectory) ) { 
        fastq.specification$output.path <- file.path(
            fastq.specification$output.path,
            output.subdirectory
        );
    }
    
    
    # Loop over samples and run each one
    bam.specification <- list();
    
    # Special handling for normal samples
    if( 'tissue' %in% names(fastq.specification) ) {
      dict <- c()
      keys <- c()
      for (i in 1:nrow(fastq.specification)) {
        if (fastq.specification$tissue[i] == 'tumour') {
          keys <- c(keys, fastq.specification$patient.id[i])
          dict <- c(dict, fastq.specification$sample.id[i])
        }
      }
      names(dict) <- keys
    }
    
    for( i in 1:nrow(fastq.specification) ) { 
        
        sample.id <- fastq.specification$sample.id[i];
        sample.output.directory <- fastq.specification$output.path[i];
        fastq.files <- fastq.specification$reads[i];
        
        if(paired.end) { 
            fastq.files <- c(fastq.files, fastq.specification$mates[i]);
        }
        
        # should the output pattern be parameterized ?
        output.filename <- paste0(sample.id, '.sorted.bam.ontarget.bam');
        output.file <- file.path(sample.output.directory, output.filename);
        job.name <- paste0('align_', sample.id);

        if( !is.null(job.name.prefix) && '' != job.name.prefix ) {
            job.name <- paste(job.name.prefix, job.name, sep = '_');
        }
        
        if( 'tissue' %in% names(fastq.specification) ) {
           if ( 'normal' == fastq.specification$tissue[i] ) {
             print(paste(sample.id, 'goes with', dict[fastq.specification$patient.id[i]]))
             tumour.sample <- dict[fastq.specification$patient.id[i]]
             sample.id <- paste0(tumour.sample, '-NORMAL')
           }
        }
        
        # should the config file be passed through?
        run.alignment.sample(
            sample.id = sample.id,
            fastq.files = fastq.files,
            output.directory = sample.output.directory,
            output.filename = output.filename,
            code.directory = file.path(output.directory, 'code'),
            log.directory = file.path(output.directory, 'log'),
            job.group = job.group,
            job.name = job.name,
            verify.options = FALSE,
            quiet = quiet
            );
        
        # TO DO: figure out how to deal with tumour/ normal issue
        sample.bam.specification <- data.frame(
            sample.id = sample.id, 
            bam.file = output.file, 
            job.dependency = job.name
            );
        
        # add patient ID and caller to BAM speciication if the information is available
        if( 'patient.id' %in% names(fastq.specification) ) {
            sample.bam.specification$patient.id <- fastq.specification$patient.id[i];
        }
        
        if( 'tissue' %in% names(fastq.specification) ) {
            sample.bam.specification$tissue <- fastq.specification$tissue[i];
        }
        
        bam.specification[[ i ]] <- sample.bam.specification;
        
    }
    
    bam.specification <- do.call(rbind, bam.specification);
    
    # If tumour and normal is present, reformat to tumour.bam and normal.bam format
    # this should probably be moved to a different 
    
    reformatted.bam.specification <- list();
    
    if( 'patient.id' %in% names(bam.specification) && 'tissue' %in% names(bam.specification) ) {
        
        tumour.bams <- bam.specification['tumour' == bam.specification$tissue, ];
        
        for( i in 1:nrow(tumour.bams) ) {

            sample.id <- tumour.bams$sample.id[i];
            patient.id <- tumour.bams$patient.id[i];
            
            patient.normal.samples <- bam.specification['normal' == bam.specification$tissue & patient.id == bam.specification$patient.id, ];
            
            if( 0 == nrow(patient.normal.samples) ) {
                # no normal BAM
                normal.bam <- NA;
                normal.job.dependency <- NA;
                
            } else if ( 1 == nrow(patient.normal.samples) ) {
                normal.bam <- patient.normal.samples$bam.file[1];
                normal.job.dependency <- patient.normal.samples$job.dependency[1];
                
            } else {
                stop('Multi-region normal samples not supported - sorry!');
            }
            
            merged.job.dependency <- tumour.bams$job.dependency[i];
            if( !is.na(normal.job.dependency) ) {
                merged.job.dependency <- paste(merged.job.dependency, normal.job.dependency);
            }
            
            reformatted.bam.specification[[ sample.id ]] <- data.frame(
                sample.id = sample.id, 
                tumour.bam = tumour.bams$bam.file[i], 
                normal.bam = normal.bam, 
                job.dependency = merged.job.dependency
            );
            
        }
        
        bam.specification <- do.call(rbind, reformatted.bam.specification);
        
    } else {
        # no need to  
        names(bam.specification)[ 'bam.file' == names(bam.specification) ] <- 'tumour.bam';
    }
    
    return(bam.specification);
    
}