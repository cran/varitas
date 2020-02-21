#' Perform sample QC by looking at target coverage.
#'
#' @inheritParams run.variant.calling
#' @param project.directory
#'  Path to project directory where code and log files should be saved
#' @param paired 
#'  Logical indicating whether the analysis is paired. This does not affect QC directly, but means normal samples get nested
#' @param output.subdirectory
#'  If further nesting is required, name of subdirectory. If no further nesting, set to FALSE
#' @param job.group
#'  Group job should be associated with on cluster
#'
#'
run.target.qc <- function(
    bam.specification,
    project.directory,
    sample.directories = TRUE,
    paired = FALSE,
    output.subdirectory = FALSE,
    quiet = FALSE,
    job.name.prefix = NULL,
    verify.options = FALSE,
    job.group = 'target_qc'
    ) { 
    
    ### INPUT TESTS ###########################################################
    
    verify.bam.specification(bam.specification);
    
    ### MAIN ##################################################################
    
    bam.specification$sample.id <- as.character(bam.specification$sample.id);
    bam.specification$tumour.bam <- as.character(bam.specification$tumour.bam);

    if( 'normal.bam' %in% names(bam.specification) ) {
        bam.specification$normal.bam <- as.character(bam.specification$normal.bam);
    }

    if( verify.options ) {
        # TO DO!
    }
    
    # create directories for log files and code
    # 	- should this be parameterized?
    if( !quiet ) {
        create.directories(
            directory.names = c('log', 'code'), 
            path = project.directory
        );
    }
    
    # add a column to bam.specification data frame with path to 
    # output direcory for that specific sample
    if( sample.directories ) { 
        
        bam.specification$output.path <- file.path(
            project.directory,
            bam.specification$sample.id
        );
        
    } else { 
        bam.specification$output.path <- project.directory;
    }
    
    # if further nesting has been requested, add to path
    if( !identical(output.subdirectory, FALSE) ) { 
        bam.specification$output.path <- file.path(
            bam.specification$output.path,
            output.subdirectory
        );
    }
    
    job.depends <- c()
    
    for( i in 1:nrow(bam.specification) ) { 
        sample.id <- bam.specification$sample.id[i];
        tumour.bam <- bam.specification$tumour.bam[i];
        sample.output.directory <- bam.specification$output.path[i];

        job.dependencies <- NULL;
        if( 'job.dependency' %in% names(bam.specification) && '' != bam.specification$job.dependency[i] && !is.na(bam.specification$job.dependency[i]) ) {
            job.dependencies <- stringr::str_split(
                bam.specification$job.dependency[i], 
                pattern = '\\s+'
            );
        }

        job.name <- paste0('target_qc_', sample.id);
        if( !is.null(job.name.prefix) && '' != job.name.prefix ) {
            job.name <- paste(job.name.prefix, job.name, sep = '_');
        }
        
    
        run.target.qc.sample(
            bam.file = tumour.bam,
            sample.id = sample.id,
            output.directory = sample.output.directory,
            job.dependencies = job.dependencies,
            code.directory = file.path(project.directory, 'code'),
            log.directory = file.path(project.directory, 'log'),
            job.group = 'target_qc',
            job.name = job.name,
            quiet = quiet
        );
        bam.specification$tumour.bam[i] <- file.path(sample.output.directory, paste0(sample.id, '.sorted.bam.ontarget.bam'))
        job.dependency <- job.name
        
        if( paired ) { 
            
            normal.bam <- bam.specification$normal.bam[i];
            
            normal.job.name <- paste0('target_qc_', sample.id, '-NORMAL');
            if( !is.null(job.name.prefix) && '' != job.name.prefix ) {
              normal.job.name <- paste(job.name.prefix, normal.job.name, sep = '_');
            }
            
            run.target.qc.sample(
                bam.file = normal.bam,
                sample.id = paste0(sample.id, '-NORMAL'),
                output.directory = sample.output.directory,
                code.directory = file.path(project.directory, 'code'),
                log.directory = file.path(project.directory, 'log'),
                job.name = normal.job.name,
                job.group = 'target_qc',
                quiet = quiet
            );	
            bam.specification$normal.bam[i] <- file.path(sample.output.directory, paste0(sample.id, '-NORMAL', '.sorted.bam.ontarget.bam'))
            job.dependency <- paste(job.dependency, normal.job.name)
        }
        job.depends <- c(job.depends, job.dependency)
        
    }
    bam.specification['job.dependency'] <- job.depends
    
    return(bam.specification)
}