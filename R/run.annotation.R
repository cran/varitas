#' Run annotation on a set of VCF files
#'
#' @description
#'	Takes a data frame with paths to VCF files, and runs ANNOVAR annotation on each file.
#' 	To allow for smooth connections with downstream pipeline steps, the function returns a variant
#'  specification data frame that can be used as input to merging steps.
#' 
#' @param vcf.specification
#' 	Data frame detailing VCF files to be processed, from \code{prepare.vcf.specification}.
#' @param output.directory 
#' 	Path to folder where code and log files should be stored in their respective subdirectories. 
#' 	If not supplied, code and log files will be stored in the directory with each VCF file.
#' @param job.name.prefix 
#' 	Prefix to be added before VCF name in job name. Defaults to 'annotate', but should be changed if 
#' 	running multiple callers to avoid 
#' @param job.group
#'  Group job should be associated with on cluster
#' @param quiet 
#'  Logical indicating whether to print commands to screen rather than submit them
#' @param verify.options 
#'  Logical indicating whether to run verify.varitas.options
#'
#' @return Data frame with details of variant files
#'
#' @examples
#' run.annotation(
#'   data.frame(
#'     sample.id = c('a', 'b'),
#'     vcf = c('a.vcf', 'b.vcf'),
#'     caller = c('mutect', 'mutect')
#'   ),
#'   output.directory = '.',
#'   quiet = TRUE
#' )
#'
#' @export
run.annotation <- function(
    vcf.specification, 
    output.directory = NULL, 
    job.name.prefix = NULL, 
    job.group = NULL,
    quiet = FALSE,
    verify.options = !quiet
    ) { 
    
    ### INPUT TESTS ###########################################################
    
    if( !is.null(output.directory) && !dir.exists(output.directory) && !quiet ) { 
        error.message <- paste('Directory', output.directory, 'does not exist or is not a directory');
        stop(error.message);
    }
    
    if( !is.null(job.name.prefix) && !is.character(job.name.prefix) ) { 
        stop('job.name.prefix must be a string');
    } 
    
    if( !is.null(job.name.prefix) && length(job.name.prefix) != 1 && length(job.name.prefix) != nrow(vcf.specification) ) {
        stop('job.name.prefix must have length 1 or length equal to vcf.specification');
    }
    
    ### MAIN ##################################################################
    
    # make sure to avoid factor errors
    vcf.specification$sample.id <- as.character(vcf.specification$sample.id);
    vcf.specification$vcf <- as.character(vcf.specification$vcf);
    
    if('caller' %in% names(vcf.specification) ) {
        vcf.specification$caller <- as.character( vcf.specification$caller );
    }

    verify.vcf.specification(vcf.specification);
    
    if( verify.options) {
        verify.varitas.options( stages.to.run = 'annotation' );
    }
    
    vcf.specification$job.name.prefix <- job.name.prefix;
    
    # create directories for log files and code
    # 	- should this be parameterized?
    
    code.directory <- NULL; 
    log.directory <- NULL;
    
    if( !is.null(output.directory) ) { 
        if( !quiet ) {
            create.directories(
                directory.names = c('log', 'code'), 
                path = output.directory
            );
        }
        
        code.directory <- file.path(output.directory, 'code');
        log.directory <- file.path(output.directory, 'log');
        
    }
    
    
    # want to store information about file paths and job dependencies to use downstream 
    variant.specification <- list();
    
    # store paths to output files, and their job dependencies
    for(i in 1:nrow(vcf.specification) ) { 
        
        # TO DO: should VCF specification contain sample ID at all? 
        sample.id <- vcf.specification$sample.id[i];
        vcf.file <- vcf.specification$vcf[i];
        
        # save to same directory as VCF
        annotation.output.directory <- dirname(vcf.file);
        vcf.filename <- basename(vcf.file);
        
        # get name of VCF file – specifying this here allows us to control it and output to 
        buildver <- get.buildver();
        annotation.filename <- paste0(vcf.filename, '.annovar.', buildver, '_multianno.vcf.txt');
        
        # give job a name

        if( 'caller' %in% names(vcf.specification) ) {
            job.name <- paste(vcf.specification$caller[i], vcf.filename, sep = '_');

            if( 'isis' == as.character(vcf.specification$caller[i]) ) {
                isis <- TRUE;
            } else {
                isis <- FALSE;
            }

        } else {
            job.name <- vcf.filename;
            isis <- FALSE;
        }
        
        # sort out job dependencies
        # (assume to be comma or semicolon separated – verify this in other scripts!)
        job.dependencies <- NULL;
        
        if( 'job.dependency' %in% names(vcf.specification) && '' != vcf.specification$job.dependency[i] && !is.na(vcf.specification$job.dependency[i]) ) { 
            job.dependencies <- stringr::str_split(
                vcf.specification$job.dependency[i], 
                pattern = ',(\\s)?|;(\\s)?|\\s'
            )[[1]];
        }

        # add prefix to job name
        if( !is.null(job.name.prefix) && '' != job.name.prefix ) {
            job.name <- paste(job.name.prefix, job.name, sep = '_');
        }

        print(job.name);
        
        run.annovar.vcf(
            vcf.file = vcf.file, 
            output.filename = annotation.filename,
            output.directory = annotation.output.directory,
            code.directory = code.directory, 
            log.directory = log.directory, 
            job.name = job.name,
            job.dependencies = job.dependencies, 
            job.group = 'annovar',
            verify.options = FALSE,
            isis = isis,
            quiet = quiet
            );
        
        sample.variant.specification <- data.frame(
            sample.id = sample.id, 
            variant.file = file.path(annotation.output.directory, annotation.filename),
            job.dependency = job.name, 
            stringsAsFactors = FALSE
        );
        
        if( 'caller' %in% names(vcf.specification) ) { 
            sample.variant.specification$caller <- vcf.specification$caller[i];
        }	
        
        variant.specification[[ i ]] <- sample.variant.specification;
        
    }
    
    # make data frame of variant specifications
    variant.specification <- do.call(rbind, variant.specification);
    
    return(variant.specification);
}
