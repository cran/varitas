#' Check that FASTQ specification data frame matches expected format, and that all files exist
#'
#' @param fastq.specification Data frame containing columns sample.id and reads, and optionally a column mates
#' @param paired.end Logical indicating whether paired end reads are used
#' @param files.ready Logical indicating if the files already exist on disk. If there are job dependencies, this should be set to FALSE.
#'
#' @return None

#'
#'
verify.fastq.specification <- function(
    fastq.specification, 
    paired.end = FALSE, 
    files.ready = FALSE
    ) {
    
    ## check type
    if( !is.data.frame(fastq.specification) ) { 
        stop('fastq.specification is not a data frame.');
    }
    
    ## check data frame dimensions – expect either 2 or 3
    if(ncol(fastq.specification) < 2) { 
        stop('fastq.specification has fewer than 2 columns.');
    }
    
    ## check column headers
    for(label in c('sample.id', 'reads')) { 
        if( !(label %in% names(fastq.specification)) ) { 
            error.message <- paste('No column named', label, 'in fastq.specification');
            stop(error.message);
        }	
    }
    
    if( paired.end && !( 'mates' %in% names(fastq.specification) ) ) { 
        stop('No column named mates in fastq.specification');
    }	
    
    # if files are supposed to be ready on disk,
    # check that files exist
    if( files.ready ) {
        read.files <- fastq.specification$reads;
        
        if(paired.end) {
            read.files <- c(read.files, fastq.specification$mates);
        }
        
        reads.exist <- file.exists(read.files);
        
        if( !all(reads.exist) ) { 
            error.message <- paste0(
                'The following FASTQ files do not exist:\n', 
                paste(read.files[!reads.exist], collapse = '\n')
            );
            
            stop(error.message);
            
        }
        
        # check for white space in sample IDs
        if( any(grepl('\\s', fastq.specification$sample.id)) ) { 
            stop('Sample IDs can not contain whitespace.');
        }
        
    }
    
}