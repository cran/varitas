#' Check that sample specification data frame matches expected format, and that all files exist
#'
#' @param bam.specification Data frame containing columns sample.id and tumour.bam, and optionally a column normal.bam.
#'
#' @return None
#'
#'
verify.bam.specification <- function(bam.specification) { 
    
    ## check type
    if( !is.data.frame(bam.specification) ) { 
        stop('bam.specification is not a data frame.');
    }
    
    ## check data frame dimensions – expect either 2 or 3
    if(ncol(bam.specification) < 2) { 
        stop('bam.specification has fewer than 2 columns.');
    }
    
    if(ncol(bam.specification) > 4) { 
        stop('bam.specification has more than 4 columns.');
    }	
    
    ## check column headers
    for(label in c('sample.id', 'tumour.bam')) { 
        if( !(label %in% names(bam.specification)) ) { 
            error.message <- paste('No column named', label, 'in bam.specification');
            stop(error.message);
        }	
    }
    
    # can't check that files exist.. since they might be created later
    
    # check for white space in sample IDs
    if( any( grepl('\\s', bam.specification$sample.id) ) ) { 
        stop('Sample IDs can not contain whitespace.');
    }
    
}
