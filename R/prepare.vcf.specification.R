#' prepare.vcf.specification 
#'
#' @description
#'  Prepare VCF specification data frame for annotation
#'
#' @param vcf.details 
#'  Data frame containing details of VCF files
#' @param sample.id.column 
#'  Identifier of column in \code{vcf.details} containing sample IDs (index or name)
#' @param vcf.column 
#'  Identifier of column in \code{vcf.details} containing VCF file (index or name)
#' @param job.dependency.column 
#'  Identifier of column in \code{vcf.details} containing job dependency (index or name)
#' @param caller.column
#'  Identifier of column in \code{vcf.details} containing caller (index or name)
#'
#' @return Properly formatted VCF details
#'
#'
prepare.vcf.specification <- function(
    vcf.details, 
    sample.id.column = 1, 
    vcf.column = 2, 
    job.dependency.column = NA, 
    caller.column = NA
    ) { 
    
    ### INPUT TESTS ###########################################################
    
    if( !is.data.frame(vcf.details) ) {
        stop('vcf.details must be a data frame');
    }
    
    required.columns <- c(
        'sample.id.column' = sample.id.column, 
        'vcf.column' = vcf.column
    );
    
    for( column in names(required.columns) ) { 
        
        column.identifier <- required.columns[column]; 
        error.message <- paste(column, 'not found in vcf.details data frame');
        
        if( is.numeric(column.identifier) && column.identifier > ncol(vcf.details) ) { 
            stop(error.message);
        }	
        
        if( is.character(column.identifier) && !(column.identifier %in% names(vcf.details)) ) { 
            stop(error.message);
        }
    }
    
    optional.columns <- c(
        'job.dependency.column' = job.dependency.column, 
        'caller.column' = caller.column
    );
    
    for( column in names(optional.columns) ) { 
        
        column.identifier <- optional.columns[column];
        
        # if NA, skip ahead 
        if( is.na(column.identifier) ) next;
        
        # if not NA, subject to same quality input tests as required columns
        error.message <- paste(column, 'not found in vcf.details data frame');
        
        if( is.numeric(column.identifier) && column.identifier > ncol(vcf.details) ) { 
            stop(error.message);
        }	
        
        if( is.character(column.identifier) && !(column.identifier %in% names(vcf.details)) ) { 
            stop(error.message);
        }
        
    }
    
    ### MAIN ##################################################################
    
    vcf.specification <- data.frame(
        sample.id = vcf.details[, sample.id.column], 
        vcf = vcf.details[, vcf.column], 
        stringsAsFactors = FALSE
    );
    
    if( !is.na(job.dependency.column) ) { 
        vcf.specification$job.dependency <- vcf.details[, job.dependency.column];
    }
    
    if( !is.na(caller.column) ) { 
        vcf.specification$caller <- vcf.details[, caller.column];
    }
    
    
    verify.vcf.specification(vcf.specification);
    
    return(vcf.specification);
}
