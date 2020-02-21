#' Prepare BAM specification data frame to standardized format for downstream analyses.
#'
#' @description
#' This function prepares a data frame that can be used to run variant callers. 
#' For matched normal variant calling, this data frame will contain three columns with names: sample.id, tumour.bam, normal.bam
#' For unpaired variant calling, the data frame will contain two columns with names: sample.id, tumour.bam
#'
#' @param sample.details Data frame where each row represents a sample to be run. 
#' 	Must contain sample ID, path to tumour BAM, and path to normal BAM.
#' @param paired Logical indicating whether the sample specification is for a paired analysis.
#' @param sample.id.column Index or string giving column of sample.details that contains the sample ID
#' @param tumour.bam.column Index or string giving column of sample.details that contains the path to the tumour BAM
#' @param normal.bam.column Index or string giving column of sample.details that contains the path to the normal BAM
#'
#' @return bam.specification Data frame with one row per sample to be run
#'
#'
prepare.bam.specification <- function(
    sample.details, 
    paired = TRUE,
    sample.id.column = 1,
    tumour.bam.column = 2,
    normal.bam.column = 3	
) { 
    
    ### INPUT TESTS ###########################################################	
    
    # Verify that tumour and normal BAM columns exist
    # They can be passed in as either an index or a string
    #	=> treat two cases separately
    
    tumour.bam.error.message <- 'tumour.bam.column not found in input data.';
    if( is.character(tumour.bam.column) && !(tumour.bam.column %in% names(sample.details)) ) { 
        stop(tumour.bam.error.message);		
    }
    
    if( is.numeric(tumour.bam.column) && tumour.bam.column > ncol(sample.details) ) { 
        stop(tumour.bam.error.message);	
    }
    
    normal.bam.error.message <- 'normal.bam.column not found in input data.';
    if( is.character(normal.bam.column) && !(normal.bam.column %in% names(sample.details)) ) {
        stop(normal.bam.error.message);
    }	
    
    if( is.numeric(normal.bam.column) && normal.bam.column > ncol(sample.details) ) {
        stop(normal.bam.error.message);
    }
    
    ### MAIN ##################################################################
    
    # make sure BAM specification is in a format we like
    
    bam.specification <- data.frame(
        sample.id = sample.details[, sample.id.column],
        tumour.bam = sample.details[, tumour.bam.column], 
        stringsAsFactors = FALSE
    );
    
    if(paired) { 
        bam.specification$normal.bam <- sample.details[, normal.bam.column];
    } 
    
    # verify sample sheet has correct format
    verify.bam.specification(bam.specification);
    
    return(bam.specification);
}

