#' prepare.fastq.specification 
#'
#' @description
#'  Prepare FASTQ specification data frame to standardized format for downstream analyses.
#'
#' @details
#'  This function prepares a data frame that can be used to run alignment.
#'  For paired-end reads, this data frame will contain three columns with names: sample.id, reads, mates
#'  For single-end reads, the data frame will contain two columns with names: sample.id, reads
#'
#' @param sample.details 
#'	Data frame where each row represents a sample to be run. Must contain sample ID, path to tumour BAM, and path to normal BAM.
#' @param sample.id.column 
#'  Index or string giving column of \code{sample.details} that contains the sample ID
#' @param fastq.columns 
#'  Index or string giving column(s) of \code{sample.details} that contain path to FASTQ files 
#' @param patient.id.column
#'  Index or string giving column of \code{sample.details} that contains the patient ID
#' @param tissue.column
#'  Index or string giving column of \code{sample.details} that contains information on tissue (tumour/ normal)
#'
#' @return  Data frame with one row per sample to be run
#'
#'
prepare.fastq.specification <- function(
    sample.details,
    sample.id.column = 1,
    fastq.columns = c(2, 3), 
    patient.id.column = NA,
    tissue.column = NA
    ) { 
    
    ### INPUT TESTS ###########################################################
    
    # Note: Need to do some input tests here to avoid errors in assembling 
    # the data frame to start with.
    # More tests are run in verification step 
    
    
    # Verify that FASTQ columns exists
    # They can be passed in as either an index or a string
    #	=> treat two cases separately
    fastq.error.message <- 'FASTQ columns not found in input data.';
    
    if( is.character(fastq.columns) && !all(fastq.columns %in% names(sample.details)) ) { 
        stop(fastq.error.message);		
    }
    
    if( is.numeric(fastq.columns) && max(fastq.columns) > ncol(sample.details)) { 
        stop(fastq.error.message);
    }
    
    if( length(fastq.columns) > 2 ) { 
        stop('Can accept at most two FASTQ columns.');
    }
    
    
    ### MAIN ##################################################################
    
    fastq.specification <- data.frame(
        'sample.id' = sample.details[, sample.id.column], 
        'reads' = sample.details[, fastq.columns[1]], 
        stringsAsFactors = FALSE
    );
    
    if( 2 == length(fastq.columns) ) {
        fastq.specification[, 'mates'] <- sample.details[, fastq.columns[2]];
    }  
    
    if( !is.na(patient.id.column) ) {
        fastq.specification[, 'patient.id'] <- sample.details[, patient.id.column];
    }
    
    if( !is.na(tissue.column) ) {
        fastq.specification[, 'tissue'] <- sample.details[, tissue.column];
    }
    
    
    verify.fastq.specification(fastq.specification);
    
    return(fastq.specification);
}
