#' Read variant calls from file and format for ease of downstream analyses.
#'
#' @param variant.file
#' 	Path to variant file.
#' @param variant.caller
#'	String indicating which variant caller was used. Needed to format the headers. 
#'
#' @return variant.calls Data frame of variant calls
#'
#'
read.variant.calls <- function(variant.file, variant.caller) { 
    
    ### INPUT TESTS ###########################################################
    
    if( is.factor(variant.file) ) variant.file <- as.character(variant.file);
    
    if( !is.character(variant.file) ) {
        stop('variant.file must be a string.');
    }
    
    if( !file.exists(variant.file) ) { 
        error.message <- paste('File', variant.file, 'does not exist');
        stop(error.message);
    }
    
    
    ### MAIN ##################################################################
    
    # get sample ID
    sample.id <- gsub('(.*)\\.passed\\.ontarget\\.vcf\\.annovar.*', '\\1', basename(variant.file) );

    uppercase.caller <- toupper(variant.caller);
    
    variant.calls <- utils::read.table(
        variant.file,
        header = TRUE,
        sep = '\t', 
        stringsAsFactors = FALSE,
        quote = ""
        );

    # if no variants in file, skip ahead			
    if( 0 == nrow(variant.calls) ) {
        warning.message <- paste('File', variant.file, 'exists, but does not contain any variants');
        warning(warning.message);
        
        return(NULL);
    }	

    # make sure REF and ALT are character vectors, don't want any accidental T -> TRUE conversion e
    variant.calls$REF <- logical.to.character(variant.calls$REF);
    variant.calls$ALT <- logical.to.character(variant.calls$ALT);

    
    names(variant.calls) <- fix.names(names(variant.calls), variant.caller, sample.id = sample.id);
    
    if(nrow(variant.calls) > 0) {
        variant.calls[, paste0(uppercase.caller, '.VARIANT.READS')] <- variant.calls[, paste0(uppercase.caller, '.TUMOUR.AD')];
    }
    
    return(variant.calls);	
    
}
