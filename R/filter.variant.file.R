#' Filter variants in file.
#' 
#' @description
#' Filter variants from file, and save to output. Wrapper function that opens the 
#' variant file, calls filter.variants, and saves the result to file
#'
#' @param variant.file Path to variant file
#' @param output.file Path to output file
#' @inheritParams filter.variants
#'
#' @return None
#'
#'
filter.variant.file <- function(
    variant.file, 
    output.file,
    config.file = NULL,
    caller = c('vardict', 'ides', 'mutect', 'pgm', 'consensus')
    ) {
    
    caller <- match.arg(caller);
    
    ### INPUT TESTS ###########################################################	
    
    if( !file.exists(variant.file) ) { 
        error.message <- paste('File', variant.file, 'does not exist');
        stop(error.message);
    }
    
    ### MAIN ##################################################################
    
    
    variants <- utils::read.table(
        variant.file, 
        sep = '\t', 
        header = TRUE
        );
    
    filtered.variants <- filter.variants(
        variants, 
        caller = caller, 
        config.file = config.file
        );
    
    utils::write.table(
        filtered.variants, 
        output.file, 
        sep = '\t',
        row.names = FALSE
        );
    
}
