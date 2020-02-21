#' verify.fasta.index
#'	
#' @description
#' 	 Verify that fasta index files exist for a given fasta file.
#'
#' @param fasta.file Fasta file to check
#' @param error Logical indicating whether to throw an (informative) error if verification fails
#'
#' @return faidx.exists Logical indicating if fasta index files were found (only returned if error set to FALSE)
#'
#'
verify.fasta.index <- function(fasta.file, error = FALSE) {
    
    faidx.file <- paste0(fasta.file, '.fai');
    
    faidx.exists <- file.exists(faidx.file);
    
    if( error ) {
        
        assertthat::assert_that(
            file.exists(faidx.file), 
            msg = paste('Fasta index file not found for file', fasta.file,  '\nTry running samtools faidx.')
        );
        
    } else {
        return(faidx.exists);
    }
    
}