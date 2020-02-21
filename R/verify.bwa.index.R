#' verify.bwa.index
#'
#' @inheritParams verify.fasta.index
#'
#' @description
#' 	 Verify that bwa index files exist for a fasta file
#'
#' @return index.files.exist Logical indicating if bwa index files were found (only returned if error set to FALSE)
#'
#'
verify.bwa.index <- function(fasta.file, error = FALSE) {
    
    bwa.index.files <- paste0(fasta.file, '.', c('amb', 'ann', 'bwt', 'pac', 'sa'));
    index.files.exist <- all( file.exists(bwa.index.files) );
    
    if( error ) {
        
        if(!index.files.exist) {
            error.message <- paste('Index files not found for reference genome file', fasta.file, '- try running bwa index.');
            stop(error.message);
        }
        
    } else {
        
        return(index.files.exist);
    }
    
    
}