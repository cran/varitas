#' verify.sequence.dictionary
#' 
#' @inheritParams verify.fasta.index
#'
#' @description
#' 	 Verify that sequence dictionary exists for a fasta file.
#'
#' @return dict.exists Logical indicating if sequence dictionary files were found (only returned if error set to FALSE)
#'
#'
verify.sequence.dictionary <- function(fasta.file, error = FALSE) {
    
    fasta.file.extension <- tools::file_ext(fasta.file);
    
    dict.file <- gsub(
        pattern = paste0('.', fasta.file.extension, '$'), 
        replacement = '.dict', 
        fasta.file
    );
    
    dict.exists <- file.exists(dict.file);
    
    if(error) {
        
        if( !dict.exists ) {
            error.message <- paste('Sequence dictionary not found for file', fasta.file,  '- try running GATK CreateSequenceDictionary.');
            stop(error.message);
        }
        
    } else {
        return(dict.exists);
    }
    
}