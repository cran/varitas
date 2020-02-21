#' Return VariTAS settings
#'
#' @param option.name Optional name of option. If no name is supplied, the full list of VariTAS options will be provided. 
#' @param nesting.character String giving Regex pattern of nesting indication string. Defaults to '\\.'
#' 
#' @return varitas.options list specifying VariTAS options
#' 
#' @examples 
#'	reference.build <- get.varitas.options('reference_build');
#'	mutect.filters <- get.varitas.options('filters.mutect');
#'
#' @export
get.varitas.options <- function(option.name = NULL, nesting.character = '\\.') { 
    
    if( is.null(option.name) ) { 
        varitas.options <- getOption('varitas');
        
    } else { 
        
        varitas.options <- get.option(
            option.name,
            nesting.character = nesting.character
        );
        
    }
    
    return( varitas.options );
}