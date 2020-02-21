#' Helper function to recursively get an VariTAS option
#' 
#' @param name Option name
#' @param varitas.options Optional list of options to search in
#' @param nesting.character String giving Regex pattern of nesting indication string. Defaults to '\\.'
#'
#'
#' @return value Requested option
get.option <- function(name, varitas.options = NULL, nesting.character = '\\.') { 
    
    # if varitas.options not specified, get from settings
    if( is.null(varitas.options) ) { 
        varitas.options <- get.varitas.options();
    }
    
    # split name of option into the immediate level and subsequent levels
    name.components <- stringr::str_split(
        name, 
        pattern = nesting.character,
        n = 2
    )[[ 1 ]];
    
    option.level.name <- name.components[1];
    
    if( 1 == length(name.components) ) { 
        # we've reached our destination – get requested option
        value <- varitas.options[[ option.level.name ]];
    } else { 
        
        if( !(option.level.name %in% names(varitas.options)) ) { 
            stop('Key not found in list.');
        }
        
        value <- get.option(
            name.components[2], 
            varitas.options = varitas.options[[ option.level.name ]]
        );
    }
    
    return(value);
}