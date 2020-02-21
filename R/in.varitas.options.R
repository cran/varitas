#' Check if a key is in VariTAS options
#'
#' @param option.name 
#'  String giving name of option (with different levels joined by \code{nesting.character})
#' @param varitas.options
#'  Ampliseq options as a list. If missing, they will be obtained from \code{get.varitas.options()} 
#' @inheritParams get.varitas.options
#'
#' @return in.options Boolean indicating if the option name exists in the current varitas options
in.varitas.options <- function(
    option.name = NULL, 
    varitas.options = NULL, 
    nesting.character = '\\.'
    ) { 
    
    option.name <- gsub(nesting.character, '\\.', option.name);
    
    if( is.null(varitas.options) ) { 
        varitas.options <- get.varitas.options();
    } 
    
    # split name of option into the immediate level and subsequent levels
    name.components <- stringr::str_split(
        option.name, 
        pattern = nesting.character,
        n = 2
    )[[ 1 ]];
    
    option.level.name <- name.components[1];
    
    if(	option.level.name %in% names(varitas.options) ) { 
        
        if( 1 == length(name.components) ) { 
            # we have reached a leaf, and all components were in the options
            #	=> The option exists!
            in.options <- TRUE;
        } else { 
            # we are not at a leaf yet, try next level
            in.options <- in.varitas.options(
                option.name = name.components[2], 
                varitas.options = varitas.options[[ option.level.name ]],
                nesting.character = nesting.character
            );			
        }
        
    } else { 
        # key not found -> return FALSE
        in.options <- FALSE;
    }
    
    
    return(in.options);
    
}