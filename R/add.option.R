#' add.option 
#'
#' @description
#'  Add option to nested list of options. Applied recursively
#' 
#' @param name 
#'  Option name. Nesting is indicated by character specified in nesting.character.
#' @param value 
#'  New value of option
#' @param old.options 
#'  Nested list the option should be added to
#' @param nesting.character
#'  String giving Regex pattern of nesting indication string. Defaults to '\\.'
#'
#' @return Nested list with updated options
#'
add.option <- function(name, value, old.options, nesting.character = '\\.') { 
    
    # split name of option into the immediate level and subsequent levels
    name.components <- stringr::str_split(
        name, 
        pattern = nesting.character,
        n = 2
    )[[ 1 ]];

    option.level.name <- name.components[1];
    new.options <- old.options;

    # if trying to add filters, need to use get.filter to expand default filters
    # treat this case separately and return
    if( 'filters' == name ) {
        new.options$filters <- get.filters(value);
        return(new.options); 
    }
    
    if( 1 == length(name.components) ) { 
        # if we're at a leaf - add to list 
        new.options[[ option.level.name ]] <- value;
    } else { 
        # not at a leaf... recursion!
        if( !(option.level.name %in% names(old.options) ) ) { 
            error.message <- paste('Key', option.level.name, 'not found in list'); 
            stop(error.message);
        }
        
        new.options[[ option.level.name ]] <- add.option(
            name.components[2], 
            value,
            new.options[[ option.level.name ]]
        );        
    }
    
    return(new.options);
}