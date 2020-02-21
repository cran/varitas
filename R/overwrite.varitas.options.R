#' overwrite.varitas.options
#'
#' @description
#'	Overwrite VariTAS options with options provided in config file.
#'
#' @param config.file Path to config file that should be used to overwrite options
#'
#' @return None
#' @examples 
#' \dontrun{
#' config <- file.path(path.package('varitas'), 'config.yaml')
#' overwrite.varitas.options(config)
#' }
#'
#' @export
overwrite.varitas.options <- function(config.file) { 
    
    ### INPUT TESTS ###########################################################
    
    if( !file.exists(config.file) ) { 
        error.message <- paste('File', config.file, 'not found');
        stop(error.message);
    }
    
    ### MAIN ##################################################################
    
    config <- yaml::yaml.load_file(config.file);
    config$pkgname <- get.varitas.options('pkgname');

    # if a mode has been set, start by setting those settings
    # they will later be overwritten by 
    if( 'mode' %in% names(config) ) {

        # convert to lower case to allow users to specify ctDNA rather than ctdna
        config$mode <- tolower(config$mode);

        # only ctDNA and tumour supported
        if( !( config$mode %in% c('ctdna', 'tumour') ) ) {
            stop('mode must be either ctDNA or tumour');
        }

        # read mode defaults from file
        mode.default.file <- system.file(
            paste0(config$mode, '_defaults.yaml'),
            package = get.varitas.options('pkgname')
            );

        mode.defaults <- yaml::yaml.load_file( mode.default.file );

        # update mode defaults with specified values
        config <- utils::modifyList(mode.defaults, config);
    }

    # update filters so defaults are considered as baseline for each caller
    if( 'filters' %in% names(config) ) {
        config$filters <- get.filters(config$filters);
    }
       
    options(varitas = config);   
}
