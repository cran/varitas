#' Set options for varitas pipeline.
#'
#' @description
#' 	Set or overwrite options for the VariTAS pipeline. Nested options should be separated by a dot. 
#' 	For example, to update the reference genome for grch38, use reference_genome.grch38
#'
#' @param \dots options to set
#' @return None
#'
#' @examples
#' \dontrun{
#'	set.varitas.options(reference_build = 'grch38'); 	
#' 	set.varitas.options(
#'		filters.mutect.min_normal_depth = 10,
#'		filters.vardict.min_normal_depth = 10
#'		);
#' }
#' @export
set.varitas.options <- function(...) {
    # TO DO:
    # 	- implement this 'options'-style, i.e. allow for arguments passed directly
    
    updated.options <- get.varitas.options();
    
    options.to.add <- list(...);


    # currently don't support default filters - not sure how I could implement those
    if( any( grepl('^filters.default', names(options.to.add) )) ) {

        error.message <- paste(
            'Currently cannot set default filters with set.varitas.options.',
            'Please use overwrite.varitas.options instead'
            );

        stop( error.message );
    }

    # if one of the settings is mode, update that first so other settings will take precedence
    if( 'mode' %in% names(options.to.add) ) {

        option.value <- tolower( options.to.add[[ 'mode' ]] );

        # only ctDNA and tumour supported
        if( !( option.value %in% c('ctdna', 'tumour') ) ) {
            stop('mode must be either ctDNA or tumour');
        }

        # warn user that settings are being overwritten
        warning.message <- paste(
            'Setting mode to', 
            option.value, 
            '- overwriting any previous filters'
            );

        warning( warning.message );

        # read mode defaults from file
        mode.default.file <- system.file(
            paste0(option.value, '_defaults.yaml'),
            package = get.varitas.options('pkgname')
            );

        mode.defaults <- yaml::yaml.load_file( mode.default.file );

        # update settings for each 
        for(setting.name in names(mode.defaults) ) {

            updated.options <- add.option(
                name = setting.name,
                value = mode.defaults[[ setting.name ]],
                old.options = updated.options
                );
        }

        options.to.add <- options.to.add[ 'mode' != names(options.to.add) ];

    }
    
    for( i in seq_along(options.to.add) ) { 
        
        option.name <- names(options.to.add)[i];
        option.value <- options.to.add[[i]]; 

        updated.options <- add.option(
            name = option.name,
            value = option.value,
            old.options = updated.options
            );

    }
    
    options(varitas = updated.options);
}