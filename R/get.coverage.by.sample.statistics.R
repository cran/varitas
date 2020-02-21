#' Get statistics about coverage per sample
#'
#' @param project.directory Path to project directory. Each sample should have its own subdirectory
#' 
#' @return coverage.by.sample.statistics Data frame with coverage statistics per sample
#' 
#'
get.coverage.by.sample.statistics <- function(project.directory) {
    
    ### INPUT TESTS ###########################################################
    
    if( !is.character(project.directory) || length(project.directory) > 1 ) {
        stop('project.directory must be a string.');
    }
    
    if( !file.exists(project.directory) ) {
        error.message <- paste('Directory', project.directory, 'does not exist');
        stop(error.message);
    }
    
    
    ### MAIN ##################################################################
    
    # gather data from individual files
    total.coverage.statistics <- process.total.coverage.statistics(project.directory);
    coverage.report.data <- process.coverage.reports(project.directory);
    
    # merge – this has to be done pairwise
    # does the sample order have to match Ros' ? 
    coverage.by.sample.statistics <- merge(
        total.coverage.statistics,
        coverage.report.data,
        by = "sample.id",
        all = TRUE
    ); 
    
    return(coverage.by.sample.statistics);
    
}
