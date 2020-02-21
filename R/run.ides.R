#' Run iDES
#'
#' @details 
#' Run iDES step 1on each sample, to tally up calls by strand.
#' Files are output to a the sample subdirectory
#'
#' @param project.directory
#'  Directory containing files
#' @param sample.id.pattern
#'  Regex pattern to match sample IDs
#' @param sample.ids
#'  Vector of sample IDs
#' @param job.dependencies
#'  Vector of job dependencies
#'
#' @return None
#'
#' @note
#'  Deprecated function for running iDES. 
#'  Follows previous development package without specification data frames
#'
#' @references \url{https://cappseq.stanford.edu/ides/}
#'
run.ides <- function(
    project.directory, 
    sample.id.pattern = "._S\\d+$", 
    sample.ids = NULL, 
    job.dependencies = NULL 
    ) { 
    
    # TO DO:
    #   - should wrapper functions that submit jobs be named in a way that indicates this? 
    
    ### INPUT TESTS ###########################################################
    
    if( !is.character(project.directory) || length(project.directory) > 1 ) {
        stop('project.directory should be a single string giving path to project directory');
    }
    
    if( !dir.exists(project.directory) ) {
        stop( paste(project.directory, 'is not a directory.') );
    }
    
    for(sample.id in sample.ids) {
        
        if( !dir.exists(file.path(project.directory)) ) {
            stop.message <- paste('Subdirectory', sample.id, 'does not exist in', project.directory);
            stop(stop.message);
        }
        
    }
    
    ### MAIN ##################################################################
    
    # if sample IDs have not been explicitly named, find all the directories
    # matching pattern in project directory
    if( is.null(sample.ids) ) {
        sample.ids <- dir(
            path = project.directory, 
            pattern = sample.id.pattern
        );
    }
    
    # TO DO: change this if we change the name!
    iDES.script <- system.file("perl", "run_iDES_sample.pl", package = getOption('')$pkgname );
    
    config.file <- save.config();
    
    # loop over samples to submit alignment job for each sample
    for(sample.id in sample.ids) {
        
        command <- make.command.line.call(
            main.command = c("perl", iDES.script), 
            options = list(
                "project_directory" = project.directory,
                "sample_id" = sample.id, 
                "config_file" = config.file,
                "job_dependencies" = job.dependencies
            )
        );
        
        # submit to system
        system(command);
    }
    
}