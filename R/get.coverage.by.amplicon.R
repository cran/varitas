#' Process sample coverage per amplicon data
#'
#' @description
#' Parse coverageBed output to get coverage by amplicon
#'
#' @references \url{http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html}
#'
#' @inheritParams get.coverage.by.sample.statistics
#'
#' @return combined.data Data frame giving coverage per amplicon per sample.
#'
#'
get.coverage.by.amplicon <- function(project.directory) {
    
    # This parses the bedtools coverage step (run without -hist option) 
    # NOTE: might be able to use bedr for actual bedtools commands, but it's currently not working
    #	
    # After each interval in A, bedtools coverage will report:
    #	- The number of features in B that overlapped (by at least one base pair) the A interval.
    #	- The number of bases in A that had non-zero coverage from features in B.
    #	- The length of the entry in A.
    #	- The fraction of bases in A that had non-zero coverage from features in B.
    #
    # These are the last four columns. Earlier columns will depend on the format of the panel BED file.
    
    coverage.paths <- system.ls(pattern = "*/*.sort.txt", directory = project.directory, error = TRUE);
    sample.ids <- extract.sample.ids(coverage.paths, from.filename = TRUE);
    
    
    combined.data <- data.frame();
    
    for( i in 1:length(coverage.paths) ) { 
        
        path <- coverage.paths[i];
        sample.id <- sample.ids[i];
        
        coverage.data <- utils::read.delim(
            path,
            sep = "\t",
            as.is = TRUE,
            header = FALSE
            );
        
        reads.mapped.to.amplicon <- coverage.data[, ncol(coverage.data) - 3];
        
        if(1 == i) {
            combined.data <- coverage.data[, 1:(ncol(coverage.data) - 4)];
            names(combined.data)[1:3] <- c('chr', 'start', 'end');
        }
        
        #? do we need input checks here?
        #	checking that order matches will slow things down, and mismatch should not happen as long as coverage 
        combined.data[, sample.id] <- reads.mapped.to.amplicon;
        
    }
    
    return(combined.data);
    
}