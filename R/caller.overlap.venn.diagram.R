#'. Make Venn diagram of variant caller overlap
#' 
#' @param variants 
#' 	Data frame containing variants, typically from merge.variants function
#' @param file.name
#'  Name of output file
#'
#'
caller.overlap.venn.diagram <- function(variants, file.name) {
    
    
    ### INPUT TESTS ###########################################################
    
    if( !('caller' %in% names(variants) ) ) {
        stop('variant data frame must have a field callers');
    }
    
    ### MAIN ##################################################################
    
    all.callers <- stringr::str_split(variants$caller, pattern = ':');
    unique.callers <- unique( unlist(all.callers) );
    
    if( length(unique.callers) > 5) {
        stop('Cannot make a Venn diagram for more than 5 unique callers');
    }
    
    # create ID field to uniquely identify variants
    variants$id <- paste0(
        variants$sample.id, '-', 
        variants$CHROM, ':', 
        variants$POS, '-',
        variants$REF, '>', 
        variants$ALT
    );
    
    caller.results <- lapply(
        unique.callers, 
        function(caller, variants) {
            return(variants$id[ grepl(caller, variants$caller) ])
        }, 
        variants = variants
    );
    
    names(caller.results) <- capitalize.caller(unique.callers);

    # turn off log files
    if( requireNamespace('futile.logger', quietly = TRUE) ) {
        futile.logger::flog.threshold(futile.logger::ERROR, name = 'VennDiagramLogger');
    }
    
    colour.scheme <- c(
        '#0039A6', '#FF6319', '#6CBE45', '#996633', '#A7A9AC', 
        '#FCCC0A', '#B933AD', '#EE352E',  '#808183', '#00933C'
        );

    VennDiagram::venn.diagram(
        caller.results, 
        filename = file.name,
        fill = colour.scheme[ 1:length(unique.callers) ], 
        ext.text = FALSE,
        ext.percent = rep(0.01, 3),
        cat.pos = 0
        );
    
}