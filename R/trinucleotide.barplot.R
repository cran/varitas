#' Make barplot of trinucleotide substitutions
#' 
#' @param variants Data frame with variants
#' @param file.name Name of output file
#'
#' @return None
#'
#'
#'
#' @importFrom magrittr "%>%"
#'
#'
trinucleotide.barplot <- function(variants, file.name) {
    
    ### INPUT TESTS ###########################################################
    
    if( !all(c('sample.id', 'REF', 'ALT') %in% names(variants)) ) { 
        stop('Variant data frame missing required columns');
    }
    
    ### MAIN ##################################################################
    
    variants$substitution <- get.base.substitution(ref = variants$REF, alt = variants$ALT);
    
    counts <- variants %>% 
        dplyr::group_by(substitution) %>%
        dplyr::summarise(n = n());
    
    barplot.data <- as.matrix( t( stats::xtabs(n ~ ., counts )) );
    
    
	if( !is.null(file.name) ) {
		grDevices::png(
			file.name,
			width = 7,
			height = 5,
			units = 'in',
			res = 400
			);
	}

	graphics::par(
		mar = c(4, 4, 0.5, 0.5),
		cex.axis = 0.8,
		font.axis = 2,
		oma = c(0, 0, 0, 0),
		tcl = -0.2,
		las = 2,
		mgp = c(3, 0.25, 0)
		);

	graphics::barplot(
		barplot.data,
		ylab = 'Variants detected'
		);

	if( !is.null(file.name) ) grDevices::dev.off();
    
}
