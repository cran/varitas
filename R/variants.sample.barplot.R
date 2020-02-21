 #' Make barplot of variants per sample
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
variants.sample.barplot <- function(
    variants, 
    file.name
    ) { 
    
    ### INPUT TESTS ###########################################################
    
    if( !all(c('sample.id', 'Type') %in% names(variants)) ) { 
        stop('Variant data frame missing required columns');
    }
    
    ### MAIN ##################################################################
    
    variants$Type <- factor(
        variants$Type,
        levels = c('SNV', 'MNV', 'indel')
        );


    counts <- variants %>% 
        dplyr::group_by(sample.id, Type) %>%
        dplyr::summarise(n = n()) %>% 
        dplyr::ungroup() %>% 
        tidyr::complete(sample.id, Type, fill = list(n = 0));
    
    counts.wide <- counts %>% tidyr::spread(sample.id, n);
        
    barplot.data <- as.matrix( counts.wide[, 'Type' != names(counts.wide)] );
    rownames(barplot.data) <- counts.wide$Type;
    totals <- apply(barplot.data, 2, sum)
    max.value <- max(totals)

    # try to select a decent value for the axis font size
    cex.axis <- 1.2;
    if( length( unique(counts$sample.id) ) > 30) cex.axis <- 0.6;

    print(barplot.data);

    if( !is.null(file.name) ) {
		grDevices::png(
			file.name,
			width = 9,
			height = 5,
			units = 'in',
			res = 300
			);
	}
	
	graphics::par(
		mar = c(4.5, 2, 0.5, 0.5),
		cex.axis = 0.7,
		font.axis = 1,
		oma = c(0, 0, 0, 0),
		tcl = -0.2,
		las = 2,
		mgp = c(3, 0.25, 0)
		);

	graphics::barplot(
		barplot.data,
		ylim = c(0, max.value * 1.15),
		ylab = 'Variants detected',
		col = c('#E32017', '#FFD300', '#F3A9BB'),
		legend.text = rownames(barplot.data),
		args.legend = list(
		        x = 'top',
            #x = graphics::par('usr')[2],
            #y = graphics::par('usr')[4]/2,
            ncol = 3,
            xjust = 0,
            yjust = 0.5,
            xpd = TRUE,
            bty = 'n'
            )
		);

	if( !is.null(file.name) ) grDevices::dev.off();
    
}