#' Make barplot of variants per caller
#' 
#' @param variants 
#'  Data frame with variants
#' @param file.name 
#'  Name of output file
#' @param group.by
#'  Optional grouping variable for barplot
#'  
#' @return None
#'
#'
#'
#'
#' @importFrom magrittr "%>%"
#'
#'
variants.caller.barplot <- function(
    variants, 
    file.name, 
    group.by = NULL
    ) {
    
    ### INPUT TESTS ###########################################################
    
    if( !all(c('sample.id', 'REF', 'ALT') %in% names(variants)) ) { 
        stop('Variant data frame missing required columns');
    }
    
    ### MAIN ##################################################################
    
    # colour scheme from COSMIC
    trinucleotide.colours <- c(
        'C>A' = '#5CBDEB',
        'C>G' = '#050708',
        'C>T' = '#D23C32',
        'T>A' = '#CBCACB',
        'T>C' = '#ABCD72',
        'T>G' = '#E7C9C6'
        );

    if( !is.null(group.by) && 'substitution' == group.by && !('substitution' %in% names(variants)) ) {
        variants$substitution <- get.base.substitution(ref = variants$REF, alt = variants$ALT);
    } else if( !is.null(group.by) && 'type' == group.by && !('type' %in% names(variants))) {
        variants$type <- classify.variant(ref = variants$REF, alt = variants$ALT);
    }
    
    split.variants <- split.on.column(
        variants, 
        column = 'caller', 
        split.character = ':'
    );
    
    # capitalized properly
    split.variants$caller <- capitalize.caller(split.variants$caller);
    
    grouping.variables <- 'caller';
    if( !is.null(group.by) ) {
        grouping.variables <- c(grouping.variables, group.by);
    }
    
    counts <- split.variants %>% 
        dplyr::group_by_at(grouping.variables) %>%
        dplyr::summarise(n = n());
    
    barplot.data <- as.matrix( t(stats::xtabs(n ~ ., counts)) );
    

    colour.scheme <- c("#0039A6", "#FF6319", "#6CBE45", "#996633", "#A7A9AC", "#FCCC0A");
    if( !is.null(group.by) && 'substitution' == group.by ) {
        colour.scheme <- trinucleotide.colours[ rownames(barplot.data) ];
    }

	if( !is.null(file.name) ) {
		grDevices::png(
			file.name,
			width = 6,
			height = 5,
			units = 'in',
			res = 300
			);
	}
	
	graphics::par(
		mar = c(4.5, 4, 0.5, 5),
		cex.axis = 0.9,
		font.axis = 2,
		oma = c(0, 0, 0, 0),
		tcl = -0.2,
		las = 2,
		mgp = c(3, 0.25, 0)
		);

	graphics::barplot(
		barplot.data,
		ylab = 'Variants detected',
		col = colour.scheme,
		legend.text = rownames(barplot.data),
		args.legend = list(
            x = graphics::par('usr')[2],
            y = graphics::par('usr')[4]/2,
            xjust = 0,
            yjust = 0.5,
            xpd = TRUE,
            bty = 'n'
            )
		);

	if( !is.null(file.name) ) grDevices::dev.off();
    
}
