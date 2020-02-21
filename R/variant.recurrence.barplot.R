#' Make barplot of variants per caller
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
variant.recurrence.barplot <- function(variants, file.name) {

	### INPUT TESTS ###########################################################

	if( !is.data.frame(variants) ) {
		stop('variants must be a data frame');
	}

	if( !all(c('CHROM', 'REF', 'ALT', 'caller') %in% names(variants))) {
		stop('variants must contain the columns CHROM, REF, ALT, and caller');
	}

	### MAIN ##################################################################

	# add a label for each variant
	variants$position.id <- paste0(variants$CHROM, ':', variants$POS, ' | ', variants$REF, '->', variants$ALT);

	# prepare data for barplot
	position.counts <- variants %>% 
		dplyr::group_by(position.id, caller) %>% 
		dplyr::summarise(n = n()) %>%
		dplyr::ungroup() %>% 
		tidyr::complete(position.id, caller, fill = list(n = 0));

	counts.wide <- position.counts %>% tidyr::spread(position.id, n);

	counts.matrix <- as.matrix(counts.wide[, 'caller' != names(counts.wide) ]);
	rownames(counts.matrix) <- counts.wide$caller;
  
	if (nrow(counts.matrix) > 1) { # Should only be ordered if there is more than one row
	  counts.matrix <- counts.matrix[, order(colSums(counts.matrix), decreasing = TRUE) ];
	}

	# colour scheme
	# use London underground colours
	colour.scheme <- c(
		'Central' = '#E32017', 'Circle' = '#FFD300', 'Hammersmith and City' = '#F3A9BB', 
		'Jubilee' = '#A0A5A9', 'Waterloo and City' = '#95CDBA', 'Metropolitan' = '#9B0056', 
		'Northern' = '#000000', 'Piccadilly' = '#003688', 
		'DLR' = '#00A4A7', 'Overground' = '#EE7C0E', 'Victoria' = '#0098D4', 
		'Tramlink' = '#84B817', 'Cable Car' = '#E21836', 'Crossrail' = '#7156A5', 
		'District' = '#00782A', 'Bakerloo' = '#B36305'
		);

	if( !is.null(file.name) ) {
		grDevices::png(
			file.name,
			width = 8,
			height = 5,
			units = 'in',
			res = 300
			);
	}
	
	graphics::par(
		mar = c(8, 4, 0.5, 8),
		cex.axis = 0.7,
		font.axis = 2,
		oma = c(0, 0, 0, 0),
		tcl = -0.2,
		las = 2,
		mgp = c(3, 0.25, 0)
		);

	graphics::barplot(
		counts.matrix[, 1:min(25, ncol(counts.matrix)) ],
		ylab = 'Samples with variant',
		col = colour.scheme[ 1:nrow(counts.matrix) ],
		legend.text = rownames(counts.matrix),
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
