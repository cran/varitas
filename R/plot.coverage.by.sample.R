#' plot.coverage.by.sample
#'
#' @description
#'	Make a barplot of coverage per sample
#'
#' @param coverage.sample Data frame of coverage data, typically from \code{get.coverage.by.sample.statistics}
#' @param file.name Name of output file
#' @param statistic Statistic to be plotted (mean or median)
#' 
#' @return None
plot.coverage.by.sample <- function(
	coverage.sample, 
	file.name,
	statistic = c('mean', 'median')
	) {

	statistic <- match.arg(statistic);

	### INPUT TESTS ###########################################################

	statistic.column <- paste0(statistic, '.coverage');

	if( !( statistic.column %in% names(coverage.sample) ) ) {
		error.message <- paste(
			'coverage.sample data frame does not contain a column',
			statistic.column
			);

		stop(error.message);
	}

	### MAIN ##################################################################

	barplot.data <- coverage.sample[, statistic.column];
	names(barplot.data) <- coverage.sample$sample.id;

	cex.axis <- 0.7;
  if( nrow(coverage.sample) > 30) cex.axis <- 0.4;

    yaxis.label <- c(
    	'mean' = 'Mean coverage',
    	'median' = 'Median coverage'
    	);

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
		mar = c(5.3, 4, 0.5, 0.5),
		cex.axis = cex.axis,
		font.axis = 2,
		oma = c(0, 0, 0, 0),
		tcl = -0.2,
		las = 2,
		mgp = c(3, 0.25, 0)
		);

	graphics::barplot(
		barplot.data,
		ylab = yaxis.label[ statistic ]
		);

	if( !is.null(file.name) ) grDevices::dev.off();
}