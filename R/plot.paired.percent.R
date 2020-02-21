#' plot.paired.percent
#'
#' @description
#'	Make a barplot of percent paired reads per sample
#'
#' @param coverage.sample Data frame of coverage data, typically from \code{get.coverage.by.sample.statistics}
#' @param file.name Name of output file
#'
#' @return None
plot.paired.percent <- function(coverage.sample, file.name) {

	### INPUT TESTS ###########################################################

	if( !('paired.percent' %in% names(coverage.sample) ) ) {
		stop('coverage.sample data frame does not contain a column paired.percent');
	}

	### MAIN ##################################################################

	barplot.data <- 100*coverage.sample$paired.percent;
	names(barplot.data) <- coverage.sample$sample.id;

	cex.axis <- 0.7;
    if( nrow(coverage.sample) > 30) cex.axis <- 0.4;

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
		mar = c(5.1, 4, 0.5, 0.5),
		cex.axis = cex.axis,
		font.axis = 2,
		oma = c(0, 0, 0, 0),
		tcl = -0.2,
		las = 2,
		mgp = c(3, 0.25, 0)
		);

	graphics::barplot(
		barplot.data,
		ylab = 'Paired reads (%)'
		);

	if( !is.null(file.name) ) grDevices::dev.off();

}