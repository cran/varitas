#' Plot amplicon coverage by genome order
#' 
#' @description
#' Use values obtained by bedtools coverage to make a plot of coverage by genome order 
#'
#' @param coverage.data data frame with results from bedtools coverage command
plot.coverage.by.genome.order <- function(coverage.data) {
		
	# TO DO:
	#	- debug pool colour part of this function

	# make sure it is ordered by genome location
	coverage.data <- coverage.data[order(coverage.data[, 1], coverage.data[, 2]), ];

	# not clear how many columns exist in this, as that will depend on the panel BED file
	reads.mapped.to.amplicon <- coverage.data[, ncol(coverage.data) - 3];

	# if pool information is avaiable, colour by it
	pools <- get.pool.from.panel.data(coverage.data);

	if(is.null(pools)) {
		point.colours <- "black";

	} else {

		# create a colour scheme to recode pool vector by
		unique.pools <- unique(pools);

		colour.scheme <- get.colours(length(unique.pools));
		names(colour.scheme) <- unique.pools;

		point.colours <- colour.scheme[pools];
	}

	graphics::plot(
		reads.mapped.to.amplicon,
		xlab = "Genome order",
		ylab = "Coverage",
		pch = 16,
		col = point.colours
		);
}




