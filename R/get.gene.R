#' get.gene
#' 
#' @description
#' 	Use guesswork to extract gene from data frame of targeted panel data. The panel
#'  designer output can change, so try to guess what the format is.
#'
#' @param bed.data Data frame containing data from bed file
#' 
#' @return vector of gene names, one entry for each row of \code{bed.data}
#'
#'
get.gene <- function(bed.data) {

	### INPUT TESTS ###########################################################

	if( !is.data.frame(bed.data) ) {
		stop('bed.data must be a data frame');
	}

	### MAIN ##################################################################

	# figure out which column contains gene
	gene.column <- 5;
	if( any( grepl('GENE_ID', bed.data[, 6] )) ) {
		gene.column <- 6;
	}

	# this needs to be more sophisticated
	if( all( grepl('GENE_ID', bed.data[, gene.column]) ) ) {
		genes <- gsub('(.*)GENE_ID=(.+?)(_|;|$)(.*)', '\\2', bed.data[, gene.column]);
	} else {
		genes <- sapply(
			strsplit(bed.data[, gene.column], split = ' '),
			function(x) x[1]
			);
	}

	return(genes);
}