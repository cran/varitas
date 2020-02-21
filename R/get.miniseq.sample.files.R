#' get.miniseq.sample.files
#'
#' @description
#' 	Get files for a sample in a directory, ensuring there's only a single match per sample ID.
#'
#' @param sample.ids Vector of sample ids. Should form first part of file name
#' @param directory Directory where files can be found
#' @param file.suffix Regex expression for end of file name. For example, `file.suffix = '_S\\d{1,2}_.*_R1_.*'` will match R1 files.1 files.
#'
#' @return Character vector of file paths
#'
get.miniseq.sample.files <- function(sample.ids, directory, file.suffix = '_S\\d{1,2}_.*') {


	files <- sapply(
		sample.ids,
		function(sample.id, directory, file.suffix) {

			# MiniSeq sometimes converts underscores to dashes
			# check for either
			sample.id.regex <- paste0('(', sample.id, '|', gsub('_', '-', sample.id), ')');

			sample.files <- list.files(
				pattern = paste0('^', sample.id.regex, file.suffix),
				path = directory, 
				full.names = TRUE
				);

			# be strict – 
			if( length(sample.files) > 1) {
				error.message <- paste(
					'More than one file found in', directory, 'for sample', sample.id, '\n',
					paste(sample.files, collapse = ' ')
					);
				stop(error.message);
			}

			if( 0 == length(sample.files) ) {
				error.message <- paste('No file found in', directory, 'for sample', sample.id);
				stop(error.message);
			}

			return(sample.files);
			},
			directory = directory,
			file.suffix = file.suffix
		);

	return(files);
}
