#' Fix variant call column names
#'
#' @description
#'	Fix headers of variant calls to prepare for merging. This mostly 
#'	consists in making sure the column headers will be unique by prefixing
#'	the variant caller in question.
#'
#' @param column.names Character vector of column names
#' @param variant.caller String giving name of variant caller
#' @param sample.id Optional sample ID. Used to fix headers.  
#'
#'
#' @return new.column.names Vector of column names after fixing]
#'
fix.names <- function(column.names, variant.caller, sample.id = NULL) {

	### INPUT TESTS ###########################################################
	
	if( !is.character(column.names) ) { 
		stop('column.names must be a character vector');
	}
	
	### MAIN ##################################################################

	variant.caller <- toupper(variant.caller);
	
	replacement.key <- c(
		'DP' = 'DEPTH',
		'TUMOUR\\.match' = 'NORMAL',
		# TPU sample IDs
		'^[A-Z]\\d{6}' = 'TUMOUR',
		# MiniSeq IDs
		'^S\\d{1,2}\\.' = 'TUMOUR.',
		# isis variant caller IDs
		'.*\\.VF$' = 'TUMOUR.AF', 
		'.*\\.GT$' = 'TUMOUR.GT', 
		# Varscan labels
		'Sample1' = 'TUMOUR',
		# isis calls allele frequency variant frequency
		# generic formatting
		'NORMAL' = paste0(variant.caller, '.NORMAL'), 
		'TUMOUR' = paste0(variant.caller, '.TUMOUR'),
		'TUMOR' = paste0(variant.caller, '.TUMOUR'),
		'.*((t|T)umour|(m|M)etastasis)(_\\d)?' = paste0(variant.caller, '.TUMOUR'), 
		'.*(G|g)ermline(_\\d)?' = paste0(variant.caller, '.NORMAL'),
		'QUAL' = paste0(variant.caller, '.QUAL'), 
		'STATUS' = paste0(variant.caller, '.STATUS')
		);

	new.column.names <- column.names;

	if( !is.null(sample.id) ) {
	  
	  sample.id <- stringr::str_replace(sample.id, '-', '.')

		regex.sample.id <- paste0('^X?', gsub('\\.', '\\\\.', sample.id), '\\.');

		new.column.names <- stringr::str_replace(
			new.column.names, 
			pattern = regex.sample.id,
			replacement = 'TUMOUR.'
			);

	}
	
	for( pattern in names(replacement.key) ) { 
		new.column.names <- stringr::str_replace(
			new.column.names, 
			pattern = pattern,
			replacement = replacement.key[ pattern ]
			);
	}

	if( 'ISIS' == variant.caller ) {
		new.column.names <- stringr::str_replace(
			new.column.names, 
			pattern = 'DEPTH',
			replacement = 'ISIS.TUMOUR.DEPTH'
			);
	}
	
	return(new.column.names);

}

