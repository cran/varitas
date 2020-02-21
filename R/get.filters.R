#' get.filters
#' 
#' @description
#' 	Determine filters per caller, given default and caller-specific values.
#'
#' @param filters
#' 	List of filter values. These will be updated to use default as the baseline, 
#'	with caller-specific filters taking precedence if supplied.
#'  
#' @return 
#'  A list with updated filters
#' 
#' 
get.filters <- function(filters) {

	### INPUT TESTS ###########################################################

	if( !is.list(filters) ) {
		stop('filters must be a list');
	}

	### MAIN ##################################################################

	# hardcode a minimal set of callers for now... 
	# I suspect there are better ways of doing this
	callers <- union(
		c('mutect', 'pgm', 'vardict', 'isis', 'consensus'),
		names(filters)[ 'default' != names(filters) ]
		);	

	# defaults that will be updated by user-specified filters
	baseline.filters <- list();

	if( 'default' %in% names(filters) ) {
		
		default.filters <- filters$default;

		# set baseline equal to default for each caller
		baseline.filters <- lapply(
			callers,
			function(caller, default.filters) return(default.filters),
			default.filters = default.filters
			);
		names(baseline.filters) <- callers;
	}

	new.filters <- utils::modifyList(baseline.filters, filters);

	return(new.filters);
}