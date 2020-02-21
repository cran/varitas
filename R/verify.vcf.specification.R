#' verify.vcf.specification
#'
#' @description
#'	Verify that VCF specification data frame fits expected format
#'
#' @param vcf.specification 
#'	VCF specification data frame
#'
#' @return None
#' 
#'
verify.vcf.specification <- function(vcf.specification) { 

	# check type
	if( !is.data.frame(vcf.specification) ) { 
		stop('vcf.specification is not a data frame.');
	}

	# check data frame dimensions – expect at least 2 columns
	if(ncol(vcf.specification) < 2) { 
		stop('vcf.specification has fewer than 2 columns.');
	}

	# check that sample.id and vcf columns exist
	if( !( 'sample.id' %in% names(vcf.specification) ) ) {
		stop('vcf.specification is missing required column sample.id');
	}

	if( !( 'vcf' %in% names(vcf.specification) ) ) {
		stop('vcf.specification is missing required column vcf');
	}   

	# check for white space in sample IDs
	if( any( grepl('\\s', vcf.specification$sample.id) ) ) { 
		stop('Sample IDs can not contain whitespace.');
	}


}


