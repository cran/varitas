#' prepare.miniseq.specifications
#' 
#' @description
#' 	Process a MiniSeq directory and sample sheet to get specification data frames
#'	that can be used to run the VariTAS pipeline.
#'
#' 	Note: This assumes normal samples are not available.
#'
#' @param sample.sheet 
#'	Data frame containing sample information, or path to a MiniSeq sample sheet
#' @param miniseq.directory
#'	Path to directory with MiniSeq files
#'
#' @return A list with specification data frames 'fastq', 'bam', and 'vcf' (as applicable)
#' @examples 
#' miniseq.sheet <- file.path(path.package('varitas'), 'extdata/miniseq/Example_template.csv')
#' miniseq.directory <- file.path(path.package('varitas'), 'extdata/miniseq')
#' miniseq.info <- prepare.miniseq.specifications(miniseq.sheet, miniseq.directory)
#' 
#'
#' @export
prepare.miniseq.specifications <- function(
	sample.sheet,
	miniseq.directory
	) {

	# TO DO:
	# 	- make a read.miniseq.sample.sheet function (more robust version of current implementation)

	### INPUT TESTS ###########################################################

	if( !dir.exists(miniseq.directory) ) {
		stop( paste(miniseq.directory, 'does not exist or is not a directory') );
	}

	### MAIN ##################################################################

	if( is.character(sample.sheet) ) {
		# check format – if starts with [Header], skip first 29 lines
		# as it means the file contains metadata, etc.
		first.line <- readLines(sample.sheet, n = 1);

		skip <- 0;
		if( '[Header]' == first.line ) skip <- 29;

		sample.sheet <- utils::read.csv(
			sample.sheet, 
			skip = skip,
			stringsAsFactors = FALSE,
			row.names = NULL
			);
	}	

	directories <- file.path(miniseq.directory, c('fastq', 'bam_bai', 'isis'));
	names(directories) <- c('fastq', 'bam', 'vcf');

	# isis could be named either VCF of isis, make sure you check both
	if( !dir.exists(directories['vcf']) ) {
		directories['vcf'] <- file.path(miniseq.directory, 'vcf');
	}

	directories <- directories[ dir.exists(directories) ];

	if( length(directories) > 0 ) {
		cat('Found directories', paste(names(directories), collapse = ' '), '\n');
	} else {
		cat('No fastq/bam/vcf directories found\n');
		return(NULL);
	}

	# use sample ID column, take everything until the first dash
	sample.ids <- sapply(
		strsplit(sample.sheet$Sample_ID, split = '-'),
		function(x) x[1]
		);

	specifications <- list();

	if( 'fastq' %in% names(directories) ) {

		r1.files <- get.miniseq.sample.files(
			sample.ids = sample.ids,
			directory = directories[ 'fastq' ],
			file.suffix = '_S\\d{1,2}_.*_R1_.*\\.fastq(\\.gz)?'
			);
		
		r2.files <- get.miniseq.sample.files(
			sample.ids = sample.ids,
			directory = directories[ 'fastq' ],
			file.suffix = '_S\\d{1,2}_.*_R2_.*(\\.gz)?'
			);

		fastq.specification <- data.frame(
			sample.id = sample.ids,
			reads = r1.files,
			mates = r2.files
			);

		verify.fastq.specification(fastq.specification);

		specifications[[ 'fastq' ]] <- fastq.specification;

	}

	if( 'bam' %in% names(directories) ) {

		tumour.bams <- get.miniseq.sample.files(
			sample.ids = sample.ids,
			directory = directories[ 'bam' ],
			file.suffix = '_S\\d{1,2}.*\\.bam'
			);

		bam.specification <- data.frame(
			sample.id = sample.ids,
			tumour.bam = tumour.bams
			);

		specifications[[ 'bam' ]] <- bam.specification;
	}	

	if( 'vcf' %in% names(directories) ) {

		vcf.files <- get.miniseq.sample.files(
			sample.ids = sample.ids,
			directory = directories[ 'vcf' ],
			file.suffix = '_S\\d{1,2}.*\\.vcf$'
			);

		vcf.specification <- data.frame(sample.id = sample.ids, vcf = vcf.files);

		specifications[[ 'vcf' ]] <- vcf.specification;
	}


	return(specifications);

}