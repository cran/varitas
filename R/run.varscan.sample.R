#' Run VarScan for a sample
#' 
#' @inheritParams run.vardict.sample
#'
#'
#'
#'
run.varscan.sample <- function(
	tumour.bam,
	sample.id,
	paired,
	normal.bam = NULL, 
	output.directory = NULL,
	output.filename = NULL,
	code.directory = NULL,
	log.directory = NULL,
	config.file = NULL,
	job.dependencies = NULL, 
	quiet = FALSE,
	job.name = NULL,
	verify.options = !quiet,
	job.group = NULL
	) { 
	
	### INPUT TESTS ###########################################################

	# make sure no factors have been passed in
	assertthat::assert_that( is.character(tumour.bam) );
	assertthat::assert_that( is.character(sample.id) );
	assertthat::assert_that( is.null(normal.bam) || is.character(normal.bam));
	assertthat::assert_that( is.null(output.directory) || is.character(output.directory) );


	if(paired && is.null(normal.bam)) {
		stop('paired is set to true but no normal sample BAM has been supplied.');
	}	
	
	if(!paired && !is.null(normal.bam)) {
		stop('Getting mixed signals: paired is set to false but a normal sample BAM has been supplied.');
	}
	
	### MAIN ##################################################################

	if( verify.options ) {
		verify.varitas.options(
			varitas.options = config.file, 
			stages.to.run = 'calling', 
			variant.callers = 'varscan'
			);
	}


	# save temporary config file that can be read by Perl
	if(is.null(config.file)) {
		config.file <- save.config();
	}


	# sort out whether to pass paired flag to Perl
	flags <- NULL;

	if(paired) {
		flags <- 'paired';
	}


	script <- system.file('perl', 'run_varscan.pl', package = get.varitas.options('pkgname') );
			
	command <- make.command.line.call(
		main.command = c('perl', script),
		options = c(
			'tumour_bam' = tumour.bam,
			'normal_bam' = normal.bam,
			'sample_id' = sample.id,
			'config_file' = config.file,
			'output_directory' = output.directory,
			'output_filename' = output.filename,
			'code_directory' = code.directory,
			'log_directory' = log.directory,
			'job_dependencies' = job.dependencies,
			'job_name' = job.name,
			'job_group' = job.group
			), 
		flags = flags
		);
	

	if( quiet ) { 
	    cat(command, '\n');
	} else {
		system(command);
	}
			
}


