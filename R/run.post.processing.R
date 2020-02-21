#' run.post.processing
#'
#' @description
#' 	Submit post-processing job to the cluster with appropriate job dependencies
#'
#' @param variant.specification 
#'	Data frame specifying files to be processed
#' @param output.directory
#'	Path to directory where output should be saved
#' @param code.directory 
#'	Directory where code should be saved
#' @param log.directory 
#'	Directory where log files should be saved
#' @param config.file 
#'	Path to config file
#' @param job.name.prefix
#'  Prefix for job names on the cluster
#' @param quiet 
#'	Logical indicating whether to print commands to screen rather than submit the job
#' @param email 
#'	Email address that should be notified when job finishes. If NULL or FALSE, no email is sent
#' @param verify.options 
#'	Logical indicating whether \code{verify.varitas.options()} should be run. 
#'
#' @return None
#' @examples
#' run.post.processing(
#'   variant.specification = data.frame(
#'     sample.id = c('a', 'b'),
#'     vcf = c('a.vcf', 'b.vcf'),
#'     caller = c('mutect', 'mutect'),
#'     job.dependency = c('example1', 'example2')
#'   ),
#'   output.directory = '.',
#'   quiet = TRUE
#' )
#'
#' @export
run.post.processing <- function(
	variant.specification, 	
	output.directory,
	code.directory = NULL,
	log.directory = NULL,
	config.file = NULL, 
	job.name.prefix = NULL,
	quiet = FALSE, 
	email = NULL,
	verify.options = !quiet
	) {

	### INPUT TESTS ###########################################################

	if( !is.data.frame(variant.specification) ) {
		stop('variant.specification must be a data frame.');
	}

	if( !is.character(output.directory) || length(output.directory) > 1 ) {
		stop('output.directory must be a single string');
	}

	if( !dir.exists(output.directory) && !quiet) {
		error.message <- paste('Directory', output.directory, 'does not exist or is not a directory');
		stop(error.message);
	}

	if( identical(email, FALSE) ) email <- NULL;

	if( !is.null(email) && !is.character(email) ) {
		stop('email should be a character string');
	}

	### MAIN ##################################################################

	# if job dependencies in variant.specification and no job dependencies provided,
	# consider job.dependency column to be your job dependencies – hopefully less error-prone
	if( 'job.dependency' %in% names(variant.specification)) {
		job.dependencies <- unlist(stringr::str_split(
		    variant.specification$job.dependency, 
		    pattern = ',(\\s)?|;(\\s)?|\\s'
		    ));
	}

	# save variant specification to file
	variant.specification.file <- file.path(
		output.directory, 
		paste0(Sys.Date(), '_variant_specification.txt')
		);

	if( !quiet ) {
	    utils::write.table(
	        variant.specification, 
	        variant.specification.file,
	        sep = '\t', 
	        row.names = FALSE
	        );
	}

	# config file

	# NOTE: this step is currently redundant due to no config-dependency in merging step
	# 	- include step anyways in case that changes in the future
	#
	# NOTE 2: run this before saving config to file – don't want to save wrong stuff to disk!
	
	if( verify.options ) {
	    verify.varitas.options(
	        varitas.options = config.file, 
	        stages.to.run = 'merging'
	        );
	}

	if( is.null(config.file) ) {

		config.file <- file.path(
			output.directory, 
			paste0(Sys.Date(), '_config.yaml')
			);

		if( !quiet ) save.config(output.file = config.file);
	}

	job.name <- 'post_processing';
    if( !is.null(job.name.prefix) && '' != job.name.prefix ) {
        job.name <- paste(job.name.prefix, job.name, sep = '_');
    }

	script <- system.file('perl', 'run_post_processing.pl', package = get.varitas.options('pkgname') );

	post.processing.options <-  list(
		'variant_specification' = variant.specification.file, 
		'config_file' = config.file, 
		'output_directory' = output.directory, 
		'code_directory' = code.directory, 
		'log_directory' = log.directory, 
		'job_dependencies' = job.dependencies,
		'job_name' = job.name
		);

	if( !is.null(email) ) {
		post.processing.options[[ 'email' ]] <- email;
	}

	command <- make.command.line.call(
		main.command = c('perl', script), 
		options = post.processing.options
		);

	if(quiet) { 
	    cat(command, '\n');
	} else { 
		system(command);
	}

}