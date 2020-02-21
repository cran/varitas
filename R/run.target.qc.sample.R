#' Get ontarget reads and run coverage quality control
#'
#' @param bam.file Path to BAM file
#' @inheritParams run.alignment.sample
#'
run.target.qc.sample <- function(
	bam.file,
	sample.id,
	output.directory = NULL,
	code.directory = NULL,
	log.directory = NULL,
	config.file = NULL,
	job.dependencies = NULL,
	job.name = NULL,
	job.group = NULL,
	quiet = FALSE
	) { 
	
	### INPUT TESTS ###########################################################
	

	
	### MAIN ##################################################################
	
	if(is.null(config.file)) {
		config.file <- save.config();
	}

	script <- system.file('perl', 'target_qc.pl', package = getOption('varitas')$pkgname );
	
	# TO DO: 
	#	- add ontarget_bam_filename option (can name it output.filename for consistency?)
	command <- make.command.line.call(
		main.command = c('perl', script), 
		options = c(
			'bam_file' = bam.file,
			'sample_id' = sample.id,
			'config_file' = config.file,
			'output_directory' = output.directory,
			'code_directory' = code.directory,
			'log_directory' = log.directory,
			'job_dependencies' = job.dependencies,
			'job_name' = job.name,
			'job_group' = job.group
			)
		);

	if(quiet) { 
	    cat(command, '\n');
	} else { 
		system(command);
	}

}

