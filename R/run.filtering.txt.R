#' Run filtering on an ANNOVAR-annotated txt file
#' 
#' @param variant.file Path to variant file
#' @param caller String giving variant caller that was used (affects which filters were applied.
#' @inheritParams run.alignment.sample
run.filtering.txt <- function(
	variant.file, 
	caller = c('consensus', 'vardict', 'ides', 'mutect'), 
	output.directory = NULL,
	output.filename = NULL,
	code.directory = NULL,
	log.directory = NULL,
	config.file = NULL,
	job.dependencies = NULL,
	job.group = NULL,
	quiet = FALSE
	) { 
	
	caller <- match.arg(caller);
	
	
	### INPUT TESTS ###########################################################
	
	# TO DO:
	# 	- txt file extension
	#	- required headers

	### MAIN ##################################################################
	
	if( is.null(config.file) ) { 
		# need to save to non-temporary directory to avoid errors
		# 	=> perl calls R again, temp files already deleted
		
		config.output.directory <- NULL;
		config.output.filename <- paste0('config_filter_', caller, '_', basename(variant.file), '.yaml');
		
		if( !is.null(log.directory) ) { 
			config.output.directory <- code.directory; 
		} else if( !is.null(output.directory) ) { 
			config.output.directory <- output.directory;
		}
		
		if( !is.null(config.output.directory) ) { 
			config.output.path <- file.path(config.output.directory, config.output.filename);
		} else { 
			config.output.path <- config.output.filename;
		}
				
		config.file <- save.config(output.file = config.output.path);	
	}
	
	script <- system.file('perl', 'filter.pl', package = getOption('varitas')$pkgname );
	
	command <- make.command.line.call(
		main.command = c('perl', script), 
		options = c(
			'variant_file' = variant.file,
			'variant_caller' = caller,
			'config_file' = config.file, 
			'output_directory' = output.directory,
			'output_filename' = output.filename,
			'code_directory' = code.directory,
			'log_directory' = log.directory,
			'job_dependencies' = job.dependencies,
			'job_group' = job.group 
			)
		);
	
	if(quiet) { 
		print(command);
	} else { 
		system(command);
	}
	
}
