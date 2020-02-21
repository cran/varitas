#' Run ANNOVAR on a VCF file
#'
#' @param vcf.file Path to VCF file
#' @param isis Logical indicating whether VCF files are from the isis (MiniSeq) variant caller
#' @inheritParams run.vardict.sample
#'
#' @return None
#'
#'
#'
#'
run.annovar.vcf <- function(
	vcf.file, 
	output.directory = NULL,
	output.filename = NULL,
	code.directory = NULL,
	log.directory = NULL,
	config.file = NULL,
	job.dependencies = NULL,
	job.group = NULL,
	job.name = NULL,
	isis = FALSE,
	quiet = FALSE,
	verify.options = !quiet
	) { 

	### INPUT TESTS ###########################################################
	
	if( 0 == length(vcf.file) ) { 
		stop('vcf.file must be provided.');
	}

	if( length(job.name) > 1 ) {
		stop('job.name must have length 1');
	}

	### MAIN ##################################################################
		
	if( verify.options ) {
		verify.varitas.options(
		  varitas.options = config.file, 
			stages.to.run = 'annotation'
			);
	}


	if(is.null(config.file)) {
		config.file <- save.config();
	}

	if(isis) {
		script <- system.file('perl', 'run_annovar_isis.pl', package = getOption('varitas')$pkgname );
	} else {
		script <- system.file('perl', 'run_annovar.pl', package = getOption('varitas')$pkgname );
	}

	command <- make.command.line.call(
		main.command = c('perl', script), 
		options = c(
			'vcf_file' = vcf.file,
			'config_file' = config.file,
			'output_directory' = output.directory,
			'output_filename' = output.filename,
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