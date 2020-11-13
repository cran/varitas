#' tabular.median
#'
#' @description
#'	Calculate the median of data in tabular format
#'
#' @param values 
#'	Vector of values
#' @param frequencies 
#'	Frequency corresponding to each value
#' 
#' @param ...
#'	Additional parameters passed to \code{sum}
#'
#' @return calculated median
#'
tabular.median <- function(values, frequencies, ...) {

	# TO DO:
	#	- make this more robust to huge datasets
	# 	- input tests
	
	expanded.values <- rep(values, frequencies);
	median.value <- stats::median(expanded.values, ...);
	
	return( median.value );
}

#' tabular.mean
#'
#' @description
#' 	Calculate the mean of data in tabular format
#'
#' @param values 
#'	vector of values
#' @param frequencies 
#'	frequency corresponding to each value
#' @param ...
#'	Additional parameters passed to \code{sum}
#' 
#' @return calculated mean
#'
tabular.mean <- function(values, frequencies, ...) {
	
	# TO DO: 
	#   - check this is actually robust to ... parameters
	#   - figure out why Ros had as.numeric here
	mean.value <- sum(values*frequencies, ...)/sum(frequencies, ...);
	
	return(mean.value);
}

#' Extract sample IDs from file paths
#'
#' @description
#' Extract sample IDs from a set of paths to files in sample-specific subfolders
#'
#' @param paths vector of file paths
#' @param from.filename Logical indicating whether sample ID should be extracted from filename rather than path
#'
#' @return vector of extracted sample IDs
extract.sample.ids <- function(paths, from.filename = FALSE) {
	
	# get all sample IDs
	# logic: coverage reports are in subdirectories named for each sample.
	# get second last component in file path split by /
	
	split.paths <- strsplit(paths, '/');
	
	if(from.filename) { 
		
		filenames <- sapply(
			split.paths,
			function(components) {
				return( components[ length(components) ] );
			});
			
			
		filename.components <- stringr::str_split(
			filenames,
			pattern = '\\.', 
			n = 2
			);
			
		sample.ids <- sapply(
			filename.components,
			function(components) return( components[1] )
			);	
	
		
	} else { 
		sample.ids <- sapply(
			split.paths,
			function(components) {
				return( components[ length(components) - 1 ] );
			});
	
	}

		
	return(sample.ids);	
	
}

#' Run ls command
#'
#' @description
#' Runs ls command on system. This is a workaround since list.files can not match patterns based on subdirectory structure.
#'
#' @param pattern pattern to match files
#' @param directory base directory command should be run from
#' @param error logical indicating whether to throw an error if no matching founds found. Defaults to False.
#' 
#' @return paths returned by ls command
system.ls <- function(pattern = "", directory = "", error = FALSE) {
	
	### INPUT TESTS
	
	
		
	### MAIN	
		
	# if a directory has been passed, make sure it ends with "/"
	# no need to check if it already does, as "//" will not do any harm
	if("" != directory) {
		directory <- paste0(directory, "/");		
	}
		
	system.command <- paste0("ls ", directory, pattern);
	paths <- try(system(system.command, intern = TRUE), silent = TRUE);
	
	if(error && 0 == length(paths)) {
		error.message <- paste("Did not find any files matching", pattern);	
		if("" != directory) error.message <- paste(error.message, "in", directory);
			
		stop(error.message);
	}
		
	return(paths);
}


#' Get pool corresponding to each amplicon
#' 
#' @description
#' The bed files are not consistent, so it's not clear where the pool will appear. 
#' This function parses through the columns to identify where the pool 
#'
#' @param panel.data data frame pool should be extracted from
#'
#' @return pools vector of pool information
get.pool.from.panel.data <- function(panel.data) {

	# first three columns are chr, start, end (per bed file requirements)
	#	 => parse through remaining columns from right to left to see if they contain pool information

	pool.column <- NULL;

	for(i in seq( ncol(panel.data), 4 )) { 

		contains.pool <- grepl("pool", panel.data[, i], ignore.case = TRUE);

		if( all(contains.pool) ) {
			pool.column <- i;
			break;
		}
	}

	if( is.null(pool.column) ) { 
		warning("Unable to find pool information");
		return(NULL);
	}

	pools <- stringr::str_extract(
		panel.data[, pool.column], 
		"(P|p)ool(\\s)?=(\\s)?(.*)(;)?"
		);

	return(pools);

}

#' Summarise panel coverage by gene
#'
#' @param panel.file path to panel
#' @param gene.col index of column containing gene name
#' 
#' @return panel.coverage.by.gene data frame giving the number of amplicons and their total length by gene
get.panel.coverage.by.gene <- function(panel.file, gene.col = 5) { 

	# TO DO:
	#	- warn if content of gene.col does not look like a gene
	# 	- standardize whether functions take data frames or file paths
	
	panel.data <- utils::read.delim(
		panel.file, 
		sep = "\t",
		header = FALSE
		);	

	panel.data$length <- panel.data[, 3] - panel.data[, 2];

	# get gene without exon specification
	# column gives gene exX if specified at exon level -> keep anything before first space
	panel.data$gene <- stringr::str_extract(
		panel.data[, gene.col], 
		"[^\\s]+"
		);

	panel.coverage.by.gene <- panel.data %>% 
		dplyr::group_by(gene) %>%
		dplyr::summarise(
			n.amplicons = n(), 
			total.length = sum(length)
			);

	return( data.frame(panel.coverage.by.gene) );
}



#' Generate a colour scheme
#' 
#' @param n Number of colours desired
#'  
#' @return Colour.scheme generated colours
get.colours <- function(n) {

	full.colour.scheme <- c(
		'#1f78b4',
		'#e31a1c',
		'#ff7f00',
		'#6a3d9a',
		'#b15928',
		'#33a02c',
		# second half of pairs
		'#a6cee3',
		'#fb9a99',
		'#fdbf6f',
		'#cab2d6',
		'#ffff99',
		'#b2df8a'
		);

	colour.scheme <- full.colour.scheme[1:n];

	return( colour.scheme );
}

#' Make string with command line call from its individual components
#'
#' @param main.command String or vector of strings giving main part of command (e.g. "python test.py" or c("python", "test.py"))
#' @param options Named vector or list giving options
#' @param flags Vector giving flags to include.
#' @param option.prefix String to preface all options. Defaults to "--"
#' @param option.separator String to separate options form their values. Defaults to a single space.
#' @param flag.prefix String to preface all flags. Defaults to "--"
#' 
#' @return command string giving command line call
#' 
#'
make.command.line.call <- function(
		main.command, 
		options = NULL, 
		flags = NULL,
		option.prefix = "--", 
		option.separator = " ",
		flag.prefix = "--"
		) { 

	### INPUT TESTS
	if( !is.character(main.command) ) {
		stop("Argument main.command must be a string or character vector.");
	}

	### MAIN 

	# TO DO:
	#	- add support for flags 
	#	- add tests
	# 	- add input tests
	
	# Remove any elements given as NULL or empty vectors. 
	# Adding this removal step was easier than making sure the function would never be passed any empty options (for job dependencies)
	options <- options[!sapply(options, is.null) & sapply(options, length) > 0];

	# if options or main command have been passed as vector, parse into string
	main.command <- paste(main.command, collapse = " ");
	

	# merge options into a single string
	combined.options <- '';

	if( !is.null(options) ) { 
		options <- lapply(options, paste, collapse = " ");

		option.strings <- paste0(
			option.prefix, 
			names(options), 
			option.separator,
			options
			);

		combined.options <- paste(option.strings, collapse = " ");

	}

	# merge flags into a single string
	combined.flags <- '';

	if( !is.null(flags) ) {
		flag.strings <- paste0(
			flag.prefix, 
			flags
			);

		combined.flags <- paste(flag.strings, collapse = " ");
	}


	# combine main command and options
	# if no options have been passed in, will the extra space cause problems? 
	command <- paste(main.command, combined.options, combined.flags);

	return(command);
}

#' Parse job dependencies
#'
#' @description
#' Parse job dependencies to make the functions more robust to alternate inputs (e.g. people writing alignment instead of bwa)
#'
#' @param dependencies Job dependency strings to be parsed. 
#'
#' @return parsed.dependencies Vector of job dependencies after reformatting.
parse.job.dependencies <- function(dependencies) { 

	### INPUT TESTS

	if( !is.character(dependencies) ) { 
		stop('The argument dependencies must be a string or vector of strings.');
	}

	### MAIN

	# convert to lower case
	parsed.dependencies <- tolower(dependencies);
	parsed.dependencies[parsed.dependencies %in% c('bwa', 'alignment', 'bwa mem', 'bwamem', 'align')] <- 'bwa';


	return(parsed.dependencies);

}

#' save.config
#'
#' @description
#' 	Save current varitas config options to a temporary file, and return filename.
#'
#' @param output.file Path to output file. If NULL (default), the config file will be saved as a temporary file.
#'
#' @return Path to config file
#'
save.config <- function(output.file = NULL) {

	config <- getOption('varitas');

	if( is.null(output.file) ) { 
		config.file <- tempfile(fileext = '.yaml');
	} else { 
		config.file <- output.file;	
	} 
	
	# Note: writeLines does not support indentation
	#  	- use cat and capture.output instead
	utils::capture.output(
		cat(yaml::as.yaml(config, indent.mapping.sequence = TRUE)),
		file = config.file
		);

	return(config.file);
}

#' get.file.path 
#'
#' @description
#' 	Get absolute path to sample-specific file for one or more samples
#'
#' @param sample.ids 
#'	Vector of sample IDs to match filename on
#' @param directory 
#'	Path to directory containing files
#' @param extension 
#'	String giving extension of file
#' @param allow.multiple 
#'	Boolean indicating whether to allow multiple matching files. 
#'	Defaults to false, which throws an error if the query matches more than one file.
#' @param allow.none 
#'	Boolean indicating whether to allow no matching files. 
#'	Defaults to false, which throws an error if the query does not match any files.
#'
#' @return Paths to matched files
#'
get.file.path <- function(
		sample.ids, 
		directory, 
		extension = NULL,
		allow.multiple = FALSE,
		allow.none = FALSE
		) {
		
	### INPUT TESTS ###########################################################
		
	# avoid problems with identifying which sample is from 
	if(allow.multiple && length(sample.ids) > 1) { 
		stop('Cannot match multiple filepaths for more than one sample.'); 
	}
	
	### MAIN ##################################################################
	
	paths <- c();
	
	for(sample.id in sample.ids) { 
		if( !is.null(extension) ) { 
			extension.pattern <- paste0('.*\\.', extension, '$');
			pattern <- paste0(sample.id, extension.pattern);
		} else { 
			pattern <- sample.id;
		}
		
		filename <- list.files(
			pattern = pattern, 
			path = directory
			);	
			
		# error checks	
		if( 0 == length(filename)) {
			
			if(allow.none) {
				return(NULL);
			} else {
				error.message <- paste('No file found for', sample.id, 'in directory', directory);
				stop(error.message);
			}
		}
	
		if( !allow.multiple && length(filename) > 1) {
			error.message <- paste(
				'Found more than one file for', 
				sample.id, 
				'in directory', 
				directory, 
				'\n\n',
				paste(filename, collapse = '\n')
				);

			stop(error.message);
		}
	
		paths <- c(paths, file.path(directory, filename));
	}
		
	return(paths);
}



#' create.directories
#'
#' @description
#' 	Create directories in a given path
#' 
#' @param directory.names 
#' 	Vector of names of directories to be created
#' @param path 
#'	Path where directories should be created
#'  
create.directories <- function(directory.names, path) { 

	### INPUT TESTS

	if( !is.character(directory.names) ) { 
		stop('directory.names must be a character vector.');
	}

	if( !is.character(path) && 1 == length(path) ) { 
		stop('path must be a single string.');
	}

	### MAIN

	# loop over directories to create, and create them
	for( directory.name in directory.names ) { 

		directory.path <- file.path(path, directory.name);
		if( !dir.exists(directory.path) ) { 
			dir.create(directory.path);
		}
		
	}

}

 
#' get.buildver 
#'
#' @description
#' 	Get build version (hg19/hg38) based on settings.
#'
#' @description
#' Parses VariTAS pipeline settings to get the build version. When this function was first developed, the idea was to
#' be able to explicitly set ANNOVAR filenames based on the build version.
#'
#' @return String giving reference genome build version (hg19 or hg38)
#'
get.buildver <- function() { 
	
	varitas.options <- get.varitas.options();
	reference.build <- tolower( varitas.options$reference_build );
	
	if( 'grch37' == reference.build ) {
		buildver <- 'hg19';
	} else if( 'grch38' == reference.build ) { 
		buildver <- 'hg38';
	} else {
		# unrecognized reference build, throw an error
		error.message <- paste(reference.build, 'is not a recognized reference build. Please update the VariTAS settings');		
		stop(error.message);
	}
	
	return(buildver);
}


#' classify.variant
#'
#' @description
#' 	Classify a variant as SNV, MNV, or indel based on the reference and alternative alleles
#'
#' @param ref 
#'	Vector of reference bases
#' @param alt 
#'	Vector of alternate bases
#'
#' @return Character vector giving type of variant.
#'
classify.variant <- function(ref, alt) { 
	
	### INPUT TESTS ###########################################################
	
	if( length(ref) != length(alt) ) { 
		stop('ref and alt vectors must be the same length');
	}
	
	### MAIN ##################################################################
	
	is.indel <- nchar(ref) != nchar(alt);
	is.mnv <- (nchar(ref) == nchar(alt)) & (nchar(ref) > 1);
	is.snv <- (1 == nchar(ref)) & (1 == nchar(alt));

	
	# quality control checks
	if( any(is.snv & is.indel) | any(is.mnv & is.indel) | any(is.snv & is.mnv) ) { 
		stop('At least one variant has been classified in two or more categories.');
	}
	
	
	variant.types <- rep(NA, length(ref));
	variant.types[is.indel] <- 'indel';
	variant.types[is.snv] <- 'SNV';
	variant.types[is.mnv] <- 'MNV';


	return(variant.types);
}

#' mean.field.value
#'
#' @description
#'	Get mean value of a variant annotation field
#'
#'
#' @details
#' 	As part of the variant merging process, annotated variant data frames are merged into 
#' 	one, with the value from each caller prefixed by CALLER. For example, the VarDict normal allele 
#' 	freqeuncy will have header VARDICT.NORMAL.AF. This function takes the average of all callers' value 
#' 	for a given field, removing NA's. If only a single caller is present in the data frame, that value is returned.
#'
#' @param variants Data frame with variants
#' @param field String giving field of interest. 
#' @param caller String giving caller to calculate values from
#'
#'
#'
#' @return Vector of mean values.
mean.field.value <- function(
	variants, 
	field = c('TUMOUR.DP', 'NORMAL.DP', 'NORMAL.AF', 'TUMOUR.AF', 'QUAL'), 
	caller = c('consensus', 'vardict', 'pgm', 'mutect', 'isis', 'varscan', 'lofreq')
	) {
	
	# TO DO:
	#	- add support for taking average of (say) two out of four variants? 
	
	field <- match.arg(field);
	caller <- match.arg(caller);
	
	### INPUT TESTS ###########################################################
	
	if( 'QUAL' == field && sum( grepl('QUAL$', names(variants)) ) > 1 && 'consensus' == caller ) { 
		warning('Taking the average of several QUAL field. Probably not a good idea as they tend to be wildly different');
	}
	
	### MAIN ##################################################################
	
	if( 'consensus' == caller ) { 
		# add end of line suffix $ to minimize risk of errors
		values <- rowMeans(
			variants[ , grepl(paste0(field, '$'), names(variants))], 
			na.rm = TRUE
			);
		
		# if no non-NA values for a field, NaN is returned. Change back to NA to avoid confusion
		values[ is.nan(values) ] <- NA;
	
	} else { 
	
		caller.prefix <- '';
		if( any( grepl(toupper(caller), names(variants)) ) ) { 
			caller.prefix <- paste0(toupper(caller), '.');
		}
		
		values <- variants[, paste0(caller.prefix, field)];
	
	}	
	
	return(values);

}

#' Get base substitution
#'
#' @description
#' 	Get base substitution represented by pyrimidine in base pair. 
#'	If more than one base in REF/ALT (i.e. MNV or indel rather than SNV), NA will be returned
#'
#' @param ref Vector of reference bases
#' @param alt Vector of alternate bases
#'
#' @return base.substitutions
get.base.substitution <- function(ref, alt) {
	
	variant.key <- c(
    	'A>C' = 'T>G', 'A>G' = 'T>C', 'A>T' = 'T>A',
    	'C>A' = 'C>A', 'C>G' = 'C>G', 'C>T' = 'C>T',
    	'T>G' = 'T>G', 'T>C' = 'T>C', 'T>A' = 'T>A',
    	'G>T' = 'C>A', 'G>C' = 'C>G', 'G>A' = 'C>T'
    	);

	base.substitutions <- variant.key[ paste(ref, alt, sep = '>') ];
	
	return(base.substitutions);

}

#' build.variant.specification
#' 
#' @description
#' 	Build data frame with paths to variant files.
#' 
#' @details
#'	Parses through sample IDs in a project directory and returns paths to variant files based on
#' 	(theoretical) file name patterns. Useful for testing, or for entering the pipeline at non-traditional stages.
#'
#' @param sample.ids Vector of sample IDs. Must match subdirectories in project.directory.
#' @param project.directory Path to directory where sample subdirectories
#'
#' @return Data frame with paths to variant files.
#'
#'
build.variant.specification <- function(sample.ids, project.directory) { 

	# TO DO:
	#	- support for grch38

	### INPUT TESTS ###########################################################

	if( !dir.exists(project.directory) ) { 
		error.message <- paste('Directory', project.directory, 'does not exist');
		stop(error.message);
	}

	### MAIN #################################################################

	caller.file.suffixes <- c(
		'vardict' = '.passed.ontarget.vcf.annovar.hg19_multianno.vcf.txt', 
		'mutect' = '.passed.ontarget.vcf.annovar.hg19_multianno.vcf.txt', 
		'pgm' = '.snvs_and_cnvs.vcf.annovar.hg19_multianno.vcf.txt'
		);

	# store paths to variant files we come across along the way
	variant.specification <- list();

	for(sample.id in sample.ids) { 
		
		sample.directory <- file.path(project.directory, sample.id);

		if( !dir.exists(sample.directory) ) {
			error.message <- paste('No directory found for sample', sample.id);
			stop(error.message);
		}

		# loop over variant callers to see if file exists
		for(variant.caller in names(caller.file.suffixes) ) { 

			variant.file <- file.path(
				sample.directory, 
				variant.caller, 
				paste0(sample.id, caller.file.suffixes[ variant.caller ])
				);

			if( file.exists(variant.file) ) {

				variant.specification[[ paste0(sample.id, variant.caller) ]] <- data.frame(
					sample.id = sample.id, 
					variant.file = variant.file,
					caller = variant.caller
					);

			}

		} # end of caller loop
	} # end of sample loop

	# make data frame and get rid of uninformative row names
	variant.specification <- do.call(rbind, variant.specification);
	rownames(variant.specification) <- NULL; 

	return(variant.specification);
}


#' capitalize.caller 
#'
#' @description
#'	Capitalize variant caller name
#'
#' @param caller 
#'	Character vector of callers to be capitalized
#'
#' @return Vector of same length as caller where eligible callers have been capitalized
#'
#'
capitalize.caller <- function(caller) { 

	capitalization.key <- c(
		'mutect' = 'MuTect', 
		'vardict' = 'VarDict', 
		'pgm' = 'PGM',
		'isis' = 'Isis'
		);


	capitalized.caller <- capitalization.key[ caller ];
	capitalized.caller[ is.na(capitalized.caller) ] <- caller[is.na(capitalized.caller)]

	return(capitalized.caller);
}

#' @rdname capitalize.caller
#'
#'
capitalise.caller <- capitalize.caller;


#' split.on.column
#' 
#' @description 
#'  Split data frame on a concatenated column.
#'
#' @param dat Data frame to be processed
#' @param column Name of column to split on
#' @param split.character Pattern giving character to split column on
#'
#' @return Data frame after splitting on column
split.on.column <- function(dat, column, split.character) {

	### INPUT TESTS ###########################################################

	if( !is.character(column) ) {
		stop('column must be a character corresponding to the column header');
	}

	if( length(column) > 1) {
		stop('column cannot be a vector.');
	} 

	if( !(column %in% names(dat)) ) {
		error.message <- paste(column, 'not found in data frame');
		stop(error.message);
	}

	### MAIN ##################################################################

	column.components <- stringr::str_split(
		dat[, column], 
		pattern = split.character
		);


	# variable for processed output
	new.dat <- list();

	for(i in 1:nrow(dat) ) {	

		components <- column.components[[i]];

		for( j in 1:length(components) ) {

			data.row <- dat[i, names(dat) != column];
			data.row[, column] <- components[j];

			# add to list with a unique identifier to prevent overwriting previous entries
			new.dat[[ paste(i, j, sep = '-')]] <- data.row;

			}
		}

	# convert to data frame
	new.dat <- do.call(rbind, new.dat);

	# remove nonsensical row names
	rownames(new.dat) <- NULL;

	# reorder columns
	new.dat <- new.dat[, names(dat)];

	# sanity check
	if( length(unlist(column.components)) != nrow(new.dat) ) {
		stop('Unexpected number of rows in processed data frame');
	}

	return(new.dat);
}


#' get.fasta.chromosomes
#'
#' @description 
#' 	Extract chromosomes from fasta headers.
#'
#' @param fasta Path to reference fasta
#'
#' @return Vector containing all chromosomes in fasta file.
#'
#'
get.fasta.chromosomes <- function(fasta) {

	# grep for headers
	header.file <- tempfile();
	header.grep.command <- paste(
		"grep '>'", fasta
		);

	headers <- system(header.grep.command, intern = TRUE);
	chromosomes <- stringr::str_extract(
		headers, 
		pattern = '(?<=>)[^\\s]+' # match everything between > and first whitespace
		);

	return(chromosomes);
}

#' get.vcf.chromosomes
#'
#' @description
#' 	Extract chromosomes from a VCF file.
#'
#' @param vcf Path to VCF file
#' 
#' @return Vector containing all chromosomes in VCF
#'
#'
get.vcf.chromosomes <- function(vcf) { 

	chromosome.command <- paste(
		"grep -v '^#'", 
		vcf, 
		"| awk -F '\t' '{print $1}' | uniq"
		);

	chromosomes <- system(chromosome.command, intern = TRUE);

	return(chromosomes);
}

#' get.bed.chromosomes
#' 
#' @description
#' 	Extract chromosomes from bed file
#'
#' @param bed Path to BED file
#'
#' @return Vector containing all chromosomes in BED file
#'
#'
get.bed.chromosomes <- function(bed) {
	
	# TO DO: allow for header

	# NOTE: 
	# - bed file does not have to be sorted for our applications, sort here to get unique values
	chromosome.command <- paste(
		"awk -F '\t' '{print $1}'", 
		bed,
		"| sort | uniq"
		);

	chromosomes <- system(chromosome.command, intern = TRUE);

	# remove anything starting with "track name" – header
	chromosomes <- chromosomes[ !grepl('^track', chromosomes) ];

	return(chromosomes);
}


#' logical.to.character
#' 
#' @description 
#'	Convert a logical vector to a T/F coded character vector. Useful for preventing unwanted T->TRUE nucleotide conversions
#'
#' @param x Vector to be converted
#'
#' @return Character vector after converting TRUE/FALSE
#'
#'
#'
logical.to.character <- function(x) {

	character.x <- as.character(x);

	# only convert if vector was logical to start with
	# people could have strings 'TRUE' and 'FALSE' and want to keep them (ugh)
	if( is.logical(x) ) {
		character.x[ 'TRUE' == character.x ] <- 'T';
		character.x[ 'FALSE' == character.x ] <- 'F';
	}

	return(character.x);
}



#' date.stamp.file.name
#'
#' @description
#'  Prefix file name with a date-stamp.
#'
#' @param file.name File name to be date-stamped
#' @param date Date to be added. Defaults to current date.
#' @param separator String that should separate the date from the file name. Defaults to a single underscore.
#'
#' @return String giving the datestamped file name
#'
#' @examples
#' date.stamp.file.name('plot.png');
#' date.stamp.file.name('yesterdays_plot.png', date = Sys.Date() - 1);
#'
#' @aliases datestamp.file.name datestamp.filename
#'
#' @export
date.stamp.file.name <- datestamp.file.name <- datestamp.filename <- function(
    file.name,
    date = Sys.Date(),
    separator = '_'
    ) {

    ### INPUT TESTS ###########################################################

    if( !is.character(file.name) ) stop('file.name must be a string');

    if( grepl('/', file.name) ) stop('Detected forward slash in file.name. Unable to datestamp directories.');

    if(grepl('\\s', file.name) ) warning('Your file name contains whitespace - are you sure you want to use it?');

    ### MAIN ##################################################################

    datestamped.file.name <- paste(
        date,
        file.name,
        sep = separator
        );

    return(datestamped.file.name);
}

#' read.yaml
#'
#' @description Read a yaml file
#' 
#' @param file.name Path to yaml file
#' 
#' @return list containing contents of yaml file 
#' 
#' @examples
#' read.yaml(file.path(path.package('varitas'), 'config.yaml'))
#' 
#' @export
#'
read.yaml <- function(file.name) {
  
  # use calling handlers to avoid warning if last line isn't blank
  # function is perfectly able to handle it, but throws a warning
  contents <- withCallingHandlers( 
    yaml::yaml.load_file(file.name),
    warning = function(w) { 
      if( grepl('incomplete final line', w$message) ) invokeRestart('muffleWarning');
    }
  );
  
  return(contents);
}

#' alternate.gene.sort
#'
#' @description 
#'  Given a data frame containing coverage statistics and gene information, returns that frame
#'  with the rows sorted by alternating gene size (for plotting)
#'  
#' @details 
#'  Genes have varying numbers of associated amplicons and when plotting coverage statistics, 
#'  if two genes with very low numbers of amplicons are next to each other, the labels will overlap.
#'  This function sorts the coverage statistics data frame in a way that places the genes
#'  with the most amplicons (largest) next to those with the least (smallest).
#' 
#' @param coverage.statistics Data frame of coverage statistics
#' 
#' @return Coverage statistics data frame sorted by alternating gene size 
#' 
#'
#' 
alternate.gene.sort <- function(coverage.statistics) {
  genes <- get.gene(coverage.statistics);
  
  gene.start <- vapply(
    unique(genes), 
    function(x, genes) match(x, genes), 
    genes = genes,
    FUN.VALUE = 0
  );
  
  gene.end <- c(gene.start[-1], nrow(coverage.statistics) + 1);
  
  lengths <- c()
  for (i in 1:length(gene.start)) {
    lengths <- c(lengths, gene.end[i] - gene.start[i])
  }
  
  combined <- data.frame(gene=unique(genes), length=lengths, stringsAsFactors = FALSE)
  
  length.order <- order(lengths, decreasing = TRUE)
  combined <- combined[ length.order, ]
  
  total.genes <- nrow(combined)
  alternate <- data.frame()
  alternate <- rbind(alternate, combined[total.genes, ])
  big <- TRUE
  big.iter <- 1
  small.iter <- total.genes - 1
  for (i in 1:(total.genes-1)) {
    if (big) {
      alternate <- rbind(alternate, combined[big.iter, ])
      big.iter <- big.iter + 1
      big <- FALSE
    } else {
      alternate <- rbind(alternate, combined[small.iter, ])
      small.iter <- small.iter - 1
      big <- TRUE
    }
  }
  
  new.frame <- data.frame(stringsAsFactors = FALSE)
  old.frame <- coverage.statistics
  
  for (i in 1:total.genes) {
    sort.gene <- alternate[i,1]
    for (j in 1:nrow(old.frame)) {
      test.gene <- genes[j]
      if (sort.gene == test.gene){
        new.frame <- rbind(new.frame, old.frame[j, ])
      }
    }
  }
  
  return(new.frame)
}

#' fix.varscan.af
#'
#' @description 
#'  VarScan does not output allele frequencies, so this script calculates them from the
#'  DP (depth) and AD (variant allele depth) values and adds them to the annotated vcf. 
#' 
#' @param variant.specification Data frame of variant file information
#' 
#'
#' 
fix.varscan.af <- function(variant.specification) {
  for (i in 1:nrow(variant.specification)) {
    if (variant.specification[i, "caller"] != "varscan") {
      next
    }
    variant.file <- variant.specification[i, "variant.file"]
    output.dir <- dirname(variant.file)
    sample.id <- variant.specification[i, "sample.id"]
    vcf.file <- file.path(output.dir, paste0(sample.id, '.passed.ontarget.vcf'))
    vcf.somatic.line <- readLines(vcf.file)[4]
    if (grepl('SOMATIC', vcf.somatic.line, fixed=TRUE)) {
      somatic <- TRUE
    } else {
      somatic <- FALSE
    }
    vcf.df <- try(utils::read.table(vcf.file, stringsAsFactors = FALSE, header = FALSE), silent = TRUE)
    if (inherits(vcf.df, 'try-error')) { 
      next
    }
    variant.df <- utils::read.table(variant.file, stringsAsFactors = FALSE, header = TRUE)
    if (nrow(variant.df) < 1) {
      next
    }
    try(variant.df[,"NORMAL.AF"] <- variant.df[,"NORMAL.AD"] / variant.df[,"NORMAL.DP"], silent = TRUE)
    try(variant.df[,"TUMOR.AF"] <- variant.df[,"TUMOR.AD"] / variant.df[,"TUMOR.DP"], silent = TRUE)
    try(variant.df[,"TUMOR.AF"] <- variant.df[,"Sample1.AD"] / variant.df[,"Sample1.DP"], silent = TRUE)
    if (somatic) {
      variant.df <- variant.df[which(grepl('SOMATIC', vcf.df$V8, fixed=TRUE)), ]
    }
    utils::write.table(variant.df, file = variant.file, sep = '\t', quote = FALSE)
  }
}

#' fix.lofreq.af
#'
#' @description 
#'  LoFreq also does not output allele frequencies, so this script calculates them from the
#'  DP (depth) and AD (variant allele depth) values--which are also not output nicely--
#'  and adds them to the annotated vcf. 
#' 
#' @param variant.specification Data frame of variant file information
#' 
#'
#' 
fix.lofreq.af <- function(variant.specification) {
  variant.spec <- variant.specification[variant.specification["caller"] == 'lofreq', ]
  annotated.files <- variant.spec["variant.file"]
  if ( nrow(variant.spec) == 1 ) {
    single.file <- TRUE
  } else {
    single.file <- FALSE
  }
  for (i in 1:nrow(variant.spec)) {
    sample.id <- variant.spec[i, "sample.id"]
    if ( single.file ) {
      output.dir <- dirname(as.character(annotated.files[i]))
    } else {
      output.dir <- dirname(annotated.files[i, "variant.file"])
    }
    vcf.file <- file.path(output.dir, paste0(sample.id, '.passed.ontarget.vcf'))
    variant.df <- try(utils::read.table(vcf.file, stringsAsFactors = FALSE, header = FALSE), silent = TRUE)
    if (inherits(variant.df, 'try-error')) { 
      next
    }
    DP <- stringr::str_extract(variant.df[,8], "(?<=DP=)\\d+(?=;)")
    AF <- stringr::str_extract(variant.df[,8], "(?<=AF=)\\d\\.\\d+(?=;)")
    AD.str <- stringr::str_extract(variant.df[,8], "(?<=DP4=)\\d+,\\d+,\\d+,\\d+")
    AD <- sapply(AD.str, sum.dp4)
    
    if ( single.file ) {
      annovar.table <- utils::read.table(as.character(annotated.files[i]), stringsAsFactors = FALSE, header = TRUE)
    } else {
      annovar.table <- utils::read.table(annotated.files[i, "variant.file"], stringsAsFactors = FALSE, header = TRUE)
    }
    try(annovar.table["TUMOUR.DP"] <- DP, silent = TRUE)
    try(annovar.table["TUMOUR.AF"] <- AF, silent = TRUE)
    try(annovar.table["TUMOUR.AD"] <- AD, silent = TRUE)
    # Remove duplicated variants (happens when BED regions overlap)
    annovar.table <- annovar.table[!duplicated(annovar.table[,1:4]),]
    if ( single.file ) {
      utils::write.table(annovar.table, file = as.character(annotated.files[i]), sep = '\t', quote = FALSE)
    } else {
      utils::write.table(annovar.table, file = annotated.files[i, "variant.file"], sep = '\t', quote = FALSE)
    }
  }
}

#' sum.dp4
#'
#' @description 
#'  Simply calculates the depth of coverage of the variant allele given a string of
#'  DP4 values 
#' 
#' @param dp4.str String of DP4 values in the form "1234,1234,1234,1234"
#' 
#'
#' 
sum.dp4 <- function(dp4.str) {
  dp4.nums <- strsplit(dp4.str, ',')
  dp4.nums <- as.integer(dp4.nums[[1]])
  return(dp4.nums[3] + dp4.nums[4])
}
