#' Filter variant calls
#'
#' @description
#' Filter data frame of variant calls based on thresholds specified in settings.
#'
#' @param variants Data frame of variant calls with ANNOVAR annotation, or path to variant file.
#' @param caller Name of caller used (needed to match appropriate filters from settings)
#' @param config.file Path to config file to be used. If not supplied, will use the pre-existing VariTAS options.
#' @param verbose Logical indicating whether to output descriptions of filtering steps. Defaults to False, useful for debugging.
#'
#' @return filtered.variants Data frame of filtered variants
#'
#'
filter.variants <- function(
    variants, 
    caller = c('vardict', 'ides', 'mutect', 'pgm', 'consensus', 'isis', 'varscan', 'lofreq'), 
    config.file = NULL,
    verbose = FALSE
    ) { 
    
    ## TO DO:
    #	- log file outputting how many variants were filtered? 
    
    ### INPUT TESTS ###########################################################
    
    caller <- match.arg(caller);
    filters <- get.varitas.options(paste0('filters.', caller));
    
    if( !in.varitas.options( paste0('filters.', caller) ) ) { 
        error.message <- paste('No filters found in VariTAS settings for variant caller', caller);
        stop(error.message);
    } 	
    
    if( !is.data.frame(variants) ) { 
        stop('variants must be a data frame');
    }
    
    ### MAIN ##################################################################	

    print(filters);
    
    # if a config file has been passed in, overwrite config options
    # this is useful if we call R from perl.
    if( !is.null(config.file) ) {  
        config <- yaml::yaml.load_file(config.file);
        config$pkgname <- get.varitas.options('pkgname');
        
        options(varitas = config);
    }
    
    # if American spelling in variant file header, change to British one
    # Do this before saving names to a variable – American spelling should be eradicated and we do not
    # want to change back!
    names(variants) <- stringr::str_replace(
        names(variants), 
        pattern = 'TUMOR', 
        replacement = 'TUMOUR'
    );
    
    
    ## FILTERING
    
    old.names <- names(variants);
    
    # if .DEPTH instead of .DP, fix
    names(variants) <- gsub('.DEPTH', '.DP', names(variants));		
    
    # detect if paired analysis (i.e. normal provided in variants)
    paired <- FALSE;
    if( any(grepl('NORMAL', names(variants))) ) { 
        paired <- TRUE;
    }
    
    filtered.variants <- variants;
    
    if(verbose) { 
        cat('\nApplying filters to', nrow(variants), 'calls from', caller, '\n');
    }
    
    
    # NOTES: 
    #	1) Want to retain rows with NA in the column we are filtering on. There is no easy way to do this, so code it explicitly.
    # 	2) We are dropping variants at each step => need to call mean.field.value repeatedly. If not for consensus option, we could just use subset function 
    
    # reads suppporting variant
    if( 'min_tumour_variant_reads' %in% names(filters) ) {
        
        tumour.af <- mean.field.value(filtered.variants, field = 'TUMOUR.AF', caller = caller);
        tumour.dp <- mean.field.value(filtered.variants, field = 'TUMOUR.DP', caller = caller); 
        
        passed.filter <- (tumour.dp*tumour.af >= filters$min_tumour_variant_reads) | is.na(tumour.dp*tumour.af);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied min_tumour_variant_reads filter, and removed', sum(!passed.filter), 'variants\n');
        }
    }
    
    if( paired && 'max_normal_variant_reads' %in% names(filters) && caller != 'lofreq' ) {
        normal.af <- mean.field.value(filtered.variants, field = 'NORMAL.AF', caller = caller);
        normal.dp <- mean.field.value(filtered.variants, field = 'NORMAL.DP', caller = caller); 
        
        passed.filter <- (normal.dp*normal.af <= filters$max_normal_variant_reads) | is.na(normal.dp*normal.af);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied max_normal_variant_reads filter, and removed', sum(!passed.filter), 'variants\n');
        }
    }
    
    # depth
    if( 'min_tumour_depth' %in% names(filters) ) { 
        tumour.dp <- mean.field.value(filtered.variants, field = 'TUMOUR.DP', caller = caller); 
        
        passed.filter <- (tumour.dp >= filters$min_tumour_depth) | is.na(tumour.dp);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied min_tumour_depth filter, and removed', sum(!passed.filter), 'variants\n');
        }
        
        
    }
    
    if( paired && 'min_normal_depth' %in% names(filters) ) {
        normal.dp <- mean.field.value(filtered.variants, field = 'NORMAL.DP', caller = caller); 
        
        passed.filter <- (normal.dp >= filters$min_normal_depth) | is.na(normal.dp);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied min_normal_depth filter, and removed', sum(!passed.filter), 'variants\n');
        }
        
    }
    
    # allele frequency
    if( 'min_tumour_allele_frequency' %in% names(filters) ) {

        tumour.af <- mean.field.value(filtered.variants, field = 'TUMOUR.AF', caller = caller); 
        
        passed.filter <- (tumour.af >= filters$min_tumour_allele_frequency) | is.na(tumour.af);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied min_tumour_allele_frequency filter, and removed', sum(!passed.filter), 'variants\n');
        }

    } 
    
    if( paired && 'max_normal_allele_frequency' %in% names(filters) && caller != 'lofreq' ) {
        normal.af <- mean.field.value(filtered.variants, field = 'NORMAL.AF', caller = caller); 
        
        passed.filter <- (normal.af <= filters$max_normal_allele_frequency) | is.na(normal.af);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied max_normal_allele_frequency filter, and removed', sum(!passed.filter), 'variants\n');
        }
    }
    
    
    # higher threshold for indels
    if( 'indel_min_tumour_allele_frequency' %in% names(filters) ) {
        # not sure if data frame will always have this... recalculate
        mutation.type <- classify.variant(ref = filtered.variants$REF, alt = filtered.variants$ALT);
        tumour.af <- mean.field.value(filtered.variants, field = 'TUMOUR.AF', caller = caller);
        
        passed.filter <- ('indel' != mutation.type) | tumour.af >= filters$indel_min_tumour_allele_frequency | is.na(tumour.af);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied indel_min_tumour_allele_frequency filter, and removed', sum(!passed.filter), 'variants\n');
        }
    }
    
    
    
    # quality
    if( 'min_quality' %in% names(filters) ) {
        qual <- mean.field.value(filtered.variants, field = 'QUAL', caller = caller);
        
        passed.filter <- (qual >= filters$min_quality) | is.na(qual);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        
        if(verbose) { 
            cat('Applied min_quality filter, and removed', sum(!passed.filter), 'variants\n');
        }		
    }
    
    # FFPE artefact filter 
    if( 'ct_min_tumour_allele_frequency' %in% names(filters) ) { 
        
        tumour.af <- mean.field.value(filtered.variants, field = 'TUMOUR.AF', caller = caller);
        base.substitutions <- get.base.substitution(ref = filtered.variants$REF, alt = filtered.variants$ALT);
        
        # replace NA with blank to avoid selection headaches		
        base.substitutions[ is.na(base.substitutions) ] <- '';
        
        passed.filter <- (base.substitutions != 'C>T') | (tumour.af >= filters$ct_min_tumour_allele_frequency) | is.na(tumour.af);
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied ct_min_tumour_allele_frequency filter, and removed', sum(!passed.filter), 'variants\n');
        }
    }
    
    # appearance in 1000 genomes
    if( 'remove_1000_genomes' %in% names(filters) && filters$remove_1000_genomes ) { 
        
        # try not to hardcode the 1000 genomes version, in case it updates
        if( sum(grepl('1000g', names(variants) )) > 1 ) {
            stop('More than one column matching 1000g found.');
        }
        
        passed.filter <- '.' == filtered.variants[ , grepl('1000g', names(variants))];
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied remove_1000_genomes filter, and removed', sum(!passed.filter), 'variants\n');
        }		
    }

    if( 'remove_exac' %in% names(filters) && filters$remove_exac ) { 
        
        # try not to hardcode the 1000 genomes version, in case it updates
        if( sum(grepl('^ExAC', names(variants) )) > 1 ) {
            stop('More than one column matching ExAC found.');
        }
        
        # keep anything at variant allele frequency below 0.01
        passed.filter <- '.' == filtered.variants[ , grepl('^ExAC', names(variants))] | as.numeric(filtered.variants[ , grepl('^ExAC', names(variants))]) < 0.01;
        filtered.variants <- filtered.variants[passed.filter, ];
        
        if(verbose) { 
            cat('Applied remove_exac filter, and removed', sum(!passed.filter), 'variants\n');
        }       
    }

    
    if( 'remove_germline_status' %in% names(filters) && filters$remove_germline_status ) {
        
        status.columns <- grepl('STATUS$', names(filtered.variants));
        
        if( 1 == sum(status.columns) ) {
            passed.filter <- is.na(filtered.variants[, status.columns]) | filtered.variants[, status.columns] != 'Germline';			
            filtered.variants <- filtered.variants[passed.filter, ];
        } 
        
        if(verbose) { 
            cat('Applied remove_germline_status filter, and removed', sum(!passed.filter), 'variants\n');
        }				
        
    }
    
    
    # change back to original headers
    names(filtered.variants) <- old.names;
    
    return(filtered.variants);
}
