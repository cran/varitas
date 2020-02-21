#' read.all.calls
#'
#' @description
#'  Read all calls made with a certain caller 
#'
#' @param sample.ids 
#'  Vector giving sample IDs to process
#' @param caller 
#'  String indicating which caller was used
#' @param project.directory 
#'  Path to project directory
#' @param patient.ids 
#'  Optional vector giving patient ID (or other group) corresponding to each sample 
#' @param apply.filters 
#'  Logical indicating whether filters specified in VariTAS options should be applied. Defaults to TRUE. !
#' @param variant.file.pattern 
#'  Pattern indicating where the variant file can be found. Sample ID should be indicated by SAMPLE_ID  
#'
#' @return combined.variant.calls Data frame with variant calls from all patients
#'
#'
read.all.calls <- function(
    sample.ids, 
    caller = c('vardict', 'mutect', 'pgm'), 
    project.directory, 
    patient.ids = NULL, 
    apply.filters = TRUE,
    variant.file.pattern = NULL
    ) { 
    
    # TO DO: 
    #	- move fix.names option to a two-tier system (no need to add variant caller to start with)
    #	- add support for data frame specification – much more robust
    
    ### INPUT TESTS ###########################################################
    
    if( !is.character(project.directory) ) { 
        stop('project.directory must be a character string');
    }
    
    if( !dir.exists(project.directory) ) { 
        error.message <- paste('Directory', project.directory, 'not found');
        stop(error.message);
    }
    
    if( !is.character(sample.ids) ) { 
        stop('sample.ids must be character vector');
    }
    
    if( !is.null(patient.ids) && !is.character(patient.ids) ) { 
        stop('patient.ids must be a character vector');
    }
    
    ### MAIN ##################################################################
    
    caller <- match.arg(caller);
    
    # match caller to a pattern indicating where variants can be found
    # SAMPLE_ID will be replaced with sample.id variable
    file.path.patterns <- c(
        'vardict' = 'SAMPLE_ID/vardict/SAMPLE_ID.passed.ontarget.vcf.annovar.hg19_multianno.vcf.txt', 
        'mutect' = 'SAMPLE_ID/mutect/SAMPLE_ID.passed.ontarget.vcf.annovar.hg19_multianno.vcf.txt', 
        'pgm' = 'SAMPLE_ID/pgm/SAMPLE_ID.snvs_and_cnvs.vcf.annovar.hg19_multianno.vcf.txt'
    );
    
    uppercase.caller <- toupper(caller);
    
    combined.variant.calls <- list();
    
    # loop over sample IDs and read in variants
    for(i in 1:length(sample.ids) ) { 
        
        sample.id <- sample.ids[i];
        
        
        if( is.null(variant.file.pattern) ) { 
            variant.file.pattern <- file.path.patterns[ caller ];
        }
        
        variant.file.path <- file.path(
            project.directory, 
            gsub('SAMPLE_ID', sample.id, variant.file.pattern)
        );
        
        # add warning? 
        if( !file.exists(variant.file.path) ) {
            warning.message <- paste('File', variant.file.path, 'not found.');
            warning(warning.message);
            next;
        }
        
        variant.calls <- utils::read.table(
            variant.file.path,
            header = TRUE,
            sep = '\t', 
            stringsAsFactors = FALSE
        );
        
        # if no variants in file, skip ahead			
        if( 0 == nrow(variant.calls) ) {
            warning.message <- paste('File', variant.file.path, 'exists, but does not contain any variants');
            warning(warning.message);
            next;
        }	
        
        names(variant.calls) <- fix.names(names(variant.calls), caller);
        
        # add number of reads supporting variant
        variant.calls[, paste0(uppercase.caller, '.VARIANT.READS')] <- variant.calls[, paste0(uppercase.caller, '.TUMOUR.AF')]*variant.calls[, paste0(uppercase.caller, '.TUMOUR.DEPTH')];
        
        # add to list
        #	- need to do this separately based on whether patient.ids is defined (want PatientID to be first column
        if( !is.null(patient.ids) ) { 
            combined.variant.calls[[ sample.id ]] <- data.frame(
                PatientID = rep(patient.ids[i], nrow(variant.calls)), 
                SampleID = rep(sample.id, nrow(variant.calls)), 
                variant.calls
            );
            
        } else { 
            combined.variant.calls[[ sample.id ]] <- data.frame(
                SampleID = rep(sample.id, nrow(variant.calls)), 
                variant.calls
            );
        }
        
    } # end of per sample for loop 
    
    # consolidate to data frame 	
    
    combined.variant.calls <- do.call(rbind, combined.variant.calls);
    
    combined.variant.calls[, paste0('CALLED.', uppercase.caller)] <- tolower(caller);
    
    if(apply.filters) { 
        combined.variant.calls <- filter.variants(
            variants = combined.variant.calls, 
            caller = caller 
        );
    }
    
    return(combined.variant.calls);
    
}