#' Merge variants
#'
#' @description 
#'	Merge variants from multiple callers and return a data frame of merged calls. By default filtering
#' 	is also applied, although this behaviour can be turned off by setting apply.filters to FALSE. 
#'
#' @param variant.specification 
#' 	Data frame containing details of file paths, sample IDs, and caller.
#' @param apply.filters Logical indicating whether to apply filters. Defaults to TRUE.
#' @param separate.consensus.filters 
#'  Logical indicating whether to apply different thresholds to variants called by more than one caller 
#' 	(specified under consensus in config file). Defaults to FALSE.
#' @param remove.structural.variants 
#'  Logical indicating whether structural variants (including CNVs) should be removed. Defaults to TRUE. 
#' @param verbose Logical indicating whether to print information to screen
#'
#' @return Data frame 
merge.variants <- function(
    variant.specification, 
    apply.filters = TRUE,
    remove.structural.variants = TRUE, 
    separate.consensus.filters = FALSE, 
    verbose = FALSE
    ) { 
    
    # TO DO: 
    #	- file exist tests..?	
    #	- verbose output within function itself
    
    ### INPUT TESTS ###########################################################
    
    # dummy test until you figure out what to do
    if( !(all( c('sample.id', 'variant.file', 'caller') %in% names(variant.specification) ))) {
        stop('variant.specification is missing required columns');
    }
    
    ### MAIN ##################################################################
    
    filter.immediately <- apply.filters;
    if( separate.consensus.filters ) { 
        filter.immediately <- FALSE;
    }
    
    # loop over requested variant callers
    # make data frame with all variants from that caller, and merge with variants from other callers
    
    merged.variants <- NULL;
    variant.callers <- unique(variant.specification$caller);

    # keep track of variant callers that actually called mutations
    for( caller in variant.callers ) { 
        caller.variant.files <- variant.specification[caller == variant.specification$caller, ];
        
        # store caller-specific variants
        caller.variants <- list();
        
        for( i in 1:nrow(caller.variant.files) ) { 
            variant.file <- caller.variant.files$variant.file[i];
            caller <- caller.variant.files$caller[i];
            sample.id <- caller.variant.files$sample.id[i];
            
            variant.calls <- read.variant.calls(variant.file, variant.caller = caller);
            
            if( !is.null(variant.calls) ) { 
                caller.variants[[ as.character(sample.id) ]] <- data.frame(
                    sample.id = sample.id, 
                    variant.calls,
                    stringsAsFactors = FALSE
                );
            }
            
        }

        # if no variants found, skip ahead
        if( 0 == length(caller.variants) ) {
            print('Skipping ahead');
            next;
        }

        # make sure all variants have all columns
        # want to preserve order, so this may seem somewhat roundabout
        column.names <- lapply(caller.variants, names);
        unique.column.names <- unique( unlist(column.names) );
        
        full.column.names <- column.names[[ which.max( sapply(column.names, length) ) ]];

        # if any are missing from longest column names, add to it
        if( !all(unique.column.names %in% full.column.names) ) {
            full.column.names <- c(
                full.column.names, 
                unique.column.names[ !(unique.column.names %in% full.column.names) ]
                );
        }

        temp.caller.variants <- list();

        for(i in 1:length(caller.variants) ) {

            temp.data <- caller.variants[[ i ]];
            missing.columns <- full.column.names[ !(full.column.names %in% names(temp.data) ) ];

            for( column in missing.columns ) temp.data[, column] <- '.';

            temp.caller.variants[[ i ]] <- temp.data[, full.column.names];
        }

        caller.variants <- do.call(rbind, temp.caller.variants);
        
        # note which caller called the variant
        caller.variants[, paste0('CALLED.', toupper(caller))] <- caller;
        
        
        # apply filters if requested
        if( filter.immediately ) { 
            caller.variants <- filter.variants(
                caller.variants,
                caller = caller,
                verbose = verbose
            );
        }
        
        # merge with other variant callers
        if( is.null(merged.variants) ) { 
            merged.variants <- caller.variants;
        } else { 
            merged.variants <- merge(
                merged.variants, 
                caller.variants, 
                all = TRUE
            );
        }
        
    }

    print( utils::str(merged.variants) );
    print(paste0('CALLED.', toupper(variant.callers) ));

    ### POST-PROCESSING 
    if( length(variant.callers) > 1) {
        # throws an error if only one caller, handle separately
        merged.variants$caller <- apply(
            merged.variants[,  grepl('^CALLED\\.', names(merged.variants))], 
            1,
            FUN = function(x) paste(x[!is.na(x)], collapse = ':')
        );	
    } else { 
        merged.variants$caller <- variant.callers;
    }
    
    merged.variants[, paste0('CALLED.', toupper(variant.callers)) ] <- NULL;
    
    merged.variants$Type <- classify.variant(ref = merged.variants$REF, alt = merged.variants$ALT);
    
    
    ### FILTERING IF CONSENSUS RESCUING HAS BEEN REQUESTED
    
    # use filtered.variants from here on out – if merged.variants appears again it is a bug!!
    filtered.variants <- merged.variants;
    
    if( apply.filters && separate.consensus.filters ) { 
        # want to identify variants that have been called by more than one variant caller, and apply less stringent
        # filters to them
        
        for(caller in variant.callers) { 
            
            # identify variants that have only been called by this variant, and apply filters
            caller.only <- caller == filtered.variants$caller;
            
            
            filtered.caller.variants <- filter.variants(
                filtered.variants[caller.only, ], 
                caller = caller,
            );
            

            filtered.variants <- rbind(
                filtered.variants[!caller.only, ], 
                filtered.caller.variants
            );
            
        }
        
        # FILTER CONSENSUS CALLS
        
        # get rows corresponding to variants called by multiple callers
        multiple.callers <- grepl(':', filtered.variants$caller);
        
        
        filtered.consensus.variants <- filter.variants(
            filtered.variants[multiple.callers, ], 
            caller = 'consensus'
        );

        filtered.variants <- rbind(
            filtered.variants[!multiple.callers, ], 
            filtered.consensus.variants		
            );
        
    }
    
    ### REMOVE SVs IF REQUESTED
    
    # generally do not trust CNV calling from targeted panel
    # remove all entries with anything other than a specific ALT base, e.g. <CNV>
    if( remove.structural.variants ) { 
        filtered.variants <- filtered.variants[!grepl('<', filtered.variants$ALT), ];	
    }

    ## MEAN TUMOUR ALLELE FREQUENCY

    is.tumour.af.column <- grepl('TUMOUR.AF$', names(filtered.variants));
    filtered.variants$MEAN.TUMOUR.AF <- apply(
        filtered.variants[, is.tumour.af.column],
        1,
        FUN = function(AFs) {
            # remove blank values – caller doesn't have an estiwarnmate
            AFs <- AFs[ !is.na(AFs) & '' != AFs ];

            return( mean( as.numeric(AFs)) );
            }
        );
    
    return(filtered.variants);
    
}