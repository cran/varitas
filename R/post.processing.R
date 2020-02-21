#' Post-processing of variants to generate outputs
#'
#' @param variant.specification 
#'  Data frame specifying variants to be processed, or path to data frame (useful if calling from Perl)
#' @param project.directory 
#' 	Directory where output should be stored. Output files will be saved to a datestamped subdirectory
#' @param config.file
#'  Path to config file specifying post-processing options. If not provided, the current options are used (i.e. from \code{get.varitas.options()})
#' @param variant.callers
#'  Optional vector of variant callers for which filters should be included in Excel file
#' @param sleep 
#'  Logical indicating whether script should sleep for 60 seconds before starting.
#' @param verbose
#' 	Logical indicating whether to print verbose output
#' @inheritParams merge.variants
#'
#' @return None
#'
#'
post.processing <- function(
    variant.specification, 
    project.directory,
    config.file = NULL,
    variant.callers = NULL,
    remove.structural.variants = TRUE, 
    separate.consensus.filters = FALSE, 
    sleep = FALSE,
    verbose = FALSE
    ) { 


    # some problems with hard drive being too slow to save all files
    # try to get around this by sleeping and hoping this is enough to let things settle
    if( sleep ) {
        Sys.sleep(5*60);
    }


    # TO DO:
    #	- specification data frame for coverage QC? (if we could bugfix bedtools that would be ideal...)
    
    ### INPUT TESTS ###########################################################
    
    # project.directory
    assertthat::assert_that( is.character(project.directory) );
    assertthat::assert_that( length(project.directory) == 1 );
    assertthat::assert_that( 
        dir.exists(project.directory), 
        msg = paste(
            'Directory', project.directory, 
            'does not exist or is not a directory'
            )
    );
    
    # variant.specification
    assertthat::assert_that( 
        is.character(variant.specification) || is.data.frame(variant.specification), 
        msg = 'variant.specification should be either a data frame of a file path'
    );
    
    assertthat::assert_that( 
        is.data.frame(variant.specification) || file.exists(variant.specification), 
        msg = paste(
            'variant.specification file', variant.specification, 
            'does not exist'
            )
    );
    
    ### MAIN ##################################################################
    
    if( !is.null(config.file) ) overwrite.varitas.options(config.file);
    
    if( is.character(variant.specification) ) {
        
        variant.specification <- utils::read.table(
            variant.specification, 
            sep = '\t', 
            header = TRUE, 
            stringsAsFactors = FALSE
        );
    }

    # keep track of all sample IDs
    sample.ids <- unique( variant.specification$sample.id );
    
    # Get callers
    all.callers <- variant.specification$caller;
    
    unique.callers <- unique(all.callers);
    
    if( is.null(variant.callers) ) {
      variant.callers <- unique.callers
    } else {
      variant.callers <- union(variant.callers, unique.callers);
    }
    
    # Fix allele frequencies if necessary
    if ( "varscan" %in% variant.callers ) {
      fix.varscan.af(variant.specification)
    }
    
    if ( "lofreq" %in% variant.callers ) {
      fix.lofreq.af(variant.specification)
    }
    
    # create output directory and directory for plots
    # this directory contains the files that are meant to be sent to collaborators
    #   (variant calls, read statistics, PDF report)
    output.directory <- file.path(
        project.directory, 
        paste0(Sys.Date(), '-variant-data')
        );
    
    # directory for raw plots
    # not meant to be sent to collaborators, but contains plots from PDF report in
    # PNG format
    plotting.directory <- file.path(
        project.directory,
        paste0(Sys.Date(), '-plots')
        );

    if( !dir.exists(output.directory) ) dir.create(output.directory);
    if( !dir.exists(plotting.directory) ) dir.create(plotting.directory);

    ### QUALITY CONTROL DATA 
    if( length( system.ls('*/*.stats', project.directory) ) > 0  ) {

        coverage.sample <- get.coverage.by.sample.statistics(project.directory);
        coverage.amplicon <- get.coverage.by.amplicon(project.directory);

        # plot amplicon coverage per sample
        sample.scatterplot.directory <- file.path(plotting.directory, 'sample-coverage');
        if( !dir.exists(sample.scatterplot.directory) ) dir.create(sample.scatterplot.directory);
        
        plot.amplicon.coverage.per.sample(
            coverage.amplicon,
            sample.scatterplot.directory
            );

        # Ontarget reads
        ontarget.percent.file.name <- file.path(plotting.directory, 'ontarget_percent.png');
        plot.ontarget.percent(
            coverage.sample, 
            file.name = ontarget.percent.file.name
            );

        # Paired reads – not the most informative plot, but useful to detect if pairing worked
        paired.percent.file.name <- file.path(plotting.directory, 'paired_percent.png');
        plot.paired.percent(
            coverage.sample, 
            file.name = paired.percent.file.name
            );

        # Mean coverage per sample
        mean.coverage.file.name <- file.path(plotting.directory, 'mean_coverage.png');
        plot.coverage.by.sample(
            coverage.sample[order(-coverage.sample$mean.coverage),], 
            file.name = mean.coverage.file.name,
            statistic = 'mean'
            );

        # median coverage per sample
        median.coverage.file.name <- file.path(plotting.directory, 'median_coverage.png');
        plot.coverage.by.sample(
            coverage.sample[order(-coverage.sample$median.coverage),], 
            file.name = median.coverage.file.name,
            statistic = 'median'
            );
        
        # Coverage statistics Excel file
        # Should be sent to collaborators, i.e. save in output directory
        file.name <- file.path(output.directory, 'Coverage_statistics.xlsx');
        save.coverage.excel(
            project.directory = project.directory, 
            file.name = file.name,
            overwrite = TRUE
        );
        
    }
    
    
    ### GET VARIANTS
    
    # should these be parameterized? 
    filtered.variants <- merge.variants(
        variant.specification, 
        apply.filters = TRUE,
        remove.structural.variants = TRUE, 
        separate.consensus.filters = FALSE, 
        verbose = verbose
        );

    # keep track of all sample IDs so they appear in plots
    filtered.variants$sample.id <- factor(
        filtered.variants$sample.id,
        levels = sample.ids
        );
        
    
    # save to txt -> probably a good idea to keep raw data 
    utils::write.table(
        filtered.variants, 
        file.path(output.directory, 'filtered_variants.txt'), 
        sep = '\t',
        row.names = FALSE
        );
    
    
    # save to Excel
    
    filters <- get.varitas.options('filters')[ variant.callers ];
    
    save.variants.excel(
        variants = filtered.variants, 
        file.name = file.path(output.directory, 'Filtered_variants.xlsx'), 
        filters = filters
        );
    
    ### PLOTS
    
    # Caller overlap venn diagram
    caller.overlap.venn.diagram(
        filtered.variants,
        file.name = file.path(plotting.directory, 'caller_overlap.png')
        );
    
    variants.sample.barplot(
        filtered.variants, 
        file.name = file.path(plotting.directory, 'variants_per_sample.png')
        );
    
    variants.caller.barplot(
        filtered.variants, 
        file.name = file.path(plotting.directory, 'variants_caller_type.png'), 
        group.by = 'type'
        );
    
    variants.caller.barplot(
        filtered.variants, 
        file.name = file.path(plotting.directory, 'variants_caller_substitution.png'), 
        group.by = 'substitution'
        );

    trinucleotide.barplot(
        filtered.variants, 
        file.name = file.path(plotting.directory, 'trinucleotide_substitutions.png')
        );

    variant.recurrence.barplot(
        filtered.variants, 
        file.name = file.path(plotting.directory, 'variant_recurrence.png') 
        );

    ### REPORT

    report.template <- system.file('report_template.Rmd', package = 'varitas');

    rmarkdown::render(
        report.template, 
        output_format = 'pdf_document',
        output_file = file.path(output.directory, 'pipeline_report.pdf')
        );

}
