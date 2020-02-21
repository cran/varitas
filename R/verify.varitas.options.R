#' Check against common errors in the VariTAS options. 
#'
#' @description
#'	Check against common errors in the VariTAS options before launching into pipeline
#'
#' @param stages.to.run 
#'	Vector indicating which stages should be run. Defaults to all possible stages. 
#'	If only running a subset of stages, only checks corresponding to the desired stages are run
#' @param variant.callers
#'	Vector indicating which variant callers to run. Only used if calling is in \code{stages.to.run}.
#' @param varitas.options
#'	Optional file path or list of VariTAS options.  	
#'
#' @return None
#'
#'
verify.varitas.options <- function(
    stages.to.run = c('alignment', 'qc', 'calling', 'annotation', 'merging'), 
    variant.callers = c('mutect', 'vardict', 'ides', 'varscan', 'lofreq', 'muse'), 
    varitas.options = NULL
) {
    
    ## TO DO:
    #	 - ordering of reference genome and VCF files for MuTect
    #	 - chr-compatability of target panel and reference genome
    #    - check if QC needs any specific settings
    
    stages.to.run <- tolower(stages.to.run);
    
    ### INPUT TESTS ###########################################################
    
    # NOTE: 
    # 	- Merging currently does not require any specific config settings, so no checks are run. For ease-of-use, the stage is still included as an option.
    supported.stages <- c('alignment', 'qc', 'calling', 'annotation', 'merging');
    if( !all(stages.to.run %in% supported.stages) ) {
        unrecognized.stages <- stages.to.run[ !(stages.to.run %in% supported.stages) ];
        
        error.message <- paste('The following stages are not supported:', paste(unrecognized.stages, collapse = ', '));
        stop(error.message);
    }
    
    assertthat::assert_that(
        is.null(varitas.options) || is.list(varitas.options) || is.character(varitas.options), 
        msg = 'varitas.options must be a list of options or a string giving the path to the config YAML file.'
    );
    
    assertthat::assert_that(
        !is.character(varitas.options) || file.exists(varitas.options), 
        msg = paste('varitas.options file', varitas.options, 'does not exist' )
    );
    
    ### MAIN ##################################################################
    
    if( is.null(varitas.options) ) {
      varitas.options <- get.varitas.options();
    } else if( is.character(varitas.options)) {
        # NOTE: 
        # this is a bit messy, but the idea is that it should be possible to verify 
        # a possible set of options BEFORE setting them. 
        # 	=> allow both list, character, and no input
      varitas.options <- yaml::yaml.load_file(varitas.options);
    }
    ### GENERIC SETTINGS
    
    # is this needed for merging actually?
    assertthat::assert_that('reference_build' %in% names(varitas.options),msg = 'config must include reference_build');
    
    reference.build <- varitas.options$reference_build;
    assertthat::assert_that(
        reference.build %in% c('grch37', 'grch38'), 
        msg = 'reference_build must be either grch37 or grch38'
    );	
    
    ### REFERENCE GENOME
    # all stages except merging use reference genome 
    # 	=> check existence here
    if( !identical(stages.to.run, 'merging')) {
        
        # reference genome exists and has necessary derivative files
        reference.genome <- varitas.options$reference_genome[[ reference.build ]];
        assertthat::assert_that(
            file.exists(reference.genome),
            msg = paste('Reference genome file', reference.genome, 'not found')
        );
        
        # reference genome has fa or fasta extension 
        # (this is probably stricter than necessary)
        reference.genome.extension <- tools::file_ext(reference.genome);
        assertthat::assert_that(
            tolower(reference.genome.extension) %in% c('fa', 'fasta'), 
            msg = paste('Reference genome file', reference.genome, 'does not have extension .fa or .fasta')
        );
        
        # needed for future steps
        reference.genome.chromosomes <- get.fasta.chromosomes(reference.genome);
        # there's an upper limit on how long the error message can be
        # 	=> 
        if( length(reference.genome.chromosomes) > 25) {
            reference.chromosome.string <- paste(c(reference.genome.chromosomes[1:25], '...'), collapse = ' ');
        } else {
            reference.chromosome.string <- paste(reference.genome.chromosomes, collapse = ' ');
        }
        
    }
    
    ### TARGET PANEL
    # needed for both alignment and variant calling
    if( 'alignment' %in% stages.to.run || 'calling' %in% stages.to.run ) {
        
        assertthat::assert_that(
            'target_panel' %in% names(varitas.options), 
            msg = 'target_panel must be provided for alignment and variant calling stages'
        );
        
        target.panel <- varitas.options$target_panel;
        assertthat::assert_that(
            file.exists( target.panel ), 
            msg = paste('target_panel file', target.panel, 'does not exist')
        );
        
        panel.chromosomes <- get.bed.chromosomes(target.panel);	
        
        assertthat::assert_that(
            all(panel.chromosomes %in% reference.genome.chromosomes), 
            msg = paste(
                'Mismatch between reference genome and target panel.\n', 
                'Reference genome chromosomes:', reference.chromosome.string, '\n', 
                'Target panel chromosomes:', paste(panel.chromosomes, collapse = ' ')
            )
        );
        
        bed.data <- utils::read.table(target.panel)
        gene.column <- 5;
        if( !(any( grepl('GENE_ID', bed.data[, gene.column] ))) ) {
          gene.column <- 6;
        }
        if( !(any( grepl('GENE_ID', bed.data[, gene.column] ))) ) {
          stop('Target panel must contain gene/feature IDs in the 5th or 6th column\n  In the format of GENE_ID=...;etc.')
        }
        
    }
    
    
    ### ALIGNMENT-SPECIFIC OPTIONS
    
    if( 'alignment' %in% stages.to.run ) {
        
        # check that bwa index has been run on the reference genome
        verify.bwa.index(reference.genome, error = TRUE);
        verify.sequence.dictionary(reference.genome, error = TRUE);
        verify.fasta.index(reference.genome, error = TRUE);
        
        # Picard jar file
        assertthat::assert_that(
            'picard_jar' %in% names(varitas.options), 
            msg = 'config must contain picard_jar when running alignment'
        );
        
        picard.jar <- varitas.options$picard_jar;
        assertthat::assert_that(
            file.exists(picard.jar), 
            msg = paste('picard_jar file', picard.jar, 'does not exist')
        );
        
        # GATK jar file
        assertthat::assert_that(
            'gatk_jar' %in% names(varitas.options), 
            msg = 'config must contain gatk_jar if running alignment'
        );
        
        gatk.jar <- varitas.options$gatk_jar;
        assertthat::assert_that(
            file.exists(gatk.jar),
            msg = paste('gatk_jar file', gatk.jar, 'does not exist')
        );
        
        
    } 
    
    ### VARIANT CALLING 
    
    if( 'calling' %in% stages.to.run ) {
        
        verify.sequence.dictionary(reference.genome, error = TRUE);
        verify.fasta.index(reference.genome, error = TRUE);
        
        if( 'mutect' %in% variant.callers ) {
            # TO DO: 
            # 	- parameterize whether to use dbSNP and cosmic files
            
            # GATK jar file
            assertthat::assert_that(
                'gatk_jar' %in% names(varitas.options), 
                msg = 'config must contain gatk_jar if running Mutect'
            );
            
            gatk.jar <- varitas.options$gatk_jar;
            assertthat::assert_that(
                file.exists(gatk.jar),
                msg = paste('gatk_jar file', gatk.jar, 'does not exist')
            );
            
            
            
        }
        
        if( 'vardict' %in% variant.callers) {
            
            assertthat::assert_that(
                'vardict_path' %in% names(varitas.options),
                msg = 'config must include vardict_path if running VarDict'
            );
            
            vardict.path <- varitas.options$vardict_path;
            assertthat::assert_that(
                dir.exists( vardict.path ), 
                msg = paste('vardict_path directory', vardict.path, 'does not exist or is not a directory')
            );
        }
      
        if( 'varscan' %in% variant.callers ) {
          
          assertthat::assert_that(
            'varscan_path' %in% names(varitas.options),
            msg = 'config must include varscan_path if running varscan'
          );
          
          varscan.path <- varitas.options$varscan_path;
          assertthat::assert_that(
            file.exists( varscan.path ), 
            msg = paste('varscan_path directory', varscan.path, 'does not exist')
          );
        }
      
        if( 'lofreq' %in% variant.callers ) {
          
          assertthat::assert_that(
            'lofreq_path' %in% names(varitas.options),
            msg = 'config must include lofreq_path if running lofreq'
          );
          
          lofreq.path <- varitas.options$lofreq_path;
          assertthat::assert_that(
            file.exists( lofreq.path ), 
            msg = paste('lofreq_path directory', lofreq.path, 'does not exist')
          );
        }
      
        if( 'muse' %in% variant.callers ) {
          
          assertthat::assert_that(
            'muse_path' %in% names(varitas.options),
            msg = 'config must include muse_path if running muse'
          );
          
          muse.path <- varitas.options$muse_path;
          assertthat::assert_that(
            file.exists( muse.path ), 
            msg = paste('muse_path directory', muse.path, 'does not exist')
          );
        }
        
    } 
    
    
    ### ANNOTATION 
    if( 'annotation' %in% stages.to.run ) {
        
        # ANNOVAR path
        assertthat::assert_that(
            'annovar_path' %in% names(varitas.options), 
            msg = 'config must include annovar_path when running annotation'
        );
        
        annovar.path <- varitas.options$annovar_path;
        assertthat::assert_that( 
            dir.exists(annovar.path), 
            msg = paste('annovar_path directory', annovar.path, 'does not exist or is not a directory')
        );
        
        
        # ANNOVAR database
        annovar.database <- varitas.options$annovar_database[[ reference.build ]];
        assertthat::assert_that(
            dir.exists(annovar.database), 
            msg = paste('annovar_database directory', annovar.database, 'does not exist or is not a directory')
        );
        
        # TO DO: annotation-specific database
        required.databases <- c('cytoBand', 'cosmic70', 'sites\\.2015_08', 'exac03nontcga', 
                                'clinvar_20170130', 'nci60', 'icgc21', 'dbnsfp30a', 'dbnsfp31a_interpro');
        
        for( db in required.databases) {
            matched.files <- list.files(pattern = db, path = annovar.database);
            
            # should this be > 0 or > 1? Will it work without .idx files
            assertthat::assert_that(
                length(matched.files) > 0,
                msg = paste('Missing required database', db, 'in', annovar.database)
            );
        }
        
        
    } 
    
    
}
