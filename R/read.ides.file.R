#' Read iDES output
#'
#' @description
#' Read output from iDES_step1.pl and return data frame
#' 
#' @param  filename path to file
#' 
#' @return ides.data  data frame read from iDES output
read.ides.file <- function(filename) {
    
    ides.data <- utils::read.table(
        filename,
        sep = '\t', 
        as.is = TRUE, 
        header = TRUE
        );
    
    # TO DO: why inconsistency in Rplus vs Apos ? 
    colnames(ides.data) <-c(
        'Chr', 'Pos', 'Depth', 'Ref', 'Rplus', 'Rneg', 'Apos', 
        'Aneg', 'Cpos', 'Cneg', 'Tpos', 'Tneg', 'Gpos', 'Gneg'
        );
    
    return(ides.data);	
}


#' Merge potential iDES calls with variant annotation.
#' 
#' @details
#' The VarDict variant calling includes a GATK call merging the call vcf file (allele frequency information etc.) with
#' the ANNOVAR annotation, and saving the result as a table. This function is an attempt to emulate that step 
#' for the iDES calls. 
#'
#' @param ides.filename 
#'  Path to formatted iDES output (typically from convert.ides.output file)
#' @param annovar.suffix.pattern 
#'  Suffix to match ANNOAR file
#' @inheritParams convert.ides.output
#'  
#' @return annotated.calls Data frame of annotations and iDES output.
merge.ides.annotation <- function(
    ides.filename, 
    output = TRUE, 
    output.suffix = '.ann.txt', 
    annovar.suffix.pattern = '.annovar.hg(\\d{2})_multianno.txt' 
    ) { 
    
    # open ides file
    ides.calls <- utils::read.delim(
        ides.filename, 
        sep = '\t',
        header = FALSE
    );
    
    # TO DO:
    #	- make sure VarDict depth (DP) has same definition  
    names(ides.calls) <- c(
        'Chr', 'Start', 'End', 'Ref', 'Alt', 'DP', 
        'RefCalls', 'AltCalls', 'A', 'C', 'G', 'T', 'AF'
        );
    
    # check if annovar file exists
    # regex pattern to match filename of ANNOVAR annotation file
    annovar.pattern <- paste0(basename(ides.filename), annovar.suffix.pattern)
    
    # match regex in same directory as ides file
    annovar.file.matches <- list.files(
        pattern = annovar.pattern, 
        path = dirname(ides.filename)
        );
    
    if( 0 == length(annovar.file.matches) ) {
        stop('No ANNOVAR annotation file found.');
    }
    
    annovar.annotation <- utils::read.delim(
        file.path(dirname(ides.filename), annovar.file.matches[1]),
        sep = '\t', 
        header = TRUE
        );
    
    # merge ANNOVAR annotation and iDES output
    merged.data <- merge(ides.calls, annovar.annotation);
    
    # coerce to data frame (not sure why merge returns a list), and make sure columns appear with chr, start, end, ref, alt first
    annotated.calls <- data.frame(merged.data)[, union(names(ides.calls), names(annovar.annotation))];
    
    # sort by chromosome and position (start and end position are equal since only dealing with SNVs)
    # merging step treats chromosome as a character, sorting does not work as expected
    annotated.calls <- annotated.calls[order(annotated.calls$Chr, annotated.calls$Start), ];
    
    
    # rename to match GATK VariantsToTable output
    # (both isis and VarDict use that tool)
    annotated.calls$End <- NULL;
    names(annotated.calls)[1:4] <- c('CHROM', 'POS' ,'REF', 'ALT');
    
    
    # write to file if requested
    if( TRUE == output ) {
        utils::write.table(
            annotated.calls,
            paste0(ides.filename, output.suffix),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
            );	
    }
    
    return(annotated.calls);
}