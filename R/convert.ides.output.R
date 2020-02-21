
#' Convert output of iDES step 1 to variant call format
#' 
#' @param filename Path to file
#' @param output Logical indicating whether output should be saved to file. Defaults to true.
#' @param output.suffix Suffix to be appended to input filename if saving results to file
#' @param minreads Minimum numbers of reads
#' @param mindepth Minimum depth
#'
#' @return potential.calls Data frame of converted iDES calls
#'
#'
convert.ides.output <- function(
		filename, 
		output = TRUE, 
		output.suffix = '.calls.txt', 
		minreads = 5, 
		mindepth = 50 ) {

	### INPUT TESTS
	
	
	### MAIN
		
	results <- read.ides.file(filename);
    
    results <- results[results$Depth >= mindepth, ];
    ref.calls <- results$Rplus + results$Rneg; # why is this Rplus when the others are Apos, Gpos, etc. ?
       
    # determine total number of reads supporting each base
    A <- results$Apos + results$Aneg;
    C <- results$Cpos + results$Cneg;
    T <- results$Tpos + results$Tneg;
    G <- results$Gpos + results$Gneg;
    
    reads.per.base <- data.frame(A, C, T, G);
    
    alt.base <- colnames(reads.per.base)[apply(reads.per.base, 1, which.max)];
    alt.base.calls <- apply(reads.per.base, 1, max);
    
    # Note: this needs to match Annovar format, hence the duplicated position.
	# From Annovar documentation: 
	#	On each line, the first five space- or tab- delimited columns represent 
	#	chromosome, start position, end position, the reference nucleotides and the 
	#	observed nucleotides.
    all.data <- data.frame(
    	'Chr' = results$Chr,
    	'Start' = results$Pos,
   		'End' = results$Pos,
   		'Ref' = results$Ref,
    	'Alt' = alt.base, 
    	'Depth' = results$Depth,
    	ref.calls, 
    	alt.base.calls,
    	reads.per.base
    	);
    	        
    # calculate fraction of reads supporting alternate allele
    all.data$AF <- all.data$alt.base.calls/all.data$Depth; 
    
    # filter on minimum number of reads supporting alternate allele
    potential.calls <- all.data[all.data$alt.base.calls > minreads,];
    
    # write output to file if requested
    if(TRUE == output) {
    
        utils::write.table(
    		potential.calls,
    		file = paste0(filename, output.suffix),
    		sep = "\t",
    		row.names = FALSE,
    		quote = FALSE,
    		col.names = FALSE
    		);  		
    } 
    
    return(potential.calls);
    
}