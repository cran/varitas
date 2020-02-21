#' plot.amplicon.coverage.per.sample
#'
#' @description
#'	Create one scatterplot per sample, showing coverage per amplicon, and an additional plot giving the median
#'
#' @param coverage.statistics 
#'	Data frame containing coverage per amplicon per sample, typically from \code{get.coverage.by.amplicon}.
#' @param output.directory 
#'	Directory where per sample plots should be saved
#'
#' @return None
plot.amplicon.coverage.per.sample <- function(
	coverage.statistics, 
	output.directory
	) {	

	# TO DO: 
	#	- figure out a smarter way of extracting the gene
	#	- look into ordering by chromosome rather than lexicographic order

	# get first column containing a sample – first numeric one after chromosome coordinates
	first.sample.column <- 4;
	while( !is.numeric(coverage.statistics[, first.sample.column]) && first.sample.column <= ncol(coverage.statistics) ) {
		if( ncol(coverage.statistics) == first.sample.column ) {
			stop('Cannot find first sample column');
		}
		first.sample.column <- first.sample.column + 1;
	}	
	
	#config <- read.yaml(save.config())
	#target.panel <- read.table(config[['target_panel']], stringsAsFactors = FALSE)
	#first.sample.column <- ncol(target.panel) + 1
	
	#coverage.statistics <- alternate.gene.sort(coverage.statistics)
	
	genes <- get.gene(coverage.statistics);
	
	gene.start <- vapply(
	  unique(genes), 
	  function(x, genes) match(x, genes), 
	  genes = genes,
	  FUN.VALUE = 0
	);
	
	gene.end <- c(gene.start[-1], nrow(coverage.statistics) + 1);
	midpoints <- gene.start + (gene.end - gene.start)/2;

	chr.nums <- sapply(coverage.statistics$chr, function(x) substr(x, 4, nchar(x)))
	to.remain <- sapply(chr.nums, function(x) x != 'X' && x != 'Y')
	
	old.names <- names(coverage.statistics)
	coverage.statistics <- cbind(genes, chr.nums, coverage.statistics)
	names(coverage.statistics) <- c('gene', 'chr.no', old.names)
	first.sample.column <- first.sample.column + 2
	for(i in which(sapply(coverage.statistics, class) == "factor")) coverage.statistics[[i]] = as.character(coverage.statistics[[i]])
	
	colours <- c()
	shapes <- c()
	chr.palette = c(
	  '#8DD3C7', '#081D58', '#BEBADA', '#FB8072', '#CCEBC5', '#FDB462', '#999999', '#FCCDE5', '#FC8D59', '#35978F', 
	  '#F781BF', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#A65628', '#80B1D3', '#252525', '#A6761D',
	  '#B3DE69', '#F0027F', '#FFFFCC', '#FDDBC7', '#004529'
	)
	# print(length(sample.coverage))
	
	# Sort by chromosome and position (sex chromosomes at the end)
	sex.chr.rows <- coverage.statistics[!to.remain, ]
	coverage.statistics <- coverage.statistics[to.remain, ]
	coverage.statistics$chr.no <- as.integer(coverage.statistics$chr.no)
	coverage.order <- order(coverage.statistics$chr.no, coverage.statistics$start, coverage.statistics$end);
	coverage.statistics <- coverage.statistics[ coverage.order, ];
	sex.chr.order <- order(sex.chr.rows$chr.no, sex.chr.rows$start, sex.chr.rows$end);
	sex.chr.rows <- sex.chr.rows[ sex.chr.order, ];
	coverage.statistics <- rbind(coverage.statistics, sex.chr.rows)

	genes <- unique(coverage.statistics$gene)
	# chr.list <- c()
	# for (g in 1:length(genes)) {
	#   chr.list <- c(chr.list, unique(coverage.statistics[which(coverage.statistics$gene==genes[g]),2]))
	# }
	chr.list <- coverage.statistics$chr.no
	
	# Red and blue for odd/even chromosome numbers 
	# for (j in 1:length(chr.list)) {
	#   chr.ending <- chr.list[j]
	#   if (chr.ending == 'X' || chr.ending == 'Y') {
	#     colours <- c(colours, 'grey')
	#     shapes <- c(shapes, 23)
	#     next
	#   }
	#   chromosome <- as.integer(chr.ending)
	#   # print(substr(coverage.statistics[j,2], 4, nchar(coverage.statistics[j,2])))
	#   if (chromosome %% 2 == 0) {
	#     colour <- 'red'
	#     shape <- 21
	#   } else {
	#     colour <- 'blue'
	#     shape <- 22
	#   }
	#   shapes <- c(shapes, shape)
	#   colours <- c(colours, colour)
	# }
	
	# Alternating red and blue for each chromosome present
	red <- TRUE
	prev.chr <- ''
	for (j in 1:length(chr.list)) {
	  chr.ending <- chr.list[j]
	  if (chr.ending != prev.chr){
  	  if (red) {
  	    colours <- c(colours, 'red')
  	    shapes <- c(shapes, 21)
  	    red <- FALSE
  	  } else {
  	    colours <- c(colours, 'blue')
  	    shapes <- c(shapes, 22)
  	    red <- TRUE
  	  }
	  } else {
	    colours <- c(colours, colours[length(colours)])
	    shapes <- c(shapes, shapes[length(shapes)])
	  }
	  prev.chr <- chr.ending
	}
	
	# Colour for each chromosome
	# for (j in 1:length(chr.list)) {
	#   chr.ending <- chr.list[j]
	#   if (chr.ending == 'X') {
	#     colours <- c(colours, chr.palette[24])
	#   } else if (chr.ending == 'Y') {
	#     colours <- c(colours, chr.palette[25])
	#   } else {
	#     colours <- c(colours, chr.palette[as.integer(chr.ending)])
	#   }
	# }
	
	# loop over columns and plot all of them
	for(i in first.sample.column:ncol(coverage.statistics) ) {

		sample.id <- names(coverage.statistics)[i];
		sample.coverage <- coverage.statistics[, i];
		# sample.coverage <- tapply(coverage.statistics[, i], coverage.statistics$gene, sum)
		# sample.coverage <- sample.coverage[match(genes, names(sample.coverage))]
		x <- c()
		for (g in 1:length(unique(genes))) {
		  for (i in 1:length(which(coverage.statistics$gene == unique(genes)[g]))) {
		    x <- c(x, g)
		  }
		}
		grDevices::png(
			file.path(output.directory, paste0(sample.id, '.png')),
			height = 4,
			width = 7,
			units = 'in',
			res = 400
			);
		
		graphics::par(
			mar =  c(3.2, 4, 1.2, 0.2),
			cex.axis = 0.6,
			font.axis = 1,
			oma = c(0, 0, 0, 0),
			las = 2,
			tcl = -0.2
			);

		graphics::plot(
			x = jitter(x, amount = 0.15),
			y = sample.coverage,
			main = sample.id,
			cex = 0.8,
			pch = shapes,
			# pch = 21,
			bg = colours,
			col = 'black',
			xlab = '',
			ylab = 'Coverage',
			xaxt = 'n',
			xaxs = 'r'
			);

		#graphics::abline(v = gene.start[-1], col = 'grey', lty = 'dashed');
		graphics::axis(1, at = 1:length(unique(genes)), labels = unique(genes), font = 2);
		
		grDevices::dev.off();

	}

	# TO DO: consolidate this plot with the one above into a function
	# this will also solve vignette issue
	# plot median
	if (first.sample.column < ncol(coverage.statistics)) {
	  avg.coverage.stats <- stats::aggregate.data.frame(
	    coverage.statistics[, first.sample.column:ncol(coverage.statistics)], 
	    list(coverage.statistics$gene), 
	    sum
	    );
	  avg.coverage.stats <- avg.coverage.stats[match(genes, avg.coverage.stats$Group.1),]
  	median.coverage <- apply(
  		avg.coverage.stats[, 2:ncol(avg.coverage.stats)],
  		1,
  		stats::median
  		);
	} else { # Only one sample
	  median.coverage <- stats::median(coverage.statistics[, first.sample.column])
	}

	grDevices::png(
		file.path(output.directory, 'median.png'),
		height = 4,
		width = 7,
		units = 'in',
		res = 400
		);
	
	graphics::par(
	  mar =  c(3.2, 4, 1.2, 0.2),
	  cex.axis = 0.6,
	  font.axis = 1,
	  oma = c(0, 0, 0, 0),
	  las = 2,
	  tcl = -0.2
	);
	
	graphics::plot(
	  x = seq_along(median.coverage),
	  y = median.coverage,
	  main = 'Median Coverage',
	  cex = 0.8,
	  pch = 22,
	  bg = 'grey',
	  col = 'black',
	  xlab = '',
	  ylab = 'Coverage',
	  xaxt = 'n',
	  xaxs = 'r'
	);

	#graphics::abline(v = gene.start[-1], col = 'grey', lty = 'dashed');
	graphics::axis(1, at = 1:length(unique(genes)), labels = unique(genes), font = 2);
	
	grDevices::dev.off();

}
