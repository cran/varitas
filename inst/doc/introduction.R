## ----setup, include=FALSE------------------------------------------------
knitr::opts_knit$set(root.dir = '../')


## ----fastq---------------------------------------------------------------
library(varitas);

output.directory <- '';

fastq.specification <- data.frame(
    sample.id = c('A', 'B', 'C', 'D'), 
    patient.id = c('X', 'X', 'Y', 'Y'), 
    tissue = c('tumour', 'normal', 'tumour', 'normal'), 
    reads = c('A_1.fq', 'B_1.fq', 'C_1.fq', 'D_1.fq'),
    mates = c('A_2.fq', 'B_2.fq', 'C_2.fq', 'D_2.fq')
    );

print(fastq.specification);

## ----alignment, results="hide"-------------------------------------------
matched.bam.specification <- run.alignment(
    fastq.specification = fastq.specification, 
    output.directory = output.directory, 
    paired.end = TRUE,
    quiet = TRUE # only for testing, does not submit jobs to cluster
    );

## ------------------------------------------------------------------------
print(matched.bam.specification);

## ----variants1-----------------------------------------------------------
unmatched.bam.specification <- data.frame(
    sample.id = c('Z', 'Y'), 
    tumour.bam = c('Z.bam', 'Y.bam')
    );

print(unmatched.bam.specification);

## ----variants2, results = FALSE------------------------------------------
vcf.specification <- run.variant.calling(
    matched.bam.specification, 
    output.directory = output.directory, 
    variant.caller = c('vardict', 'mutect'),
    quiet = TRUE # only for testing, does not submit jobs to cluster
    );

## ------------------------------------------------------------------------
print(vcf.specification);

## ---- results = FALSE----------------------------------------------------
variant.specification <- run.annotation(
    vcf.specification, 
    output.directory = output.directory, 
    quiet = TRUE # testing only
    );

## ---- eval = FALSE-------------------------------------------------------
#  print(variant.specification);

## ---- results = FALSE----------------------------------------------------
run.post.processing(
	variant.specification = variant.specification, 
	output.directory = output.directory, 
	quiet = TRUE
	);

## ---- results = FALSE----------------------------------------------------


## ---- results = FALSE----------------------------------------------------
vcf.specification$job.dependency <- NULL;

run.varitas.pipeline(
    file.details = vcf.specification,
    output.directory = output.directory,
    start.stage = 'annotation',
    quiet = TRUE
    );

## ---- results = FALSE----------------------------------------------------
run.varitas.pipeline(
    file.details = vcf.specification,
    output.directory = output.directory,
    start.stage = 'annotation',
    email = 'Erle.Holgersen@icr.ac.uk',
    quiet = TRUE
    );

## ----wrapper1------------------------------------------------------------
library(varitas)
output.directory <- '.'

fastq.directory <- 'inst/extdata/fastq'
fastq.files <- list.files(
  pattern = 'R1.*\\.fastq', 
  path = fastq.directory, 
  full.names = TRUE
  )
fastq.mate.files <- list.files(
  pattern = 'R2.*\\.fastq', 
  path = fastq.directory, 
  full.names = TRUE
  )

fastq.specification <- data.frame(
  # Extract the sample ID from the filename
  sample.id = gsub('.*Sample0(\\d\\d).*', '\\1', basename(fastq.files)),
  reads = fastq.files,
  mates = fastq.mate.files,
  stringsAsFactors = FALSE
  )

print(fastq.specification)

## ----wrapper2, eval=FALSE, results=FALSE---------------------------------
#  set.varitas.options(filters.vardict.min_tumour_depth = 10)

## ----wrapper3, eval=FALSE, results=FALSE---------------------------------
#  config <- 'inst/extdata/varitas_config.yaml'
#  overwrite.varitas.options(config)

## ----wrapper4, eval=FALSE, results=FALSE---------------------------------
#  run.varitas.pipeline(
#      file.details = fastq.specification,
#      output.directory = output.directory,
#      variant.callers = c('mutect', 'vardict'),
#      quiet = FALSE,
#      run.name = 'EXAMPLE',
#      email = 'adam.mills@icr.ac.uk'
#      )

## ----wrapper5, eval=FALSE, results=FALSE---------------------------------
#  ###############################################################################
#  ## VariTAS Wrapper Script
#  ##
#  ###############################################################################
#  ## Author:
#  ## Adam Mills
#  ###############################################################################
#  ## Libraries:
#  library(varitas)
#  ###############################################################################
#  ## Main
#  
#  output.directory <- '.'
#  
#  fastq.directory <- 'inst/extdata/fastq'
#  fastq.files <- list.files(
#    pattern = 'R1.*\\.fastq',
#    path = fastq.directory,
#    full.names = TRUE
#    )
#  fastq.mate.files <- list.files(
#    pattern = 'R2.*\\.fastq',
#    path = fastq.directory,
#    full.names = TRUE
#    )
#  
#  fastq.specification <- data.frame(
#    sample.id = gsub('.*Sample0(\\d\\d).*', '\\1', basename(fastq.files)),
#    reads = fastq.files,
#    mates = fastq.mate.files,
#    stringsAsFactors = FALSE
#    )
#  
#  config <- 'inst/extdata/varitas_config.yaml'
#  overwrite.varitas.options(config)
#  
#  run.varitas.pipeline(
#    file.details = fastq.specification,
#    output.directory = output.directory,
#    variant.callers = c('mutect', 'vardict'),
#    quiet = FALSE,
#    run.name = 'EXAMPLE',
#    email = 'adam.mills@icr.ac.uk'
#    )
#  

## ----matched1------------------------------------------------------------
fastq.specification <- data.frame(
  sample.id = gsub('.*Sample0(\\d\\d).*', '\\1', basename(fastq.files)),
  patient.id = c('X', 'X', 'Y', 'Y'),
  tissue = c('tumour', 'normal', 'tumour', 'normal'),
  reads = fastq.files,
  mates = fastq.mate.files,
  stringsAsFactors = FALSE
  )

print(fastq.specification)

## ----hybrid1-------------------------------------------------------------
bam.directory <- 'inst/extdata/bam'
bam.files <- list.files(
  pattern = 'Sample.*\\.bam', 
  path = bam.directory, 
  full.names = TRUE
  )
vcf.directory <- 'inst/extdata/vcf'
vcf.files <- list.files(
  pattern = 'Sample.*\\.vcf', 
  path = vcf.directory, 
  full.names = TRUE
  )

bam.specification <- data.frame(
  sample.id = gsub('^Sample_(\\d+).*', '\\1', basename(bam.files)),
  tumour.bam = bam.files,
  stringsAsFactors = FALSE
  )
vcf.specification <- data.frame(
  sample.id = gsub('^Sample_(\\d+).*', '\\1', basename(vcf.files)),
  vcf = vcf.files,
  caller = rep('pgm', length(vcf.files)),
  stringsAsFactors = FALSE
  )
print(bam.specification)
print(vcf.specification)

## ----hybrid2, eval=FALSE, results=FALSE----------------------------------
#  run.varitas.pipeline.hybrid(
#   	bam.specification = bam.specification,
#   	vcf.specification = vcf.specification,
#  	output.directory = 'inst/extdata/output/',
#  	proton = TRUE,
#  	run.name = 'EXAMPLE',
#  	quiet = FALSE,
#   	email = 'adam.mills@icr.ac.uk'
#   	);

## ----hybrid3-------------------------------------------------------------
miniseq.sheet <- 'inst/extdata/miniseq/Example_template.csv'
miniseq.directory <- 'inst/extdata/miniseq'

miniseq.info <- prepare.miniseq.specifications(miniseq.sheet, miniseq.directory)
fastq.specification <- miniseq.info[[ 1 ]]
vcf.specification <- miniseq.info[[ 2 ]]
vcf.specification['caller'] <- rep('miniseq', nrow(vcf.specification))

print(fastq.specification)
print(vcf.specification)

## ----hybrid4, eval=FALSE, results=FALSE----------------------------------
#  run.varitas.pipeline.hybrid(
#   	fastq.specification = fastq.specification,
#   	vcf.specification = vcf.specification,
#  	output.directory = 'inst/extdata/output/',
#  	run.name = 'EXAMPLE',
#  	quiet = FALSE,
#   	email = 'adam.mills@icr.ac.uk'
#   	)

