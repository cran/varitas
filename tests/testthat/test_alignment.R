library(varitas);


context('Run alignment on a single sample');

test_that(
	'Valid FASTQ files supplied', {

		expect_error(
			run.alignment.sample(
				fastq.files = rep('A', 3), 
				sample.id = 'TEST', 
				job.dependencies = 'TEST',
				quiet = TRUE
				),
			'Cannot accept more than two FASTQ files'
			);

		expect_error(
			run.alignment.sample(
				fastq.files = c(), 
				sample.id = 'TEST', 
				job.dependencies = 'TEST', 
				quiet = TRUE
				),
			'Must supply at least one FASTQ file'
			);

		expect_error(
			run.alignment.sample(
				fastq.files = c('A', 'B'), 
				sample.id = 'TEST' 
				),
			'No job dependency supplied, yet FASTQ files do not exist'
			);
				
		});
	

test_that(
	'Valid sample ID supplied', {

		expect_error(
			run.alignment.sample(
				fastq.files = 'TEST', 
				sample.id = c('one', 'two'), 
				job.dependencies = 'TEST',,
				quiet = TRUE
				), 
			'sample.id must have length 1'
			);

		expect_error(
			run.alignment.sample(
				fastq.files = 'TEST', 
				sample.id = 42,
				job.dependencies = 'TEST', 
				quiet = TRUE
				), 
			'sample.id is not a character vector'
			);

	});


# test_that(
# 	'Minimal run works', {

# 		alignment.script <- system.file('perl', 'run_alignment.pl');
# 		expected.message <- paste(
# 			'perl', 
# 			alignment.script, 
# 			'--fastq FASTQ1 FASTQ2 --sample_id SAMPLE_ID --job_dependency JOB_DEPENDENCY'
# 			);

# 		expect_output(
# 			run.alignment.sample(
# 				fastq.files = c('FASTQ1', 'FASTQ2'),
# 				sample.id = 'SAMPLE_ID', 
# 				job.dependencies = 'JOB_DEPENDENCY', 
# 				quiet = TRUE
# 				),
# 			expected.message
# 			);

# 	});