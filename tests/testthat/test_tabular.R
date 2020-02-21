library(varitas);

context("Median from tabular data")

test_that(
	"Median of constant vector is constant value", {
		expect_equal(
			tabular.median(
				values = rep(5, 10),
				frequencies = rep(1, 10)
				),
				5
			);
			
		expect_equal(
			tabular.median(
				values = rep(20, 100),
				frequencies = rep(1, 100)
				),
			20
			);
			
		expect_equal(
			tabular.median(
				values = 1,
				frequencies = 10
				),
			1
			);
				
		});
	
