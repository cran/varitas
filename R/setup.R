# Set config.yaml as VariTAS options
.onLoad <- function(libname, pkgname) {

	# add package name to enable loading it later
	# can't use add.option for this since we haven't added any varitas options yet!
	options(varitas = list(pkgname = pkgname));

	# set options according to default config file
	config.file <- system.file('config.yaml', package = pkgname);

	overwrite.varitas.options(config.file);
}
