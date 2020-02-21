#' Save coverage statistics to multi-worksheet Excel file.
#' 
#' @param project.directory Path to project directory
#' @param file.name Name of output file  
#' @param overwrite Logical indicating whether to overwrite existing file if it exists.
#'
#' @return None
#' 
#'
save.coverage.excel <- function(project.directory, file.name, overwrite = TRUE) {
    
    ### INPUT TESTS ###########################################################
    
    assertthat::assert_that( is.character(project.directory) );
    assertthat::assert_that( 1 == length(project.directory) );
    assertthat::assert_that( dir.exists(project.directory) );
    
    ### MAIN ##################################################################
    
    # get data
    coverage.sample <- get.coverage.by.sample.statistics(project.directory);
    coverage.amplicon <- get.coverage.by.amplicon(project.directory);
    
    # define header style for bolding
    header.style <- openxlsx::createStyle(textDecoration = 'Bold');
    
    workbook <- openxlsx::createWorkbook();
    
    # add coverage by sample
    sheet.name <- 'Coverage by sample';
    
    openxlsx::addWorksheet(workbook, sheetName = sheet.name);
    openxlsx::setColWidths(
        workbook, 
        sheet = sheet.name, 
        cols = 1:ncol(coverage.sample), 
        widths = pmax(10, nchar(names(coverage.sample)) + 1)
    );
    
    openxlsx::writeData(
        workbook, 
        sheet = sheet.name, 
        coverage.sample, 
        headerStyle = header.style
    );	
    
    # add coverage by amplicon
    sheet.name <- 'Coverage by amplicon';
    
    openxlsx::addWorksheet(workbook, sheetName = sheet.name);
    openxlsx::setColWidths(
        workbook, 
        sheet = sheet.name, 
        cols = 1:ncol(coverage.amplicon), 
        widths = pmax(10, nchar(names(coverage.amplicon)) + 1)
    );
    
    openxlsx::writeData(
        workbook, 
        sheet = sheet.name, 
        coverage.amplicon, 
        headerStyle = header.style
    );	
    
    ### SAVE TO FILE
    
    openxlsx::saveWorkbook(
        workbook, 
        file.name,
        overwrite = overwrite
    );
    
}
