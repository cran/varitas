#' Save variants to Excel.
#'
#' @description
#'	Makes an Excel workbook with variant calls. If filters are provided, these will 
#'	be saved to an additional worksheet within the same file.
#'
#' @param variants 
#' 	Data frame containing variants
#' @param file.name 
#'	Name of output file
#' @param filters 
#' 	Optional list of filters to be saved
#' @param overwrite 
#'	Logical indicating whether to overwrite exiting file if it exists. Defaults to TRUE for consistency with other R functions.
#'
#'
save.variants.excel <- function(
    variants, 
    file.name,
    filters = NULL, 
    overwrite = TRUE
) { 
    
    ### INPUT TESTS ###########################################################
    
    assertthat::assert_that( is.data.frame(variants) );
    
    ### MAIN ##################################################################
    
    # define header style for bolding
    header.style <- openxlsx::createStyle(textDecoration = 'Bold');
    
    workbook <- openxlsx::createWorkbook();
    
    ## ADD VARIANTS
    sheet.name <- 'Variants';
    
    openxlsx::addWorksheet(workbook, sheetName = sheet.name);
    openxlsx::setColWidths(
        workbook, 
        sheet = sheet.name, 
        cols = 1:ncol(variants), 
        widths = pmax(10, nchar(names(variants)) + 1)
    );
    
    openxlsx::writeData(
        workbook, 
        sheet = sheet.name, 
        variants, 
        headerStyle = header.style
    );	
    
    ## ADD FILTERS (IF REQUESTED)
    
    if( !is.null(filters) ) { 
        sheet.name <- 'Filters';
        
        # force conversion to character vector
        filters$DUMMY <- 'DUMMY';
        filter.vector <- unlist(filters);
        filter.vector <- filter.vector[ 'DUMMY' != filter.vector ];
        
        filter.data <- data.frame(
            'Filter' = names(filter.vector), 
            'Value' = filter.vector, 
            stringsAsFactors = FALSE
        );
        
        # add to workbook
        sheet.name <- 'Filters';
        
        openxlsx::addWorksheet(workbook, sheetName = sheet.name);
        openxlsx::setColWidths(
            workbook, 
            sheet = sheet.name, 
            cols = 1:2, 
            widths = c(max(nchar(filter.data$Filter)), 10)
        );
        
        openxlsx::writeData(
            workbook, 
            sheet = sheet.name, 
            filter.data, 
            headerStyle = header.style
        );	
    }
    
    ### SAVE TO FILE
    
    openxlsx::saveWorkbook(
        workbook, 
        file.name,
        overwrite = overwrite
    );
    
}