#Load file function
#adapted from H Wickham Mastering Shiny, section 9.1.3

load_file <- function(NAME, PATH, SHEET){
  
  ext <- tools::file_ext(NAME)
  switch(ext,
         "xlsx"= readxl::read_excel(PATH, SHEET),
         csv = vroom::vroom(PATH, delim = ",", show_col_types = FALSE),
         tsv = vroom::vroom(PATH, delim = "\t"),
         txt = vroom::vroom(PATH, show_col_types = FALSE),
         validate("Invalid file. Please upload a data file")
        )
    }

#note: for large files DataTable fread is faster for downloading csv files