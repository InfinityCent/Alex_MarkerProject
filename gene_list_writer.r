# SCRIPT FOR COMBINING GENE NAMES FROM MULTIPLE FILES
# JULY 14, 2022
# ALEX DAIEJAVAD

# My first legit script outside of coursework!!

# -----------------------------------------------------------------------------
if (! requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}

if (! requireNamespace("naturalsort", quietly = TRUE)) {
  install.packages("naturalsort")
}

library(plyr)
library(naturalsort)

GENE_LIST_FOLDER = "C:/Users/peree/OneDrive/Desktop/CompBio/gene_lists"
markers = c('cdc11', 'dad2', 'heh2', 'hta2', 'nop10', 'nuf2', 'om45')

for (marker in markers) {
  setwd("C:/Users/peree/OneDrive/Desktop/CompBio")
  MARKER_NAME = marker
  ZIPPED_FOLDER = paste(marker, ".zip", sep = "")
  OUTPUT_FOLDER = paste("C:/Users/peree/OneDrive/Desktop/CompBio/", marker, sep = "")
  
  # Unzip file contents into new folder
  unzip(ZIPPED_FOLDER, exdir = OUTPUT_FOLDER)
  
  # Get all csv files in unzipped folder and read them
  csv_files = list.files(path = OUTPUT_FOLDER, pattern = "*.csv")
  csv_files = naturalsort(csv_files)
  setwd(OUTPUT_FOLDER)
  
  # Create list containing each gene names column
  gene_names_list = list()
  for (item in csv_files) {
    csv_file = read.csv(item)
    names_list = data.frame(csv_file$Name)
    gene_names_list = append(gene_names_list, names_list)
  }
  
  # Find length of longest column
  longest_col = max(sapply(gene_names_list, length))
  
  # Combine gene names from different clusters into one dataframe
  # R dataframes need to have columns with equal number of rows, so shorter columns
  # will be made equal to longer columns by doing NA padding
  gene_names_list = lapply(gene_names_list, function(x) { x[1:longest_col] })
  gene_names = data.frame(gene_names_list)
  
  headers = c()
  for (i in 0:(ncol(gene_names) - 1)) {
    header = paste("Cluster", i)
    headers = c(headers, header)
  }
  
  colnames(gene_names) = headers
  
  # And export!!
  setwd(GENE_LIST_FOLDER)
  write.csv(gene_names, file = paste(MARKER_NAME, ".csv", sep = ""), row.names = FALSE)
}
