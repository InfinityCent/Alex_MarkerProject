# Visualizing functional enrichment outputs
# Alex Daiejavad
# August 03, 2022

#==========================================================
if (! requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (! requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (! requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}
if (! requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)

setwd('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted')
directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/per_marker_outputs'
plots_directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/plots/per_marker/'

split_AS = c('with_AreaShape', 'without_AreaShape')

# https://www.geeksforgeeks.org/how-to-read-a-xlsx-file-with-multiple-sheets-in-r/
multiplesheets = function(path_name) {
  sheets = readxl::excel_sheets(path_name)
  tibble = lapply(sheets, function(x) readxl::read_excel(path_name, sheet = x))
  data_frame = lapply(tibble, as.data.frame)
  
  names(data_frame) = sheets
  
  return(data_frame)
}

clean_file_name = function(file_name) {
  break_up = stringr::str_split(file_name, "-")[[1]]
  
  transformed = c()
  for (s in break_up) {
    if (startsWith(s, "normal")) {
      nc_num = stringr::str_split(stringr::str_split(s, "of")[[1]][1], "normal")[[1]][2]
      nc_num = paste("Normal", nc_num, "of", stringr::str_split(s, "of")[[1]][2])
      transformed = c(transformed, nc_num)
    }
    if (s == "go_slim_p.xlsx") {
      transformed = c(transformed, "Biological Process")
    }
    if (s == "go_slim_c.xlsx") {
      transformed = c(transformed, "Cellular Compartment")
    }
  }
  return(transformed)
}

plot_info = function(AS) {
  transformed = c()
  if (AS == "with_AreaShape") {
    transformed = "With AreaShape"
  }
  if (AS == "without_AreaShape") {
    transformed = "Without AreaShape"
  }
  return(transformed)
}

for (AS in split_AS) {
  files_path = paste(directory, "/", AS, sep = "")
  markers = list.files(files_path)
  
  for (marker in markers) {
    marker_directory = paste(files_path, "/", marker, sep = "")
    files = list.files(marker_directory)
    
    for (f in files) {
      sheet_path = paste(marker_directory, "/", f, sep = "")
      sheet_df = readxl::read_excel(path = sheet_path)
      sheet_df = dplyr::filter(sheet_df, sheet_df$'P-value' < 0.01)
      
      if (nrow(sheet_df) == 0) {
        next
      }
      
      sheet_df_subset = data.frame(sheet_df$'Ontology', sheet_df$'Fold enrichment')
      colnames(sheet_df_subset) = c('Term Name', 'Fold Enrichment')
      sheet_df_subset$'Term Name' = factor(stringr::str_wrap(sheet_df_subset$'Term Name', 30))
      
      cluster_num = clean_file_name(f)
      
      plot_title = paste(marker, " (", cluster_num[1], ")", sep = "")
      
      plot_subtitle = plot_info(AS)
      
      plot_color = c()
      if (cluster_num[2] == "Biological Process") {plot_color = '#6586CA'}
      if (cluster_num[2] == "Cellular Compartment") {plot_color = '#FFB516'}
      
      size = c()
      if (cluster_num[2] == "Biological Process") {size = NULL}
      if (cluster_num[2] == "Cellular Compartment") {size = 13}
      
      plot_file = paste(marker, "_", AS, "_", stringr::str_split(f, "\\.")[[1]][1], ".png", sep = "")
      plot_path = paste(plots_directory, plot_file, sep = "")
      
      png(file = plot_path, width = 948, height = 874)
      
      print(ggplot(sheet_df_subset, aes(sheet_df_subset[ , 2], sheet_df_subset[ , 1])) + 
              geom_bar(stat='identity', fill=plot_color, width = 0.5) +
              labs(x='Fold Enrichment', y=cluster_num[2], title=plot_title, subtitle=plot_subtitle) + 
              theme(axis.text = element_text(size=10), 
                    axis.title = element_text(face="bold", size=13), 
                    plot.title = element_text(face="bold", hjust = 0.5, size=15), 
                    plot.subtitle = element_text(hjust = 0.5),
                    axis.text.x = element_text(size = size),
                    axis.text.y = element_text(size = size)))
      dev.off()
    }
  }
}
