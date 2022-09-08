# Functional enrichment values compiler
# Alex Daiejavad
# August 23, 2022
# -----------------------------------------------------------------------------

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(stringr)
library(naturalsort)
library(readxl)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyverse)

setwd('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/per_marker_outputs')
dir = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/per_marker_outputs'
splitAS = c('wAS', 'woAS')

compartments = c("cell cortex", "cellular bud", "cellular_component", "chromosome", 
                 "cytoplasm", "cytoplasmic vesicle", "cytoskeleton", 
                 "endomembrane system", "extracellular region", "Golgi apparatus", 
                 "membrane", "microtubule organizing center", 
                 "mitochondrial envelope", "mitochondrion", "nucleolus", "nucleus", "other", 
                 "peroxisome", "ribosome", "site of polarized growth")

f_identifier = function(marker, file_name) {
  s = stringr::str_split(file_name, "-")[[1]][1]
  
  nc_num = stringr::str_split(stringr::str_split(s, "of")[[1]][1], "normal")[[1]][2]
  nc_num = paste("Normal ", nc_num, "/", stringr::str_split(s, "of")[[1]][2], sep = "")
  
  id = paste(marker, ": ", nc_num, sep = "")
  
  return(id)
}

fun_enrich_values = function(comps, fun_enrich_df) {
  values = c()
  i = 1
  
  for (c in comps) {
    if (c %in% fun_enrich_df[, 1]) {
      values = append(values, as.numeric(fun_enrich_df[i, 2]))
      i = i + 1
    } else {
      values = append(values, NA)
    }
  }
  return(values)
}

fun_enrich_table = function(directory, AS) {
  comp_df = data.frame()
  
  marker_dir = paste(directory, "/", AS, sep = "")
  markers = list.files(marker_dir)
  
  for (marker in markers) {
    file_dir = paste(marker_dir, "/", marker, sep = "")
    all_files = naturalsort::naturalsort(list.files(file_dir))
    c_files = c()
    
    for (f in all_files) {
      if (stringr::str_split(f, "-")[[1]][2] == "go_slim_c.xlsx") {
        c_files = append(c_files, f)
      }
    }
    
    for (f in c_files) {
      f_path = paste(file_dir, "/", f, sep = "")
      f_df = readxl::read_excel(path = f_path)
      f_df = f_df[order(f_df$'Ontology'), ]
      f_df = dplyr::filter(f_df, f_df$'P-value' < 0.01)
      
      if (nrow(f_df) == 0) {
        next
      }
      file_id = f_identifier(marker, f)
      f_df_fun_enrich = data.frame(f_df$'Ontology', f_df$'Fold enrichment')
      colnames(f_df_fun_enrich) = c('Ontology', 'Fold Enrichment')
      
      fun_enrich_row = c(file_id, fun_enrich_values(compartments, f_df_fun_enrich))
      comp_df = rbind(comp_df, fun_enrich_row)
    }
  }
  columns = c("Cluster", compartments)
  colnames(comp_df) = columns
  
  return(comp_df)
}

fun_enrich_plots = function(fe_table, AS, plots_path) {
  fe_table_noclusters = fe_table[ , 2:ncol(fe_table)]
  terms = colnames(fe_table_noclusters)
  
  for (i in 1:ncol(fe_table_noclusters)) {
    term = terms[i]
    term_df = data.frame(fe_table[, 1], fe_table_noclusters[, i])
    colnames(term_df) = c('Cluster', term)
    
    term_df = term_df[order(-as.numeric(term_df[,2]), na.last = NA), ]
    term_df[ , 2] = as.numeric(term_df[ , 2])
    
    if (nrow(term_df) == 0) {next}
    if (nrow(term_df) > 20) {term_df = term_df[1:20, ]}
    
    plot_file = paste(term, "_", AS, ".png", sep = "")
    plot_path = paste(plots_path, "/", plot_file, sep = "")
    
    plot_subtitle = ''
    if (AS == 'wAS') {plot_subtitle = 'With AreaShape'}
    if (AS == 'woAS') {plot_subtitle = 'Without AreaShape'}
    
    png(file = plot_path, width = 948, height = 874)
    
    print(ggplot(term_df, aes(x = term_df[ , 2], y = fct_rev(fct_inorder(term_df[ , 1])))) + 
            geom_bar(stat='identity', fill='#FFB516', width = 0.5) +
            labs(x='Functional Enrichment', y='Marker: Cluster', title=term, 
                 subtitle=plot_subtitle) + 
            theme(axis.text = element_text(size=10), 
                  axis.title = element_text(face="bold", size=13), 
                  plot.title = element_text(face="bold", hjust = 0.5, size=15), 
                  plot.subtitle = element_text(hjust = 0.5)))
    
    
    dev.off()
  }
}


wAS = fun_enrich_table(dir, "with_AreaShape")
woAS = fun_enrich_table(dir, "without_AreaShape")

sheet_names = list('with_AreaShape' = wAS, 'without_AreaShape' = woAS)
write.xlsx(sheet_names, file = 'fun_enrich_per_marker_c.xlsx')


wAS_plots_directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/per_marker_outputs/wAS_plots_c'
woAS_plots_directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/per_marker_outputs/woAS_plots_c'
fun_enrich_plots(wAS, "wAS", wAS_plots_directory)
fun_enrich_plots(woAS, "woAS", woAS_plots_directory)



