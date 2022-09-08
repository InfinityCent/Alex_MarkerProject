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

setwd('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/combined_cluster_outputs')
dir = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/combined_cluster_outputs'
split_CC = c('no_cellcycle', 'with_cellcycle')
splitAS = c('with_AreaShape', 'without_AreaShape')
split_gene = c('all_genes', 'ess_only', 'no_vesicle', 'noness_only')

compartments = c("cell cortex", "cellular bud", "cellular_component", "chromosome", 
                 "cytoplasm", "cytoplasmic vesicle", "cytoskeleton", 
                 "endomembrane system", "extracellular region", "Golgi apparatus", 
                 "membrane", "microtubule organizing center", 
                 "mitochondrial envelope", "mitochondrion", "nucleolus", "nucleus", "other", 
                 "peroxisome", "ribosome", "site of polarized growth")

f_identifier = function(file_name) {
  s = stringr::str_split(file_name, "-")[[1]][1]
  
  nc_num = stringr::str_split(stringr::str_split(s, "of")[[1]][1], "Cluster")[[1]][2]
  nc_num = paste("Cluster ", nc_num, "/", stringr::str_split(s, "of")[[1]][2], sep = "")
  
  return(nc_num)
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

fun_enrich_table = function(root_directory, CC, AS, gene_dir) {
  comp_df = data.frame()
  
  cluster_dir = paste(root_directory, "/", CC, "/", AS, "/", gene_dir, sep = "")
  all_files = naturalsort::naturalsort(list.files(cluster_dir))
  
  c_files = c()
  
  for (f in all_files) {
    if (stringr::str_split(f, "-")[[1]][2] == "go_slim_c.xlsx") {
      c_files = append(c_files, f)
    }
  }
  
  for (f in c_files) {
    f_path = paste(cluster_dir, "/", f, sep = "")
    f_df = readxl::read_excel(path = f_path)
    f_df = f_df[order(f_df$'Ontology'), ]
    f_df = dplyr::filter(f_df, f_df$'P-value' < 0.01)
    
    if (nrow(f_df) == 0) {
      next
    }
    
    file_id = f_identifier(f)
    f_df_fun_enrich = data.frame(f_df$'Ontology', f_df$'Fold enrichment')
    colnames(f_df_fun_enrich) = c('Ontology', 'Fold Enrichment')
    
    fun_enrich_row = c(file_id, fun_enrich_values(compartments, f_df_fun_enrich))
    comp_df = rbind(comp_df, fun_enrich_row)
  }

  columns = c("Cluster", compartments)
  colnames(comp_df) = columns
  
  return(comp_df)
}

fun_enrich_plots = function(fe_table, CC, AS, gene_dir, plots_path) {
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
    
    plot_file = paste(term, "_", CC, "_", AS, "_", gene_dir, ".png", sep = "")
    plot_path = paste(plots_path, "/", plot_file, sep = "")
    
    plot_subtitle = paste(CC, "|", AS, "|", gene_dir)
    
    plot_fill = ''
    if (gene_dir == "All Genes") {plot_fill = '#53EE97'}
    if (gene_dir == "Essential Only") {plot_fill = '#EEBA53'}
    if (gene_dir == "No Vesicle") {plot_fill = '#A07FF0'}
    if (gene_dir == "Nonessential Only") {plot_fill = '#7FB4F0'}
    
    png(file = plot_path, width = 948, height = 874)
    
    print(ggplot(term_df, aes(x = term_df[ , 2], y = fct_rev(fct_inorder(term_df[ , 1])))) + 
            geom_bar(stat='identity', fill=plot_fill, width = 0.5) +
            labs(x='Functional Enrichment', y='Cluster', title=term, 
                 subtitle=plot_subtitle) + 
            theme(axis.text = element_text(size=10), 
                  axis.title = element_text(face="bold", size=13), 
                  plot.title = element_text(face="bold", hjust = 0.5, size=15), 
                  plot.subtitle = element_text(hjust = 0.5)))
    
    
    dev.off()
  }
}


wCC_wAS_allgenes = fun_enrich_table(dir, "with_cellcycle", "with_AreaShape", "all_genes")
wCC_wAS_essonly = fun_enrich_table(dir, "with_cellcycle", "with_AreaShape", "ess_only")
wCC_wAS_novesicle = fun_enrich_table(dir, "with_cellcycle", "with_AreaShape", "no_vesicle")
wCC_wAS_nonessonly = fun_enrich_table(dir, "with_cellcycle", "with_AreaShape", "noness_only")

wCC_woAS_allgenes = fun_enrich_table(dir, "with_cellcycle", "without_AreaShape", "all_genes")
wCC_woAS_essonly = fun_enrich_table(dir, "with_cellcycle", "without_AreaShape", "ess_only")
wCC_woAS_novesicle = fun_enrich_table(dir, "with_cellcycle", "without_AreaShape", "no_vesicle")
wCC_woAS_nonessonly = fun_enrich_table(dir, "with_cellcycle", "without_AreaShape", "noness_only")

woCC_wAS_allgenes = fun_enrich_table(dir, "no_cellcycle", "with_AreaShape", "all_genes")
woCC_wAS_essonly = fun_enrich_table(dir, "no_cellcycle", "with_AreaShape", "ess_only")
woCC_wAS_novesicle = fun_enrich_table(dir, "no_cellcycle", "with_AreaShape", "no_vesicle")
woCC_wAS_nonessonly = fun_enrich_table(dir, "no_cellcycle", "with_AreaShape", "noness_only")

woCC_woAS_allgenes = fun_enrich_table(dir, "no_cellcycle", "without_AreaShape", "all_genes")
woCC_woAS_essonly = fun_enrich_table(dir, "no_cellcycle", "without_AreaShape", "ess_only")
woCC_woAS_novesicle = fun_enrich_table(dir, "no_cellcycle", "without_AreaShape", "no_vesicle")
woCC_woAS_nonessonly = fun_enrich_table(dir, "no_cellcycle", "without_AreaShape", "noness_only")


sheet_names_wCC_wAS = list('all_genes' = wCC_wAS_allgenes, 'ess_only' = wCC_wAS_essonly, 
                           'no_vesicle' = wCC_wAS_novesicle, 'noness_only' = wCC_wAS_nonessonly)

sheet_names_wCC_woAS = list('all_genes' = wCC_woAS_allgenes, 'ess_only' = wCC_woAS_essonly, 
                           'no_vesicle' = wCC_woAS_novesicle, 'noness_only' = wCC_woAS_nonessonly)

sheet_names_woCC_wAS = list('all_genes' = woCC_wAS_allgenes, 'ess_only' = woCC_wAS_essonly, 
                            'no_vesicle' = woCC_wAS_novesicle, 'noness_only' = woCC_wAS_nonessonly)

sheet_names_woCC_woAS = list('all_genes' = woCC_woAS_allgenes, 'ess_only' = woCC_woAS_essonly, 
                            'no_vesicle' = woCC_woAS_novesicle, 'noness_only' = woCC_woAS_nonessonly)

write.xlsx(sheet_names_wCC_wAS, file = 'fun_enrich_cluster_wCCwAS_c.xlsx')
write.xlsx(sheet_names_wCC_woAS, file = 'fun_enrich_cluster_wCCwoAS_c.xlsx')
write.xlsx(sheet_names_woCC_wAS, file = 'fun_enrich_cluster_woCCwAS_c.xlsx')
write.xlsx(sheet_names_woCC_woAS, file = 'fun_enrich_cluster_woCCwoAS_c.xlsx')


wCCwAS_plots = paste(dir, "/", 'wCCwAS_plots_c', sep = "") 
wCCwoAS_plots = paste(dir, "/", 'wCCwoAS_plots_c', sep = "") 
woCCwAS_plots = paste(dir, "/", 'woCCwAS_plots_c', sep = "") 
woCCwoAS_plots = paste(dir, "/", 'woCCwoAS_plots_c', sep = "")

fun_enrich_plots(wCC_wAS_allgenes, "With CellCycle", "With AreaShape", "All Genes", wCCwAS_plots)
fun_enrich_plots(wCC_wAS_essonly, "With CellCycle", "With AreaShape", "Essential Only", wCCwAS_plots)
fun_enrich_plots(wCC_wAS_novesicle, "With CellCycle", "With AreaShape", "No Vesicle", wCCwAS_plots)
fun_enrich_plots(wCC_wAS_nonessonly, "With CellCycle", "With AreaShape", "Nonessential Only", wCCwAS_plots)

fun_enrich_plots(wCC_woAS_allgenes, "With CellCycle", "Without AreaShape", "All Genes", wCCwoAS_plots)
fun_enrich_plots(wCC_woAS_essonly, "With CellCycle", "Without AreaShape", "Essential Only", wCCwoAS_plots)
fun_enrich_plots(wCC_woAS_novesicle, "With CellCycle", "Without AreaShape", "No Vesicle", wCCwoAS_plots)
fun_enrich_plots(wCC_woAS_nonessonly, "With CellCycle", "Without AreaShape", "Nonessential Only", wCCwoAS_plots)

fun_enrich_plots(woCC_wAS_allgenes, "Without CellCycle", "With AreaShape", "All Genes", woCCwAS_plots)
fun_enrich_plots(woCC_wAS_essonly, "Without CellCycle", "With AreaShape", "Essential Only", woCCwAS_plots)
fun_enrich_plots(woCC_wAS_novesicle, "Without CellCycle", "With AreaShape", "No Vesicle", woCCwAS_plots)
fun_enrich_plots(woCC_wAS_nonessonly, "Without CellCycle", "With AreaShape", "Nonessential Only", woCCwAS_plots)

fun_enrich_plots(woCC_woAS_allgenes, "Without CellCycle", "Without AreaShape", "All Genes", woCCwoAS_plots)
fun_enrich_plots(woCC_woAS_essonly, "Without CellCycle", "Without AreaShape", "Essential Only", woCCwoAS_plots)
fun_enrich_plots(woCC_woAS_novesicle, "Without CellCycle", "Without AreaShape", "No Vesicle", woCCwoAS_plots)
fun_enrich_plots(woCC_woAS_nonessonly, "Without CellCycle", "Without AreaShape", "Nonessential Only", woCCwoAS_plots)