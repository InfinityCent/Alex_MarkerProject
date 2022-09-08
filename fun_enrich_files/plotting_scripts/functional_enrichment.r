# FUNCTIONAL ENCRICHMENT SCRIPT
# JULY 18, 2022
# ALEX DAIEJAVAD

# -----------------------------------------------------------------------------

if (! requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (! requireNamespace("gprofiler2", quietly = TRUE)) {
  install.packages("gprofiler2")
}

if (! requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (! requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
BiocManager::install(c("clusterProfiler", "enrichplot", "DOSE"))

library(dplyr)
library(gprofiler2)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(DOSE)

setwd("C:/Users/peree/OneDrive/Desktop/CompBio/gene_lists")
DIRECTORY = "C:/Users/peree/OneDrive/Desktop/CompBio/gene_lists"
MARKERS = list.files(DIRECTORY)

for (marker in MARKERS) {
  path = paste(DIRECTORY, "/", marker, sep = "")
  clusters = read.csv(file = path)
  
  for (i in 1:ncol(clusters)) {
    enrich_cluster = gprofiler2::gost(query = clusters[,i], organism = 'scerevisiae', 
                                      significant = TRUE)$result
    enrich_cluster = dplyr::filter(enrich_cluster, source %in% 'GO:MF') # only care about GO molecular function
    
    enrich_cluster_results = data.frame(enrich_cluster$term_name, enrich_cluster$precision)
    colnames(enrich_cluster_results) = c('Term Name', 'Proportion')
    enrich_cluster_results = enrich_cluster_results[order(enrich_cluster_results$Proportion, decreasing = TRUE), ]
    
    cluster_data = enrich_cluster_results[1:20, ]
    cluster_data$'Term Name' = factor(stringr::str_wrap(cluster_data$'Term Name', 39), 
                                   levels=stringr::str_wrap(cluster_data$'Term Name', 39))
    
    plot_title = paste("Top 20 Gene Ontology Molecular Functions (Cluster ", i-1, ")", sep = "")
    
    print(ggplot(cluster_data, aes(cluster_data[ , 2], cluster_data[ , 1])) + 
            geom_bar(stat = 'identity') +
            labs(x = 'Proportion', y = 'Molecular Function', 
                 title = plot_title))
  }
}

for (marker in MARKERS) {
  path = paste(DIRECTORY, "/", marker, sep = "")
  clusters = read.csv(file = path)
  
  for (i in 1:ncol(clusters)) {
    enrich_cluster = gprofiler2::gost(query = clusters[,i], organism = 'scerevisiae', 
                                      significant = TRUE)$result
    enrich_cluster = dplyr::filter(enrich_cluster, source %in% 'GO:MF') # only care about GO molecular function
    
    enrich_cluster_results = data.frame(enrich_cluster$term_name, enrich_cluster$precision)
    colnames(enrich_cluster_results) = c('Term Name', 'Proportion')
    enrich_cluster_results = enrich_cluster_results[order(enrich_cluster_results$Proportion, decreasing = TRUE), ]
    
    cluster_data = enrich_cluster_results[1:20, ]
    cluster_data$'Term Name' = factor(stringr::str_wrap(cluster_data$'Term Name', 39), 
                                      levels=stringr::str_wrap(cluster_data$'Term Name', 39))
    
    plot_title = paste("Top 20 Gene Ontology Molecular Functions (Cluster ", i-1, ")", sep = "")
    
    print(ggplot(cluster_data, aes(cluster_data[ , 2], cluster_data[ , 1])) + 
            geom_bar(stat = 'identity') +
            labs(x = 'Proportion', y = 'Molecular Function', 
                 title = plot_title))
  }
}

