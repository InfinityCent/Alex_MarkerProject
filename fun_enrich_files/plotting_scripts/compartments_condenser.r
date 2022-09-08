# Functionally Enriched Compartments Condenser
# Alex Daiejavad
# Sept. 06, 2022

#------------------------------------------------------------------------------

setwd('C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted')
directory = 'C:/Users/peree/OneDrive/Desktop/CompBio_Code/fun_enrich_files/to_be_plotted/per_marker_outputs'

library(naturalsort)
library(dplyr)
library(readxl)

compartments = c("cell cortex", "cellular bud", "cellular_component", "chromosome", 
                 "cytoplasm", "cytoplasmic vesicle", "cytoskeleton", 
                 "endomembrane system", "extracellular region", "Golgi apparatus", 
                 "membrane", "microtubule organizing center", 
                 "mitochondrial envelope", "mitochondrion", "nucleolus", "nucleus", "other", 
                 "peroxisome", "ribosome", "site of polarized growth")

nuc = c('chromosome', 'nucleolus', 'nucleus')
growth = c('cellular bud', 'site of polarized growth')
membrane = c('cell cortex', 'membrane')
er_golgi = c('endomembrane system', 'Golgi apparatus')
cytoskeleton = c('cytoskeleton', 'microtubule organizing center')
ves_pero = c('cytoplasmic vesicle', 'peroxisome')
mito = c('mitochondrial envelope', 'mitochondrion')
other = c('other', 'cellular_component', 'cytoplasm', 'extracellular region', 
          'ribosome')

clean_file_name = function(marker, file_name) {
  break_up = stringr::str_split(file_name, "-")[[1]]
  
  transformed = c()
  for (s in break_up) {
    if (startsWith(s, "normal")) {
      nc_num = stringr::str_split(stringr::str_split(s, "of")[[1]][1], "normal")[[1]][2]
      nc_num = paste(marker, "Normal", nc_num)
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

fun_enrich_comps = function(fun_enrich_col) {
  comps = c()
  for (comp in fun_enrich_col) {
    if (comp %in% nuc) {comps = append(comps, 'nuclear')}
    else if (comp %in% growth) {comps = append(comps, 'growth')}
    else if (comp %in% membrane) {comps = append(comps, 'membrane')}
    else if (comp %in% er_golgi) {comps = append(comps, 'ER/Golgi')}
    else if (comp %in% cytoskeleton) {comps = append(comps, 'cytoskeleton')}
    else if (comp %in% ves_pero) {comps = append(comps, 'vesicle/peroxisome')}
    else if (comp %in% mito) {comps = append(comps, 'mitochondria')}
    else if (comp %in% other) {comps = append(comps, 'other')}
  }
  comps = unique(comps)
  comps = paste(comps, collapse = ', ')
  
  return(comps)
}

split_AS = c('with_AreaShape', 'without_AreaShape')

comps_df = data.frame(c('Cluster'), c('Terms'))
for (AS in split_AS) {
  files_path = paste(directory, "/", AS, sep = "")
  markers = list.files(files_path)
  
  for (marker in markers) {
    marker_directory = paste(files_path, "/", marker, sep = "")
    files = naturalsort::naturalsort(list.files(marker_directory))
    
    c_files = c()
    for (f in files) {
      if (stringr::str_split(f, "-")[[1]][2] == "go_slim_c.xlsx") {
        c_files = append(c_files, f)
      }
    }
    
    for (f in c_files) {
      sheet_path = paste(marker_directory, "/", f, sep = "")
      sheet_df = readxl::read_excel(path = sheet_path)
      sheet_df = dplyr::filter(sheet_df, sheet_df$'P-value' < 0.01)
      
      if (nrow(sheet_df) == 0) {next}
      
      cluster_num = clean_file_name(marker, f)[1]
      comps = fun_enrich_comps(sheet_df$Ontology)
      cluster_comps = c(cluster_num, comps)
      
      comps_df = rbind(comps_df, cluster_comps)
    }
  }
}

sheet_names = list('Sheet1' = comps_df)
write.xlsx(sheet_names, file = 'fun_enrich_comps.xlsx')
