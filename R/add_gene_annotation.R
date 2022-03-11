add_gene_annotation = function(edgeR_output,gene_id_column="Geneid") {
  # Load prerequisite libraries and generate ensemble plants biomaRt object 
  library(org.At.tair.db)
  library(biomaRt)
  library(tidyverse)
  if(!exists("ensembl")){
    ensembl = useMart("plants_mart",host="https://plants.ensembl.org")
    ensembl = useDataset("athaliana_eg_gene",mart=ensembl)
  }
  edgeR_output = as.data.frame(edgeR_output)
  gene_list = edgeR_output %>% pull(gene_id_column)
  print(paste0(length(gene_list)," genes detected"))
  
  # Get gene symbol from org.At.tair.db (separated by comma if there are multiple symbols)
  print("Retrieving gene symbol")
  gene_symbol = AnnotationDbi::mapIds(org.At.tair.db,
                                      column = "SYMBOL",
                                      keys = gene_list,
                                      keytype = "TAIR",
                                      multiVals = "list") %>%
    map_chr(., ~str_c(.x, collapse=", "))

  # Get long gene description from org.At.tair.db (only one description is retrieved)
  print("Retrieving long gene description")
  long_description = AnnotationDbi::mapIds(org.At.tair.db,
                                           column = "GENENAME",
                                           keys = gene_list,
                                           keytype = "TAIR",
                                           multiVals = "first")
  
  # Get short gene description from ensemble plants and remove [Source:IDXXX]
  print("Retrieving short gene description")
  short_description = getBM(attributes = c("tair_locus","description"),
                            filters = "tair_locus",
                            values = gene_list,
                            mart = ensembl) %>%
    dplyr::rename(short_description = description) %>%
    mutate(short_description = gsub(" \\[Source.*","",short_description))
  
  # Write output table, if no gene symbol, add gene id
  print("Writing output table")
  annotated_edgeR_output = edgeR_output %>%
    mutate(gene_symbol = coalesce(gene_symbol,Geneid), .after=Geneid,
           long_description = long_description) %>%
    left_join(short_description, by = c("Geneid" = "tair_locus"),) %>%
    relocate(short_description, .after = gene_symbol) %>%
    dplyr::select(-c(Chr,Start,End,Strand,Length))

  # Return output
  return(annotated_edgeR_output)
}
