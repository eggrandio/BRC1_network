get_coex_data = function(gene_ids,
                         data_dir = "D:/ATTED/Ath-m.v21-01.G20819-S12686.combat_pca_subagging.ls.d/",
                         use_biomart = FALSE,
                         as_matrix = TRUE) {
  
  # Convert TAIR IDs to Entrez IDs (files with coexpression data are named by Entrez ID)
  # Sort gene_ids and convert to Entrez ID
  gene_ids = sort(gene_ids)
  
  entrez_ids = AnnotationDbi::select(org.At.tair.db,
                                     keys = gene_ids,
                                     columns = c("ENTREZID"),
                                     keytype = "TAIR") %>% suppressMessages() 
  
  # Sometimes BiomaRt can retrieve additional Entrez IDs
  if(use_biomart){
    ensembl = biomaRt::useMart(biomart = "plants_mart",
                               dataset = "athaliana_eg_gene",
                               host = "https://plants.ensembl.org")
    
    biomart_entrez_ids = getBM(attributes = c("tair_locus", "entrezgene_id"),
                               filters = "tair_locus",
                               values = gene_ids,
                               mart = ensembl) %>% 
      group_by(tair_locus) %>% dplyr::slice(1) %>% ungroup %>% 
      mutate_if(is.integer, as.character) %>% 
      dplyr::rename(TAIR = tair_locus)
    
    entrez_ids = left_join(entrez_ids, biomart_entrez_ids, by = "TAIR") %>% 
      arrange(TAIR) %>% 
      mutate(ENTREZID = coalesce(ENTREZID, entrezgene_id))
  }
  
  # Check if all input genes have an Entrez ID
  missing_genes = is.na(entrez_ids$ENTREZID) %>% sum
  #entrez_ids$TAIR[is.na(entrez_ids$ENTREZID)]
  if(missing_genes > 0) {
    print(paste0(missing_genes, " AGI IDs out of ", length(gene_ids)," do not have a corresponding ENTREZ ID"))
  }
  
  # Check if there is coex data available for each Entrez ID
  # Get list of files with coexpression data for each gene
  data_files = list.files(data_dir)
  print(paste0("Coexpression data available for ", length(data_files)," genes"))
  
  if(!all(entrez_ids$ENTREZID %in% data_files)) {
    missing_data = gene_ids[!(entrez_ids$ENTREZID %in% data_files)]
    cat("There is no coexpression data for", length(missing_data), "gene(s):\n", paste0(missing_data,"\n"))
  }
  
  # Read files with coexpression data info and merge them into a data.frame
  ids_with_coex = entrez_ids$ENTREZID[entrez_ids$ENTREZID %in% data_files] %>% sort()
  
  coex_data_list = sapply(ids_with_coex, function(x) { fread(paste0(data_dir, x), 
                                                             col.names = c("gene_id", x),
                                                             key = "gene_id") },
                          simplify = FALSE, 
                          USE.NAMES = TRUE)
  
  # Merge data, filter coexpression data for input genes only and reorder columns
  agi_ids = mapIds(org.At.tair.db,
                   keys = ids_with_coex,
                   column = "TAIR",
                   keytype = "ENTREZID") %>% suppressMessages()
  
  coex_data = purrr::reduce(coex_data_list, dplyr::full_join, by = "gene_id") %>% 
    .[gene_id %in% ids_with_coex, ] %>% 
    dplyr::mutate(gene_id = agi_ids[as.character(gene_id)]) %>% 
    `colnames<-`(agi_ids[colnames(.)]) %>% 
    setcolorder(., c("gene_id", sort(colnames(.)[-1]))) %>% 
    setkey(., NULL) %>% setkey(., gene_id) %>% # if data.table key is not reassigned, gene_id column is not sorted (?)
    {if(as_matrix) column_to_rownames(., "gene_id") %>% as.matrix %>% `diag<-`(0) else .}
  
  return(coex_data)
  
}