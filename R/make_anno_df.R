make_anno_df = function(de_list,
                        peak_list,
                        target_region = NULL,
                        gene_id_column = "Geneid",
                        peak_name_column = "name")
  {
  de_gene_ids = sapply(de_list, dplyr::pull, get(gene_id_column)) %>% unlist %>% unique
  
  direct_target_ids = mapply(get_direct_targets, peak_list, de_list, 
                             MoreArgs = list(target_region = target_region,
                                             output = "vector",
                                             peak_name_column = peak_name_column), 
                             SIMPLIFY = FALSE)

  anno_df = lapply(direct_target_ids, 
                   function(x,y) {as.integer(x %in% y)}, 
                   x = de_gene_ids) %>% 
    as.data.frame(row.names = de_gene_ids)
  
  return(anno_df)
  }