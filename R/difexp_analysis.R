difexp_analysis = function(feature_counts,
                           exp_name=format(Sys.Date(),"%Y_%m_%d_"),
                           sample_names,
                           sample_groups,
                           design,
                           my_contrasts,
                           RPKM_cutoff = 1,
                           FDRts = 0.05,
                           FCts = log2(2)) {
  # Split feature_counts data to create DGEL object -------------------------
  counts = feature_counts %>% dplyr::select(-c(Geneid,Chr,Start,End,Strand,Length)) %>% `colnames<-`(sample_names)
  gene_info = feature_counts %>% dplyr::select(c(Geneid,Chr,Start,End,Strand,Length))

  # Create DGEL and RPKM objects --------------------------------------------
  print("Generating DGEL object for edgeR and RPKM object")
  DGEL_object = DGEList(counts = counts, 
                        genes = gene_info,
                        group = sample_groups) %>% calcNormFactors()
  RPKM_object = rpkm(DGEL_object, gene.length = "Length") %>% as.data.frame %>% cbind("gene_id"=DGEL_object$genes$Geneid,.)

  # Remove low expression genes from DGEL and estimate dispersion -----------
  print("Performing statistical tests")
  drop = DGEL_object %>% rpkm %>% rowMeans %>% `<`(RPKM_cutoff) %>% which #get index of genes below threshold
  DGEL_object_filtered = DGEL_object[-drop, ,keep.lib.sizes=FALSE] %>% calcNormFactors %>% estimateDisp(design, robust=TRUE)
  
  # Statistical test by my_comparisons, filter by FDR and logFC thresholds and merge results --------------------
  # TO DO change loop to lapply (vectorize contrasts?)
  fit = glmQLFit(DGEL_object_filtered, design, robust=TRUE)
  
  test_list = list()
  for(contrast in colnames(my_contrasts)) {
    print(contrast)
    test_list[[contrast]] = glmQLFTest(fit, contrast=my_contrasts[,contrast]) %>% topTags(n = Inf) %>% as.data.frame() %>% 
      rename_with( ~ paste(.x, contrast, sep = "_"), .cols = c("logFC","logCPM","F","PValue","FDR"))
  }
  
  merged_tests = test_list %>% purrr::reduce(full_join, map2, by = c("Geneid","Chr","Start","End","Strand","Length"))
 
  # Filter results by FDR and logFC thresholds ------------------------------
  print(paste0("Filtering resutls by FDR lower than ",FDRts," and abs(log2FC) higher than ",FCts))

  fdr_colnames = sapply(test_list, function(x) colnames(x)[grepl("^FDR",colnames(x))])
  contrast_colnames = sapply(test_list, function(x) colnames(x)[grepl("^logFC",colnames(x))])
  
  filter_func = function(x, fdr_col, fc_col) {
    dplyr::filter(x, (get(fdr_col) < FDRts) & (abs(get(fc_col)) > FCts))
  }
  
  filtered_test_list = mapply(filter_func,
                              test_list,
                              fdr_colnames,
                              contrast_colnames,
                              SIMPLIFY = FALSE)

  filtered_tests = filtered_test_list %>% purrr::reduce(full_join, map2, by = c("Geneid","Chr","Start","End","Strand","Length"))
  
  # Add gene name and description to output file and remove extra columns
  extra_columns = colnames(merged_tests)[grepl(c("^logCPM|^F_|^PValue"),colnames(merged_tests))]
  filtered_tests_output = filtered_tests %>% add_gene_annotation %>% dplyr::select(-all_of(extra_columns))

  # Return list with generated objects --------------------------------------
  result_list = list("rpkm" = RPKM_object,
                     "DGEL" = DGEL_object_filtered,
                     "list_by_contrast" = filtered_test_list,
                     "merged_tests" = merged_tests,
                     "filtered_tests" = filtered_tests_output)
  return(result_list)
}