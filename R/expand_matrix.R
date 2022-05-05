expand_matrix = function(m){
  un1 = unique(sort(c(colnames(m), rownames(m))))
  m2 = matrix(0, length(un1), length(un1), dimnames = list(un1, un1))
  m2[row.names(m), colnames(m)] = m
  m2
}