T_test = function(data, group1, group2, adj.method="fdr"){
  # the first column of data is identifier
  # the others columns cotain gene expression, samples' name should be colnames
  # group1, group2 are used to compare two groups of data
  # adj.method defines which method used to adjust p values
  
  colnames(data)[1] = "geneSymbol"
  # only get samples we wanna to compare, be sure 
  # all names in group1 and group2 can match colnames of data 
  if (!all(c(group1, group2)%in%colnames(data))){
    stop("Error! The names in group1 and group2 must be in column names of data.")
  }
  data = data[, c("geneSymbol", group1, group2)]
  n1 = length(group1)
  n2 = length(group2)
  # define a function to get p values according to if data rows are var.equal=TURE
  get_pvalue = function(x, n1, n2){
    x1 = as.numeric(x[2:n1+1])
    x2 = as.numeric(x[(n1+2):(n1+n2+1)])
    if_var = var.test(x1, x2)$p.value
    if (if_var > 0.05){
      p_value = t.test(x1, x2, var.equal = T)$p.value
    }else{
      p_value = t.test(x1, x2, var.equal = F)$p.value
    }
    return(p_value)
  }
  pvalue = apply(data, 1, get_pvalue, n1=n1, n2=n2)
  names(pvalue) = data$geneSymbol
  adj_pvalue = p.adjust(pvalue, method = adj.method)
  return(adj_pvalue)
}