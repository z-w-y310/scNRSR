find_hv_genes <- function(count, I, J){
  count <- t(count)
  count_nzero <- lapply(1:I, function(i) setdiff(count[i,], log10(1.01)))
  mu <- sapply(count_nzero, mean)
  mu[is.na(mu)] <- 0
  sd <- sapply(count_nzero, sd)
  sd[is.na(sd)] <- 0
  cv <- sd/mu
  cv[is.na(cv)] <- 0
  high_var_genes <- which(mu >= 1 & cv >= quantile(cv, 0.25))#分位数
  if(length(high_var_genes) < 500) high_var_genes = 1:I
  count_hv <- count[high_var_genes,]
  return(t(count_hv))
}

# Data processing
process <- function(data){
  if(max(data)>100) data <- log2(data+1) 
  J = nrow(data) # number of cell
  I = ncol(data) # number of gene
  count_hv <- find_hv_genes(data, I, J)
  non.zero <- colSums(count_hv != 0)
  count_hv <- count_hv[,non.zero > 0]
  return(count_hv)
}

# Idetifying of gene subspaces
gene_subspace <- function(data, n, Q){
  name <- colnames(data)
  d <- list()
  data_col <- dim(data)[2]
  for(i in 1:data_col){
    a <- density(as.numeric(data[,i]),n)
    d[[i]] <- a$y
  }
  d <- as.data.frame(d)
  names(d) <- name
  gene_label <- as.data.frame(kmeans(t(d), centers=Q, nstart=25, iter.max=2000)$cluster)
  return(gene_label)
}



scNRSR <- function(data, Q, p_low_cut_theta=0, p_high_cut_theta=0.5)
{
  if(is.null(colnames(data))){colnames(data) <- paste0('col',seq(ncol(data)))}
  colnames_save <- colnames(data)
  p_low_cut <- ((1-1/(sqrt(2)-p_low_cut_theta)^2)+0.5)/2
  p_high_cut <- ((1-1/(sqrt(2)-p_high_cut_theta)^2)+0.5)/2
  if( nrow(data) > 10000 ){ n = 2000 }else{ n = 512 }
  
  # Data processing
  data0 <- process(data) # cell X gene
  
  # Identifying of gene subspaces
  gene_subspace <- gene_subspace(data0, n, Q)
  
  # Identifying of cell subsets
  
  result <- purrr::quietly(scDHA)(data0,seed = 1) #,sparse = FALSE
  cell_subset <- result[["result"]][["cluster"]]
  P <- length(unique(cell_subset))
  
  #Data recovery
  imputed <- as.data.frame(data0)
  
  fun1 <- function(x, b_gene, g_gene, data2, dist){
    xx <- t(data2[-x,g_gene])
    yy <- t(data2[x,g_gene])
    ximpute <- t(data2[-x,b_gene])
    
    if(ncol(xx) >= min(200, nrow(xx))){
      filterid = order(dist[x,-x])[1: min(nrow(xx),200)]
      #filterid = colnames(dist[,filterid])
      xx = xx[,filterid, drop = FALSE]
      ximpute = ximpute[, filterid, drop = FALSE]
    }
    if (sum(yy!=0) > length(yy)/20){
      if (nrow(ximpute) > 2){
        nnls_res <- penalized::penalized(response = yy, penalized = xx, unpenalized = ~0,
                                         positive = TRUE, lambda1 = 0, lambda2 = 0, maxiter = 1e+3, trace = TRUE)
        if(length(penalized::coefficients(nnls_res)) > 0){
          yimpute_res <- penalized::predict(nnls_res, penalized = ximpute, unpenalized = ~0)[,1]
          data2[x,b_gene] <- yimpute_res
        }
      }
    }
    y = data2[x,]
    return(y)
  }
  
  cl = makeCluster(12)
  registerDoParallel(cl)
  
  for ( i in 1:P ){
    cell_index <- which(cell_subset == i)
    if (length(cell_index) > 1){ 
      data1 <- imputed[cell_index,]
      if (length(cell_index) > 500){ 
        dd <- data0[cell_index,]
        dd <- irlba::prcomp_irlba(dd, n = 50)$x 
        dist <- dist(dd)
      }else { dist <- dist(data0[cell_index,]) }
      dist <- as.matrix(dist)
      m_data <- apply(data1, 2, mean)
      v_data <- apply(data1, 2, var)
      s_data <- t(apply(data1, 2, function(x) (x - mean(x)) ^ 2))
      p_data <- 1 - (s_data / ((sqrt(2) ^ 2) * v_data))
      p_data[is.na(p_data)] = 0
      dropout_candidate <- ((data1 < m_data) & (t(p_data) < p_low_cut)) |(( data1> m_data) & (t(p_data) < p_high_cut))
      
      
      for ( j in 1 :Q ){
        print(c("i=",i,"j=",j))
        f_gene<- which(gene_subspace == j)
        if (length(f_gene) > 1){
          data2 <- data1[,f_gene]
          #data2 <- as.data.frame(data2)
          dropout_candidate2 <- dropout_candidate[,f_gene]
          im <- foreach(x = 1:nrow(data2), .packages = c("penalized"), .combine = rbind) %dopar% {
            # bad data(noisy data)
            b_gene <- which(dropout_candidate2[x,]=="TRUE") 
            # good data
            g_gene <- which(dropout_candidate2[x,]=="FALSE")
            if (length(b_gene) > 0){
              y = fun1(x, b_gene, g_gene, data2, dist)
            }else { y = data2[x,] }
            y = as.data.frame(y)
            return(y)
          }
          data1[,f_gene] <- im
        }
      }
      imputed[cell_index,] <- data1
    }
  }
  
  stopCluster(cl)
  
  imputed <- as.matrix(imputed) 
  return(imputed)
}

library(parallel)
library(doParallel)
library(scDHA)
imputed <- scNRSR(data,Q=5)



