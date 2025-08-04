rm(list=ls())
options(stringsAsFactors = F)
setwd("./")
getwd() 


Gaussian_Kernel <- function(x)
{
  euclidean_dist <- as.matrix(dist(x, method = "euclidean"))
  median_dist <- median(euclidean_dist[upper.tri(euclidean_dist)])
  gamma <- 1 / (2 * median_dist^2)
  rbf_kernel <- exp(-gamma * euclidean_dist^2)
  return(rbf_kernel)
}
g_Gaussian <- Gaussian_Kernel(data_dropout)
library(proxy)
g_cosine <- proxy::simil(data_dropout, method = "cosine")
g_cosine <- as.matrix(g_cosine)
diag(g_cosine) <- 1
g_spearman <- cor(t(data_dropout), method = "spearman")

c_Gaussian <- Gaussian_Kernel(t(data_dropout))
c_cosine <- proxy::simil(t(data_dropout), method = "cosine")
c_cosine <- as.matrix(c_cosine)
diag(c_cosine) <- 1
c_spearman <- cor(data_dropout, method = "spearman")

K <- 20;
alpha <- 0.5; 
T <- 20;
dif_threshold <- 1e-7;
p_normalize <- function(X)
{
  row.sum.mdiag <- rowSums(X) - diag(X)
  row.sum.mdiag[row.sum.mdiag == 0] <- 1
  X <- X/(2 * (row.sum.mdiag))
  diag(X) <- 0.5
  return(X)
}
s_normalize <- function (xx, KK = 20) 
{
  zero <- function(x) {
    s = sort(x, index.return = TRUE)
    x[s$ix[1:(length(x) - KK)]] = 0
    return(x)
  }
  normalize <- function(X) X/rowSums(X)
  A = matrix(0, nrow(xx), ncol(xx))
  for (i in 1:nrow(xx)) {
    A[i, ] = zero(xx[i, ])
  }
  return(normalize(A))
}
p_update <- function (p,s)
{
  it <- 0
  dif <- 1
  LW <- length(p)
  nextW <- vector("list", LW)
  while(dif > dif_threshold) 
  {
    p_low <- matrix(0, nrow(p[[1]]), ncol(p[[1]]))
    for (i in 1:LW) {
      p_low <- p_low + p[[i]]
    }
    p_low <- p_low/LW
    #p_low <- p_normalize(p_low)
    p_low <- (p_low + t(p_low))/2
    for (j in 1:LW) {
      sumWJ <- matrix(0, nrow(p[[j]]), ncol(p[[j]]))
      for (k in 1:LW) {
        if (k != j) {
          sumWJ <- sumWJ + p[[k]]
        }
      }
      nextW[[j]] <- s[[j]] %*% (sumWJ/(LW - 1)) %*% t(s[[j]])
    }
    for (j in 1:LW) {
      p[[j]] <- p_normalize(nextW[[j]])
      #p[[j]] <- (p[[j]] + t(p[[j]]))/2
    }
    it <- it + 1
    p_new <- matrix(0, nrow(p[[1]]), ncol(p[[1]]))
    for (i in 1:LW) {
      p_new <- p_new + p[[i]]
    }
    p_new <- p_new/LW
    #p_new <- p_normalize(p_new)
    p_new <- (p_new + t(p_new))/2
    
    dif <- norm(p_new-p_low,type = "F")/norm(p_low,type = "F")
    cat("<<<Number of iterations: ", it, "\n")
    cat("<<<dif is: ",dif,"\n")
  }
  return(p_new)
}
g_dist <- list(g_Gaussian,g_cosine,g_spearman)
LW <- length(g_dist)
g_p <- vector("list", LW)
for (i in 1:LW) {
  g_p[[i]] <- p_normalize(g_dist[[i]])
  #g_p[[i]] <- (g_p[[i]] + t(g_p[[i]]))/2
}
g_s <- vector("list", LW)
for (i in 1:LW) {
  g_s[[i]] <- (s_normalize(g_dist[[i]], K))
}
g_pnew <- p_update(g_p,g_s)
c_dist <- list(c_Gaussian,c_cosine,c_spearman)
LW <- length(c_dist)

c_p <- vector("list", LW)
for (i in 1:LW) {
  c_p[[i]] <- p_normalize(c_dist[[i]])
  #g_p[[i]] <- (g_p[[i]] + t(g_p[[i]]))/2
}
c_s <- vector("list", LW)
for (i in 1:LW) {
  c_s[[i]] <- (s_normalize(c_dist[[i]], K))
}
c_pnew <- p_update(c_p,c_s)



library(Seurat)
library(FNN)
library(Matrix)
prepare_data <- function(count_matrix) {
  counts <- as.matrix(count_matrix)
  rownames(counts) <- iconv(rownames(counts), to = "ASCII//TRANSLIT", sub = " ")
  rownames(counts) <- gsub("|", "-", rownames(counts), fixed = TRUE)
  rownames(counts) <- gsub(" +", "-", rownames(counts))
  rownames(counts)[rownames(counts) == ""] <- "UnknownGene" 
  rownames(counts) <- make.unique(rownames(counts))
  counts <- counts[!is.na(rownames(counts)), ]
  sce <- CreateSeuratObject(counts = counts)
  sce <- NormalizeData(sce)   
  sce <- ScaleData(sce)       
  sce <- FindVariableFeatures(object = sce, selection.method = "vst", nfeatures = 2000)
  sce <- RunPCA(sce, npcs=20)  
  return(sce)
}
detect_technical_zeros <- function(sce, k=20, quantile_thresh) {
  norm_data <- as.matrix(GetAssayData(sce, slot = "data"))
  counts <- as.matrix(GetAssayData(sce, slot = "counts"))
  pca_emb <- as.matrix(Embeddings(sce, "pca"))

  knn <- get.knn(pca_emb, k=k)
  knn_index <- knn$nn.index

  technical_zeros <- matrix(FALSE, nrow(counts), ncol(counts))
  rownames(technical_zeros) <- rownames(counts)
  colnames(technical_zeros) <- colnames(counts)

  gene_activity <- rowMeans(norm_data > 0)

  for(gene_idx in 1:nrow(counts)) {
    zero_cells <- which(counts[gene_idx, ] == 0)
    if(length(zero_cells) == 0) next
    gene_active <- gene_activity[gene_idx] > quantile(gene_activity, quantile_thresh)
    neighbor_expr <- sapply(zero_cells, function(cell) {
      neighbors <- knn_index[cell, ]
      mean(norm_data[gene_idx, neighbors] > 0)
    })
    tech_flags <- gene_active & (neighbor_expr > quantile_thresh)
    technical_zeros[gene_idx, zero_cells] <- tech_flags
  }
  
  return(technical_zeros)
}
count_matrix <- data_dropout
count_matrix <- as.matrix(count_matrix)
sce <- prepare_data(count_matrix)
tech_zeros <- detect_technical_zeros(sce, k=20, quantile_thresh=0.2)
PP1 <- matrix(as.numeric(tech_zeros), nrow = nrow(tech_zeros))
M <- 1-PP1

r <- 400
nIter <- 101
error_threshold <- 1e-5
parameter <- c(10,10,1,1)
lambda1 <- parameter[1]
lambda2 <- parameter[2]
lambda3 <- parameter[3]
lambda4 <- parameter[4]
lambda5 <- 1e-7
k_n <- length(unique(data_label))
X <- as.matrix(data_dropout)
m <- nrow(X)
n <- ncol(X)
row_sum <- rowSums(X)

W <- matrix(runif(m*r,min=0,max=1),nrow=m,ncol=r)
H <- matrix(runif(r*n,min=0,max=1),nrow=r,ncol=n)
α <- matrix(runif(m*r,min=0,max=1),nrow=m,ncol=r)
β <- matrix(runif(r*n,min=0,max=1),nrow=r,ncol=n)
WH0 <- matrix(1,m,n)
X_new <- matrix(1,m,n)
iter <- 1
error <- 1
obj <- numeric(nIter)

while((iter < nIter) & (error > error_threshold)){
  cat("<<<Number of iterations: ", iter, "\n")
  w1 <- (-2) * (M * X) %*% t(H) + 2 * (M * (W %*% H)) %*% t(H) + 2 * lambda1 * Ug %*% W + 2 * lambda3 * W - α
  H1 <- (-2) * t(W) %*% (M * X) + 2 * t(W) %*% (M * (W %*% H)) + 2 * lambda2 * H %*% Uc + 2 * lambda4 * H - β
  α <- α + lambda5 * (-W)
  α[α < 0] <- 0
  β <- β + lambda5 * (-H)
  β[β < 0] <- 0
  W <- W - lambda5 * w1
  H <- H - lambda5 * H1
  
  Y <- W %*% H
  Q <- PP1 * Y
  X_new <- X + Q 
  X_list <- list(W,H,Y,Q,X_new)

  #compute object function
  obj1 <- obj1 <- norm(M * (X - (W %*% H)),type = "F")^2
  obj2 <- lambda1 * sum(diag(t(W) %*% Ug %*% W))
  obj3 <- lambda2 * sum(diag(H %*% Uc %*% t(H)))
  obj4 <- lambda3 * norm(W,type = "F")^2
  obj5 <- lambda4 * norm(H,type = "F")^2
  obj[iter] <- obj1 + obj2 + obj3 + obj4 + obj5
  cat("<<<objective function: ",obj[iter],"\n")
  error <- norm(Y - WH0,type = "F") / norm(WH0,type = "F")
  cat("<<<error is: ",error,"\n")
  iter <- iter+1
  WH0 <- Y
  cat("####################################","\n")
}

library(cluster)
library(tidyverse)

evaluation <- function(truelabel, predlabel) {
  if (length(truelabel) != length(predlabel))
    stop("truelabel and predlabel must have the same length")
  total = length(truelabel)
  x_ids = unique(truelabel)
  y_ids = unique(predlabel)

  MI = 0.0
  for (idx in x_ids) {
    for (idy in y_ids) {
      idxOccur = which(truelabel == idx)
      idyOccur = which(predlabel == idy)
      idxyOccur = intersect(idxOccur, idyOccur)
      if (length(idxyOccur) > 0) {
        MI = MI + (length(idxyOccur) / total) * log2((length(idxyOccur) * total) /
                                                       (length(idxOccur) * length(idyOccur)))
      }
    }
  }
  Hx = 0

  for (idx in x_ids) {
    idxOccurCount = length(which(truelabel == idx))
    
    Hx = Hx - (idxOccurCount / total) * log2(idxOccurCount / total)
  }
  Hy = 0

  for (idy in y_ids) {
    idyOccurCount = length(which(predlabel == idy))
    Hy = Hy - (idyOccurCount / total) * log2(idyOccurCount / total)
  }
  nmi = 2 * MI / (Hx + Hy)
  tab = table(truelabel, predlabel)
  conv_df = as.data.frame.matrix(tab)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  n2 <- choose(n, 2)
  nis2 <- sum(choose(ni[ni > 1], 2))
  njs2 <- sum(choose(nj[nj > 1], 2))
  ri = 1 + (sum(tab ^ 2) - (sum(ni ^ 2) + sum(nj ^ 2)) / 2) / n2
  ari = c(sum(choose(tab[tab > 1], 2)) - (nis2 * njs2) / n2) / ((nis2 + njs2) /
                                                                  2 - (nis2 * njs2) / n2)
  out = c(nmi,ari)
  names(out)=c("NMI","ARI")
  return(out)
}

k_n <- length(unique(data_label))
NMI <- c()
ARI <- c()
for(j in c(1:20)){
  Kmeans_cluster_result <- kmeans(t(X_new),k_n)$cluster
  NMI <-
    c(NMI, round(evaluation(data_label, Kmeans_cluster_result)[[1]], 4))
  ARI <-
    c(ARI, round(evaluation(data_label, Kmeans_cluster_result)[[2]], 4))
}
Kmeans_NMI <- round(mean(NMI), 3)
Kmeans_ARI <- round(mean(ARI), 3)

Kmeans_ARI
Kmeans_NMI



