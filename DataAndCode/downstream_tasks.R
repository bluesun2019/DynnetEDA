library(igraph)
library(fda)
library(fda.usc)
library(MASS)
library(ggplot2)
library(reshape2)
library(patchwork)
library(lubridate)
library(factoextra)
library(fpc)
library(dplyr)
library(tidyr)
library(readr)
library(MFPCA)
library(roahd)
library(ecp)
library(cluster)
library(ClusterR)
library(kselection)
library(lowmemtkmeans)
library(doSNOW)
library(robustbase)
library(Rtsne)

fPCA <- function(embedding_fASE, timestamp_vec, PCA_dim, active_points = NULL){
  m_layers = dim(embedding_fASE)[1]
  embedding_dim = dim(embedding_fASE)[3]
  if(!is.null(active_points)){
    active_number = sapply(1:dim(embedding_fASE)[2], function(x){length(active_points[[x]])})
    active_number[active_number==0] = 1
    embedding_fASE = embedding_fASE/sqrt(aperm(array(active_number, dim = c(dim(embedding_fASE)[2],dim(embedding_fASE)[1],dim(embedding_fASE)[3])),c(2,1,3)))
  }
  
  # the third dim of embedding_fASE should be consistent with the length of timestamp_vec.
  
  if(length(timestamp_vec)!=m_layers){
    cat("Parameters do not fit with each other!\n")
    break
  }else{
    ZZ = vector(mode = "list", length = embedding_dim)
    uniExpansions = vector(mode = "list", length = embedding_dim)
    for(i in 1:embedding_dim){
      if(is.matrix(embedding_fASE[,,i])){
        Zi = funData::funData(argvals = timestamp_vec, X = t(embedding_fASE[,,i]))
      }else{
        Zi = funData::funData(argvals = timestamp_vec, X = as.matrix(embedding_fASE[,,i]))
      }
      
      ZZ[[i]] = Zi
      uniExpansions[[i]] = list(type = "fda")
    }
    
    ZZ = funData::multiFunData(ZZ)
    # pca_fASE = MFPCA::MFPCA(ZZ, PCA_dim, uniExpansions, weights = cumprod(rep(0.8,times=length(ZZ)))/sum(cumprod(rep(0.8,times=length(ZZ)))))
    pca_fASE = MFPCA::MFPCA(ZZ, PCA_dim, uniExpansions)
    
    list(pca_fASE$scores, pca_fASE$functions, summary(pca_fASE)["Proportion of variance explained",])
  }
}
fPCA_angle <- function(embedding_fASE, timestamp_vec, PCA_dim, active_points = NULL){
  sintegral_use <- function(fx, x){
    Bolstad2::sintegral(x, fx)$int
  }
  subset <- function(i, ls, x){
    if(length(ls[[i]])!=0){
      h = matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
      h[which((1:length(timestamp_vec))%in%ls[[i]]),] = x[which((1:length(timestamp_vec))%in%ls[[i]]),i,]
    }else{
      h = matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
    }
    h
  }
  
  embedding_fASE2 = embedding_fASE
  if(!is.null(active_points)){
    active_points = sapply(1:length(active_points), function(i){pmatch(active_points[[i]],timestamp_vec)})
    embedding_fASE2 = sapply(1:dim(embedding_fASE)[2],subset,active_points,embedding_fASE2,simplify = "array")
    embedding_fASE2 = aperm(embedding_fASE2, c(1,3,2))
  }
  embedding_fASE3 = apply(embedding_fASE2, c(1,2), crossprod)
  embedding_fASE3[embedding_fASE3==0] = 1
  
  embedding_fASE = embedding_fASE2/sqrt(array(embedding_fASE3, dim = c(dim(embedding_fASE)[1],dim(embedding_fASE)[2],dim(embedding_fASE)[3])))
  if(!is.null(active_points)){
    active_number = sapply(1:dim(embedding_fASE)[2], function(x){length(active_points[[x]])})
    active_number[active_number==0] = 1
    embedding_fASE = embedding_fASE/sqrt(aperm(array(active_number, dim = c(dim(embedding_fASE)[2],dim(embedding_fASE)[1],dim(embedding_fASE)[3])),c(2,1,3)))
  }
  
  m_layers = dim(embedding_fASE)[1]
  embedding_dim = dim(embedding_fASE)[3]
  # the third dim of embedding_fASE should be consistent with the length of timestamp_vec.
  
  if(length(timestamp_vec)!=m_layers){
    cat("Parameters do not fit with each other!\n")
    break
  }else{
    ZZ = vector(mode = "list", length = embedding_dim)
    uniExpansions = vector(mode = "list", length = embedding_dim)
    for(i in 1:embedding_dim){
      if(is.matrix(embedding_fASE[,,i])){
        Zi = funData::funData(argvals = timestamp_vec, X = t(embedding_fASE[,,i]))
      }else{
        Zi = funData::funData(argvals = timestamp_vec, X = as.matrix(embedding_fASE[,,i]))
      }
      
      ZZ[[i]] = Zi
      uniExpansions[[i]] = list(type = "fda")
    }
    
    ZZ = funData::multiFunData(ZZ)
    pca_fASE = MFPCA::MFPCA(ZZ, PCA_dim, uniExpansions)
    
    list(pca_fASE$scores, pca_fASE$functions, summary(pca_fASE)["Proportion of variance explained",])
  }
}

active_calculation <- function(node, net, timestamp_vec = NULL){
  if(is.null(timestamp_vec)){
    which(apply(net[node,,],2,sum)!=0)
  }else{
    timestamp_vec[which(apply(net[node,,],2,sum)!=0)]
  }
}

cosine_similarity <- function(embedding_fASE, timevec){
  sintegral_use <- function(fx, x){
    Bolstad2::sintegral(x, fx)$int
  }
  cross_mul <- function(i, arr1, arr2){
    res = arr1[i,,]%*%t(arr2[i,,])
  }
  similarity_int = mapply(cross_mul, i=1:dim(embedding_fASE)[1], MoreArgs = list(arr1=embedding_fASE,arr2=embedding_fASE))
  similarity_mat = apply(similarity_int, 1, sintegral_use, timevec)
  similarity_mat = matrix(similarity_mat, nrow=dim(embedding_fASE)[2])
  similarity_mat = similarity_mat/sqrt(diag(similarity_mat)%*%t(diag(similarity_mat)))
  similarity_mat
}
cosine_similarity_partial <- function(embedding_fASE, nodes_calc, timevec_list = NULL, timevec = NULL){
  # timevec_list: how to divide the timestamp intervals
  if(is.null(timevec_list)){
    if(!is.null(timevec))  timevec_list = as.list(timevec)
    else  timevec_list = as.list(1:dim(embedding_fASE)[1])
  }
  sintegral_use <- function(fx, x){
    if(length(x)>1){
      Bolstad2::sintegral(x, fx)$int
    }else{
      fx
    }
  }
  cross_mul <- function(i, arr1, arr2){
    res = arr1[i,,]%*%t(arr2[i,,])
  }
  similarity_mat = list()
  for(s in 1:length(timevec_list)){
    if(length(timevec_list[[s]])!=1){
      arr1 = array(embedding_fASE[timevec_list[[s]],nodes_calc,], dim = c(length(timevec_list[[s]]),length(nodes_calc),dim(embedding_fASE)[3]))
      arr2  =embedding_fASE[timevec_list[[s]],,]
    }else{
      arr1 = aperm(array(embedding_fASE[timevec_list[[s]],nodes_calc,], dim = c(length(nodes_calc),dim(embedding_fASE)[3],length(timevec_list[[s]]))),c(3,1,2))
      arr2 = aperm(array(embedding_fASE[timevec_list[[s]],,], dim = c(dim(embedding_fASE)[2],dim(embedding_fASE)[3],length(timevec_list[[s]]))),c(3,1,2))
    }
    similarity_int = mapply(cross_mul, i=1:dim(arr1)[1], MoreArgs = list(arr1=arr1, arr2 = arr2))
    similarity_int = array(similarity_int, dim = c(dim(arr1)[1], dim(arr1)[2], dim(arr2)[2]))
    norm_int = apply(arr2, c(1,2), function(x){sum(x^2)})
    similarity_mat[[s]] = apply(similarity_int, 1, sintegral_use, timevec_list[[s]])
    norm_vec = apply(norm_int, 1, sintegral_use, timevec_list[[s]])
    similarity_mat[[s]] = matrix(similarity_mat[[s]], nrow=length(nodes_calc))
    similarity_mat[[s]] = similarity_mat[[s]]/sqrt(norm_vec[nodes_calc]%*%t(norm_vec))
  }
  
  similarity_mat
}
minimal_clusters <- function(similarity_mat, nodes_calc, cluster, timevec_list = NULL,  timevec = NULL, active_points = NULL, active_cluster = NULL){
  # timevec_list: how to divide the timestamp intervals
  # active_cluster: a vector given the active time for each cluster (names are cluster names)
  
  cluster = as.numeric(cluster)
  cluster_list = vector("list",length(nodes_calc))
  cluster_frame = matrix(0, nrow = length(nodes_calc), ncol = length(similarity_mat))
  if(is.null(timevec_list)){
    if(!is.null(timevec))  timevec_list = as.list(timevec)
    else  timevec_list = as.list(1:length(similarity_mat))
  }
  if(is.null(active_cluster)){
    if(!is.null(timevec)){
      active_cluster = rep(timevec[1], times = length(max(cluster)))
      names(active_cluster) = 1:length(max(cluster))
    } 
    else{
      active_cluster = rep(1, times = length(max(cluster)))
      names(active_cluster) = 1:length(max(cluster))
    } 
  }
  
  for(s in 1:length(timevec_list)){
    if(length(which(active_cluster<=timevec_list[[s]][1]))==0){
      cluster_frame[,s] = NA
    }else{
      cluster_frame[,s] = apply(1-similarity_mat[[s]],1,function(x,cluster){names(active_cluster[active_cluster<=timevec_list[[s]][1]])[
        which.min(tapply(x,cluster, median)[active_cluster<=timevec_list[[s]][1]])]},cluster=cluster)
    }
  }
  rownames(cluster_frame) = nodes_calc
  colnames(cluster_frame) = 1:length(timevec_list)
  
  for(s in 1:length(nodes_calc)){
    positions <- unlist(lapply(timevec_list, function(x) {
      any(x %in% active_points[[nodes_calc[s]]])
    }))
    name_clus = names(cluster_frame[s,(1:length(timevec_list))[positions]])
    cluster_list[[s]][[1]] = as.numeric(cluster_frame[s,(1:length(timevec_list))[positions]])
    names(cluster_list[[s]][[1]]) = name_clus
    cluster_list[[s]][[2]] = timevec_list[positions]
  }
  
  cluster_list
}
metric_depth <- function(dist_mat){
  n_nodes = dim(dist_mat)[1]
  depth_record = rep(0, times = n_nodes)
  
  for (a in 1:n_nodes){
    count_mat = pmax(matrix(dist_mat[,a],nrow = dim(dist_mat)[1], ncol = dim(dist_mat)[2]), matrix(dist_mat[a,],nrow = dim(dist_mat)[2], ncol = dim(dist_mat)[1], byrow = TRUE))
    count = sum(dist_mat>count_mat) - sum(diag(dist_mat>count_mat))
    depth_record[a] = count/(choose(n_nodes ,2)*2)
  }
  depth_record
}

spherical_transform <- function(X){
  m = length(X)
  r = sqrt(sum(X^2))
  phi = rep(0, times = m)
  phi[1] = r
  temp = 1
  for(i in rev(1:(m-1))){
    phi[i+1] = atan(X[i+1]/(X[i]*temp))
    temp = temp * cos(phi[i+1])
  }
  phi
}
radius_transform <- function(X){
  Y = apply(X, c(1,2), function(x){sqrt(sum(x^2))})
  Y[Y==0] = 1
  X/array(Y, dim = dim(X))
}

cluster_kmeans <- function(X, k){
  set.seed(1)
  calculate_wcss <- function(X, cluster_centers, cluster_assignments) {
    wcss <- 0
    for (i in 1:nrow(X)) {
      # 获取数据点i所属的簇中心和簇分配
      cluster_center <- cluster_centers[cluster_assignments[i], ]
      # 计算欧氏距离的平方
      distance_squared <- sum((X[i, ] - cluster_center)^2)
      wcss <- wcss + distance_squared
    }
    return(wcss)
  }
  Kmeans_result = KMeans_rcpp(X, k)
  clusters = Kmeans_result$clusters
  withinss = calculate_wcss(X, Kmeans_result$centroids, clusters)
  
  list(cluster = clusters, withinss = withinss)
}
# cluster_tkmeans <- function(X, k){
#   calculate_wcss <- function(X, cluster_centers, cluster_assignments) {
#     wcss <- 0
#     for (i in 1:nrow(X)) {
#       # 获取数据点i所属的簇中心和簇分配
#       cluster_center <- cluster_centers[cluster_assignments[i], ]
#       # 计算欧氏距离的平方
#       distance_squared <- sum((X[i, ] - cluster_center)^2)
#       wcss <- wcss + distance_squared
#     }
#     return(wcss)
#   }
#   tkmeans_center = lowmemtkmeans::tkmeans(X, k, alpha = 0.1, nstart = 1000, iter = 300)
#   clusters = c(lowmemtkmeans::nearest_cluster(X, tkmeans_center),recursive=TRUE)
#   withinss = calculate_wcss(X, tkmeans_center, clusters)
#   
#   list(cluster = clusters, withinss = withinss)
# }

cluster_tkmeans <- function(X, k){
  set.seed(1)
  calculate_wcss <- function(X, cluster_centers, cluster_assignments) {
    wcss <- 0
    for (i in 1:nrow(X)) {
      # 获取数据点i所属的簇中心和簇分配
      cluster_center <- cluster_centers[cluster_assignments[i], ]
      # 计算欧氏距离的平方
      distance_squared <- sum((X[i, ] - cluster_center)^2)
      wcss <- wcss + distance_squared
    }
    return(wcss)
  }
  tryCatch({
    if(k==1){
      withinss_record = c()
      clusters_record = matrix(nrow=50,ncol=dim(X)[1])
      for(iter in 1:50){
        clusters = rep(1, times = dim(X)[1])
        tkmeans_center = matrix(apply(X, 2, mean),nrow=1)
        withinss = calculate_wcss(X, tkmeans_center, clusters)
        withinss_record = c(withinss_record, withinss)
        clusters_record[iter,] = clusters
      }
      clusters = clusters_record[which.min(withinss_record),]
      withinss = withinss_record[which.min(withinss_record)]
    }else{
      withinss_record = c()
      clusters_record = matrix(nrow=50,ncol=dim(X)[1])
      for(iter in 1:50){
        init_points <- flexclust::kcca(X, k, flexclust::kccaFamily("kmeans"), control=list(initcent="kmeanspp"))@centers
        tkmeans_center = trimcluster::trimkmeans(X, k, trim = 0.1, runs = 1, points = init_points)$means
        clusters = c(lowmemtkmeans::nearest_cluster(X, tkmeans_center),recursive=TRUE)
        withinss = calculate_wcss(X, tkmeans_center, clusters)
        withinss_record = c(withinss_record, withinss)
        clusters_record[iter,] = clusters
      }
      clusters = clusters_record[which.min(withinss_record),]
      withinss = withinss_record[which.min(withinss_record)]
    }
  }, error = function(e) {
    k = k+1
    if(k==1){
      withinss_record = c()
      clusters_record = matrix(nrow=50,ncol=dim(X)[1])
      for(iter in 1:50){
        clusters = rep(1, times = dim(X)[1])
        tkmeans_center = matrix(apply(X, 2, mean),nrow=1)
        withinss = calculate_wcss(X, tkmeans_center, clusters)
        withinss_record = c(withinss_record, withinss)
        clusters_record[iter,] = clusters
      }
      clusters = clusters_record[which.min(withinss_record),]
      withinss = withinss_record[which.min(withinss_record)]
    }else{
      withinss_record = c()
      clusters_record = matrix(nrow=50,ncol=dim(X)[1])
      for(iter in 1:50){
        init_points <- flexclust::kcca(X, k, flexclust::kccaFamily("kmeans"), control=list(initcent="kmeanspp"))@centers
        tkmeans_center = trimcluster::trimkmeans(X, k, trim = 0.1, runs = 1, points = init_points)$means
        clusters = c(lowmemtkmeans::nearest_cluster(X, tkmeans_center),recursive=TRUE)
        withinss = calculate_wcss(X, tkmeans_center, clusters)
        withinss_record = c(withinss_record, withinss)
        clusters_record[iter,] = clusters
      }
      clusters = clusters_record[which.min(withinss_record),]
      withinss = withinss_record[which.min(withinss_record)]
    }
  })
  
  list(cluster = clusters, withinss = withinss)
}
cluster_no_selection <- function(X, max_k, cluster_method = "tk_means", K_selection = "Silhouette",  weight = NULL, parallel = TRUE){
  # X: FPC for K-means or trimmed K-means; fdata for FADP clustering
  set.seed(1)
  if(cluster_method == "k_means" && K_selection == "Silhouette"){
    if(!is.null(weight)){
      X = t(weight*t(X))
    }
    result = fviz_nbclust(X, cluster_kmeans, method = "silhouette", k.max = max_k)
    optim_k = result$data[which.max(result$data$y),"clusters"]
    optim_k = as.numeric(optim_k)
    optim_result = list(optim_k, cluster_kmeans(X, optim_k), result$data$y)
  }else if(cluster_method == "k_means" && K_selection == "Pham"){
    if(!is.null(weight)){
      X = t(weight*t(X))
    }

    result = kselection::kselection(X, fun_cluster = cluster_kmeans, max_centers = max_k, parallel =  parallel, k_threshold = 1)
    optim_k = result$k
    optim_result = list(optim_k, cluster_tkmeans(X, optim_k), result$f_k)
    if(parallel){
      stopCluster(cl)
    }
  }else if(cluster_method == "tk_means" && K_selection == "Silhouette"){
    if(!is.null(weight)){
      X = t(weight*t(X))
    }
    result = fviz_nbclust(X, cluster_tkmeans, method = "silhouette", k.max = max_k)
    optim_k = result$data[which.max(result$data$y),"clusters"]
    optim_k = as.numeric(optim_k)
    optim_result = list(optim_k, cluster_tkmeans(X, optim_k), result$data$y)
  }else if(cluster_method == "tk_means" && K_selection == "Pham"){
    if(!is.null(weight)){
      X = t(weight*t(X))
    }
    if(parallel){
      cl <- makeCluster(25)
      registerDoSNOW(cl)
    }
    result = kselection::kselection(X, fun_cluster = cluster_tkmeans, max_centers = max_k, parallel = parallel, k_threshold = 1)
    optim_k = result$k
    optim_result = list(optim_k, cluster_tkmeans(X, optim_k), result$f_k)
    if(parallel){
      stopCluster(cl)
    }
  }
  
  optim_result
}

adjacency_visualize <- function(adjacency_matrix, label=NULL, permuted = FALSE, guide_break = NULL, n = NULL){
  reperm<-function(accumulated_matrix,index,s,k){
    if(s>1&&k>1){
      if(s==k){
        acc1=accumulated_matrix[index==s,index==s]
        acc2=reperm(accumulated_matrix,index,s,s-1)
        acc3=reperm(accumulated_matrix,index,s-1,s)
        temp=cbind(acc1,acc2)
        temp2=cbind(acc3,reperm(accumulated_matrix,index,s-1,s-1))
        accumulated_matrix_new=rbind(temp,temp2) 
      }else{
        if(s<k){
          acc=reperm(accumulated_matrix,index,s-1,k)
          if(length(which(index==k))==1){
            accumulated_matrix_new=c(accumulated_matrix[index==s,index==k],acc)
          }else{
            accumulated_matrix_new=rbind(accumulated_matrix[index==s,index==k],acc)
          }
        }else{
          acc=reperm(accumulated_matrix,index,s,k-1)
          if(length(which(index==s))==1){
            accumulated_matrix_new=cbind(t(accumulated_matrix[index==s,index==k]),acc)
          }else{
            accumulated_matrix_new=cbind(accumulated_matrix[index==s,index==k],acc)
          }
        }
      }
    }else{
      acc1=accumulated_matrix[index==s,index==k]
      if(is.null(dim(acc1))&&length(which(index==s))==1){
        acc1=t(acc1)
      }
      accumulated_matrix_new=acc1
    }
    accumulated_matrix_new
  }
  
  
  if(!is.null(label)){
    k=range(label)[2]
    t=rev(table(label))
    t=cumsum(t)
    if(permuted) adjacency_matrix=reperm(adjacency_matrix,label,k,k)
  }
  adjacency_matrix[adjacency_matrix==0]=NA
  colnames(adjacency_matrix)=NULL
  permute=melt(adjacency_matrix)
  node_num = dim(adjacency_matrix)[1]
  permute$Var1=node_num-permute$Var1
  permute = permute[!is.na(permute$value),]
  if(length(setdiff(1:n,c(permute[,1:2],recursive = TRUE)))!=0){
    h = setdiff(1:n,c(permute[,1:2],recursive = TRUE))
    d = rbind(cbind(h,n,NA),cbind(h,1,NA),cbind(n,h,NA),cbind(1,h,NA))
    colnames(d) = c("Var1", "Var2", "value")
    permute = rbind(permute,d)
  }
  
  p1<-ggplot(permute,aes(x=Var1,y=Var2,fill=value))
  if(is.null(guide_break)){
    p2<-p1+geom_tile()+scale_x_discrete(name="Nodes")+scale_y_discrete(name="Nodes")+scale_fill_gradient2(low="midnightblue",mid="cornflowerblue",high="white",na.value = "transparent")
  }else{
    if(max(guide_break)<=3){
      p2<-p1+geom_tile()+scale_x_discrete(name="Nodes")+scale_y_discrete(name="Nodes")+
        scale_fill_gradientn(colours = c("cornflowerblue","royalblue","midnightblue"), breaks = guide_break, labels = guide_break, values = scales::rescale(c(1,max(1+guide_break)/2,max(guide_break))), na.value = "transparent")
    }else{
      p2<-p1+geom_tile()+scale_x_discrete(name="Nodes")+scale_y_discrete(name="Nodes")+
        scale_fill_gradientn(colours = c("cornflowerblue","royalblue","midnightblue"), breaks = guide_break, labels = guide_break, values = scales::rescale(c(1,3,max(guide_break))),na.value = "transparent")
    }
  }
  if(!is.null(label)){
    p3<-p2+geom_segment(x=node_num-t[1]-0.5,xend=node_num-t[1]-0.5,y=0,yend=t[2],color="pink")+geom_segment(x=node_num-t[2]-0.5,xend=node_num-0.5,y=t[1],yend=t[1],color="pink")
    if(k>2){
      for(kk in 2:(k-1)){
        p3<-p3+geom_segment(x=node_num-t[kk]-0.5,xend=node_num-t[kk]-0.5,y=t[kk-1],yend=t[kk+1],color="pink")+geom_segment(x=node_num-t[kk+1]-0.5,xend=node_num-t[kk-1]-0.5,y=t[kk],yend=t[kk],color="pink")
      }
    }
    p3<-p3+theme(legend.position = "left")
  }else{
    p3<-p2
  }
  p3
}

function_sort <- function(X, timestamp_vec,active_points=NULL){
  node_embedding_frame = c()
  node_embedding_frame_name = c()
  
  if(is.matrix(X)){
    n_nodes = dim(X)[2]
    m_layers = dim(X)[1]
    
    if(length(timestamp_vec)!=m_layers){
      cat("Parameters do not fit with each other!\n")
      break
    }else{
      for(i in 1:n_nodes){
        node_embedding_frame = rbind(node_embedding_frame , cbind(X[,i],t = timestamp_vec,i = i))
      }
      colnames(node_embedding_frame)[1]="r"
      node_embedding_frame = as.data.frame(node_embedding_frame)
      
      node_embedding_frame
    }
  }else{
    n_nodes = dim(X)[2]
    m_layers = dim(X)[1]
    embedding_dim = dim(X)[3]
    
    if(length(timestamp_vec)!=m_layers){
      cat("Parameters do not fit with each other!\n")
      break
    }else{
      for(i in 1:n_nodes){
        node_embedding_frame = rbind(node_embedding_frame , cbind(X[,i,],t = timestamp_vec,i = i))
      }
      for(j in 1:embedding_dim){
        node_embedding_frame_name = c(node_embedding_frame_name, paste0("x",j))
      }
      
      colnames(node_embedding_frame)[1:embedding_dim]=node_embedding_frame_name 
      node_embedding_frame = as.data.frame(node_embedding_frame)
      
      node_embedding_frame
    }
  }
}
function_sort_angle <- function(X, timestamp_vec,active_points=NULL){
  node_embedding_frame = c()
  node_embedding_frame_name = c()
  subset <- function(i, ls, x){
    if(length(ls[[i]])!=0){
      h = matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
      h[which((1:length(timestamp_vec))%in%ls[[i]]),] = x[which((1:length(timestamp_vec))%in%ls[[i]]),i,]
    }else{
      h = matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
    }
    h
  }
  subset_mat <- function(i, ls, x){
    if(length(ls[[i]])!=0){
      h = matrix(0, nrow = dim(x)[1], ncol = 1)
      h[which((1:length(timestamp_vec))%in%ls[[i]]),] = x[which((1:length(timestamp_vec))%in%ls[[i]]),i]
    }else{
      h = matrix(0, nrow = dim(x)[1], ncol = 1)
    }
    h
  }
  
  if(!is.null(active_points)){
    active_points = lapply(active_points, function(x){match(x,timestamp_vec)})
  }
  
  if(is.matrix(X)){
    n_nodes = dim(X)[2]
    m_layers = dim(X)[1]
    
    embedding_fASE2 = X
    if(!is.null(active_points)){
      embedding_fASE2 = sapply(1:dim(X)[2],subset_mat,active_points,embedding_fASE2,simplify = "array")
    }
    X = embedding_fASE2
    
    if(!is.null(active_points)){
      active_number = sapply(1:dim(X)[2], function(x){length(active_points[[x]])})
      active_number[active_number==0] = 1
      X = t(t(X)/sqrt(active_number))
    }
    
    
    if(length(timestamp_vec)!=m_layers){
      cat("Parameters do not fit with each other!\n")
      break
    }else{
      for(i in 1:n_nodes){
        node_embedding_frame = rbind(node_embedding_frame , cbind(X[,i],t = timestamp_vec,i = i))
      }
      colnames(node_embedding_frame)[1]="r"
      node_embedding_frame = as.data.frame(node_embedding_frame)
      
      node_embedding_frame
    }
  }else{
    n_nodes = dim(X)[2]
    m_layers = dim(X)[1]
    embedding_dim = dim(X)[3]
    
    embedding_fASE2 = X
    if(!is.null(active_points)){
      embedding_fASE2 = sapply(1:dim(X)[2],subset,active_points,embedding_fASE2,simplify = "array")
      embedding_fASE2 = aperm(embedding_fASE2, c(1,3,2))
    }
    embedding_fASE3 = apply(embedding_fASE2, c(1,2), crossprod)
    embedding_fASE3[embedding_fASE3==0] = 1
    
    X = embedding_fASE2/sqrt(array(embedding_fASE3, dim = c(dim(X)[1],dim(X)[2],dim(X)[3])))
    if(!is.null(active_points)){
      active_number = sapply(1:dim(X)[2], function(x){length(active_points[[x]])})
      active_number[active_number==0] = 1
      X = X/sqrt(aperm(array(active_number, dim = c(dim(X)[2],dim(X)[1],dim(X)[3])),c(2,1,3)))
    }
    
    if(length(timestamp_vec)!=m_layers){
      cat("Parameters do not fit with each other!\n")
      break
    }else{
      for(i in 1:n_nodes){
        node_embedding_frame = rbind(node_embedding_frame , cbind(X[,i,],t = timestamp_vec,i = i))
      }
      for(j in 1:embedding_dim){
        node_embedding_frame_name = c(node_embedding_frame_name, paste0("x",j))
      }
      
      colnames(node_embedding_frame)[1:embedding_dim]=node_embedding_frame_name 
      node_embedding_frame = as.data.frame(node_embedding_frame)
      
      node_embedding_frame
    }
  }
}

get_fASE_angle <- function(X, timestamp_vec,active_points=NULL){
  node_embedding_frame = c()
  node_embedding_frame_name = c()
  subset <- function(i, ls, x){
    if(length(ls[[i]])!=0){
      h = matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
      h[which((1:length(timestamp_vec))%in%ls[[i]]),] = x[which((1:length(timestamp_vec))%in%ls[[i]]),i,]
    }else{
      h = matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
    }
    h
  }
  subset_mat <- function(i, ls, x){
    if(length(ls[[i]])!=0){
      h = matrix(0, nrow = dim(x)[1], ncol = 1)
      h[which((1:length(timestamp_vec))%in%ls[[i]]),] = x[which((1:length(timestamp_vec))%in%ls[[i]]),i]
    }else{
      h = matrix(0, nrow = dim(x)[1], ncol = 1)
    }
    h
  }
  
  if(!is.null(active_points)){
    active_points = lapply(active_points, function(x){match(x,timestamp_vec)})
  }
  
  if(is.matrix(X)){
    X = node_fASE
    n_nodes = dim(X)[2]
    m_layers = dim(X)[1]
    k = max(as.numeric(cluster))
    
    embedding_fASE2 = X
    if(!is.null(active_points)){
      embedding_fASE2 = sapply(1:dim(X)[2],subset_mat,active_points,embedding_fASE2,simplify = "array")
    }
    X = embedding_fASE2
    
    if(!is.null(active_points)){
      active_number = sapply(1:dim(X)[2], function(x){length(active_points[[x]])})
      active_number[active_number==0] = 1
      X = t(t(X)/sqrt(active_number))
    }
  }else{
    X = node_fASE
    n_nodes = dim(X)[2]
    m_layers = dim(X)[1]
    embedding_dim = dim(X)[3]
    k = max(as.numeric(cluster))
    
    embedding_fASE2 = X
    if(!is.null(active_points)){
      embedding_fASE2 = sapply(1:dim(X)[2],subset,active_points,embedding_fASE2,simplify = "array")
      embedding_fASE2 = aperm(embedding_fASE2, c(1,3,2))
    }
    embedding_fASE3 = apply(embedding_fASE2, c(1,2), crossprod)
    embedding_fASE3[embedding_fASE3==0] = 1
    
    X = embedding_fASE2/sqrt(array(embedding_fASE3, dim = c(dim(X)[1],dim(X)[2],dim(X)[3])))
    if(!is.null(active_points)){
      active_number = sapply(1:dim(X)[2], function(x){length(active_points[[x]])})
      active_number[active_number==0] = 1
      X = X/sqrt(aperm(array(active_number, dim = c(dim(X)[2],dim(X)[1],dim(X)[3])),c(2,1,3)))
    }
  }
  X
}
get_centroids_angle <- function(X, timestamp_vec,cluster){
  node_embedding_frame = c()
  node_embedding_frame_name = c()
  subset <- function(i, ls, x){
    if(length(ls[[i]])!=0){
      h = matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
      h[which((1:length(timestamp_vec))%in%ls[[i]]),] = x[which((1:length(timestamp_vec))%in%ls[[i]]),i,]
    }else{
      h = matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
    }
    h
  }
  subset_mat <- function(i, ls, x){
    if(length(ls[[i]])!=0){
      h = matrix(0, nrow = dim(x)[1], ncol = 1)
      h[which((1:length(timestamp_vec))%in%ls[[i]]),] = x[which((1:length(timestamp_vec))%in%ls[[i]]),i]
    }else{
      h = matrix(0, nrow = dim(x)[1], ncol = 1)
    }
    h
  }
  
  if(!is.null(active_points)){
    active_points = lapply(active_points, function(x){match(x,timestamp_vec)})
  }
  
  if(is.matrix(X)){
    n_nodes = dim(X)[2]
    m_layers = dim(X)[1]
    k = max(as.numeric(cluster))
    
    X = get_fASE_angle(X, timestamp_vec)
    
    if(length(timestamp_vec)!=m_layers){
      cat("Parameters do not fit with each other!\n")
      break
    }else{
      centroid_cluster = matrix(0, nrow = dim(X)[1], ncol = k)
      for(i in 1:k){
        centroid_cluster[,i] = apply(X[,cluster==i],1,mean)
      }
    }
  }else{
    n_nodes = dim(X)[2]
    m_layers = dim(X)[1]
    embedding_dim = dim(X)[3]
    k = max(as.numeric(cluster))
    
    X = get_fASE_angle(X, timestamp_vec)
    
    if(length(timestamp_vec)!=m_layers){
      cat("Parameters do not fit with each other!\n")
      break
    }else{
      centroid_cluster = array(0, dim = c(dim(X)[1], k, dim(X)[3]))
      for(i in 1:k){
        centroid_cluster[,i,] = apply(X[,cluster==i,],c(1,3),mean)
      }
    }
  }
  centroid_cluster
}
active_clusters <- function(node_fASE, cluster_fASE, timestamp_vec, cluster, minimal_active_nodes = 10, active_points=NULL){
  find_sublist <- function(lst, element) {
    indices <- which(sapply(lst, function(sublst) element %in% sublst))
    return(indices)
  }
  find_nearest_cluster <- function(cluster_matrix, data_matrix) {
    # calculate distance
    distances <- sapply(1:nrow(data_matrix), function(i) {
      apply(cluster_matrix, 1, function(cluster_center) {
        sqrt(sum((data_matrix[i, ] - cluster_center)^2))
      })
    })
    
    
    # the nearest index
    nearest_clusters <- apply(distances, 2, which.min)
    
    return(nearest_clusters)
  }
  # 若该聚类中已激活节点到聚类质心最小的个数>=5，则认为聚类激活
  
  active_cluster_list = vector(length(timestamp_vec), mode="list")
  for(i in 1:length(timestamp_vec)){
    nodes_active = find_sublist(active_points, timestamp_vec[i])
    nodes_active_positions = node_fASE[i,nodes_active,]
    cluster_positions = cluster_fASE[i,,]
    active_cluster_now = find_nearest_cluster(cluster_positions, nodes_active_positions)
    active_cluster_overall = cluster[nodes_active]
    active_cluster_now = diag(table(active_cluster_now, active_cluster_overall))
    names(active_cluster_now) = 1:length(active_cluster_now)
    active_cluster_now = as.numeric(names(active_cluster_now)[active_cluster_now>=minimal_active_nodes])
    active_cluster_now = sort(unique(active_cluster_now), decreasing = FALSE)
    active_cluster_list[[i]] = active_cluster_now
  }
  active_cluster_list
}

get_list_fd <- function(W, B_splines){
  W_T = aperm(W, c(3,1,2))
  embedding_dim = dim(W_T)[3]
  
  Z_list = vector(mode="list", embedding_dim)
  for(s in 1:embedding_dim){
    Z_list[[s]] = fd(coef = W_T[,,s], basisobj = B_splines)
  }
  Z_list
}
get_fData_fd <- function(embedding_fASE, timestamp_vec){
  m_layers = dim(embedding_fASE)[1]
  embedding_dim = dim(embedding_fASE)[3]
  
  fData_list = vector(mode="list", embedding_dim)
  for(i in 1:embedding_dim){
    fData_list[[i]] = t(embedding_fASE[,,i])
  }
  mfData_object = roahd::mfData(timestamp_vec, fData_list)
  
  mfData_object
}

# membership_visualize <- function(cluster_list, node_label, active_points_list, timestamp_vec, cluster_no, color_set = NULL, timemap = NULL, facet_control = NULL){
#   # timemap is a list including:
#   # mapping vector (POSIXct vector)
#   # unit
#   # date_labels
#   # limits
#   if(length(node_label)!=length(cluster_list)){
#     cat("The length of cluster_list and node_label does not match each other!")
#     break
#   }else{
#     basic_colors = c("#ef7a82","#C82423","#F5502A","#F89100","#FFDA9F","forestgreen","#A5CE9D","#82B0D2","#126085","#c9a2c6","#992f87")
#     if(is.null(color_set)){
#       color_set = basic_colors
#     }
#     
#     
#     cluster_frame = c()
#     for(i in 1:length(cluster_list)){
#       cluster_thisnode = rep(NA, times = length(timestamp_vec))
#       cluster_thisnode[active_points_list[[i]]] = c(cluster_list[[i]][[1]], recursive=TRUE)
#       
#       cluster_frame = rbind(cluster_frame, data.frame(clusters = cluster_thisnode, nodes = node_label[i], t = timestamp_vec))
#     }
#     cluster_frame$clusters = factor(cluster_frame$clusters, levels = cluster_no)
#     
#     
#     if(is.null(timemap)){
#       if(is.null(facet_control)){
#         cluster_changeplot = ggplot(cluster_frame, aes(color = factor(nodes), x = t, y = clusters, group = factor(nodes))) +
#           geom_line() + facet_wrap(~factor(nodes)) +
#           labs(y = "cluster") + # Y axis label
#           scale_color_manual(values = color_set)+
#           theme_bw() +  # 使用简洁主题
#           theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black')) 
#       }else{
#         cluster_changeplot = ggplot(cluster_frame, aes(color = factor(nodes), x = t, y = clusters, group = factor(nodes))) +
#           geom_line() + facet_wrap(~factor(nodes), nrow = facet_control$nrow, ncol = facet_control$ncol) +
#           labs(y = "cluster") + # Y axis label
#           scale_color_manual(values = color_set)+
#           theme_bw() +  # 使用简洁主题
#           theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black')) 
#       }
#     }else{
#       if(is.null(facet_control)){
#         cluster_frame$time = timemap$mappingvec[cluster_frame$t]
#         cluster_changeplot = ggplot(cluster_frame, aes(color = factor(nodes), x = time, y = clusters, group = factor(nodes))) +
#           geom_line() + facet_wrap(~factor(nodes)) +
#           scale_x_datetime(name = timemap$unit, date_labels = timemap$date_labels, limits = timemap$limits) +  # 格式化X轴时间显示
#           labs(y = "cluster") + # Y axis label
#           scale_color_manual(values = color_set)+
#           theme_bw() +  # 使用简洁主题
#           theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black')) 
#       }else{
#         cluster_frame$time = timemap$mappingvec[cluster_frame$t]
#         cluster_changeplot = ggplot(cluster_frame, aes(color = factor(nodes), x = time, y = clusters, group = factor(nodes))) +
#           geom_line() + facet_wrap(~factor(nodes), nrow = facet_control$nrow, ncol = facet_control$ncol) +
#           scale_x_datetime(name = timemap$unit, date_labels = timemap$date_labels, limits = timemap$limits) +  # 格式化X轴时间显示
#           labs(y = "cluster") + # Y axis label
#           scale_color_manual(values = color_set)+
#           theme_bw() +  # 使用简洁主题
#           theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black')) 
#       }
#     }
#     cluster_changeplot 
#   }
# }
membership_visualize <- function(cluster_list, node_label, active_points_list, timestamp_vec, cluster_no, color_set = NULL, ltype_set = NULL, timemap = NULL, facet_control = NULL){
  # timemap is a list including:
  # mapping vector (POSIXct vector)
  # unit
  # date_labels
  # limits
  if(length(node_label)!=length(cluster_list)){
    cat("The length of cluster_list and node_label does not match each other!")
    break
  }else{
    basic_colors = c("#ef7a82","#C82423","#F5502A","#F89100","#FFDA9F","forestgreen","#A5CE9D","#82B0D2","#126085","#c9a2c6","#992f87")
    if(is.null(color_set)){
      color_set = basic_colors
    }
    
    basic_ltype = c("solid", "dashed", "dotted", "longdash")
    if(is.null(ltype_set)){
      ltype_set = basic_ltype
    }

    cluster_frame = c()
    for(i in 1:length(cluster_list)){
      cluster_thisnode = rep(NA, times = length(timestamp_vec))
      cluster_thisnode[active_points_list[[i]]] = c(cluster_list[[i]][[1]], recursive=TRUE)

      cluster_frame = rbind(cluster_frame, data.frame(clusters = cluster_thisnode, nodes = node_label[i], t = timestamp_vec))
    }
    cluster_frame = cluster_frame[!is.na(cluster_frame$clusters),]
    for(node in unique( cluster_frame$nodes)){
      ltp = ifelse(abs(diff(cluster_frame$clusters[cluster_frame$nodes==node])) < 0.1, "solid", "dashed")
      cluster_frame$linetype[cluster_frame$nodes==node] = c("dashed", ltp)
      cluster_frame$lag_t[cluster_frame$nodes==node] =  lag(cluster_frame$t[cluster_frame$nodes==node])
      cluster_frame$lag_clusters[cluster_frame$nodes==node] =  lag(cluster_frame$clusters[cluster_frame$nodes==node])
    }
    cluster_frame$linetype = factor(cluster_frame$linetype , levels = c("solid", "dashed", "dotted", "longdash"))
    cluster_frame$clusters = factor(cluster_no[cluster_frame$clusters], levels = cluster_no)
    cluster_frame$lag_clusters = factor(cluster_no[cluster_frame$lag_clusters], levels = cluster_no)
    cluster_frame$nodes = factor(cluster_frame$nodes)
    cluster_frame2 = cluster_frame
    cluster_frame = cluster_frame[!is.na(cluster_frame$lag_clusters),]
    
    if(is.null(timemap)){
      if(is.null(facet_control)){
        cluster_changeplot = ggplot(cluster_frame, aes(color = factor(nodes), x = t, y = clusters, group = factor(nodes))) +
          geom_segment(aes(x = lag_t, xend = t, y = lag_clusters, yend = clusters, linetype = linetype))  + 
          geom_point(data = cluster_frame2)+
          facet_wrap(~factor(nodes)) + guides(linetype = FALSE) +
          labs(y = "cluster") + # Y axis label
          scale_color_manual(values = color_set)+
          scale_linetype_manual(values = ltype_set)+
          scale_x_continuous(limits = c(min(timestamp_vec), max(timestamp_vec)))+
          theme_bw() +  # 使用简洁主题
          theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
      }else{
        cluster_changeplot = ggplot(cluster_frame, aes(color = factor(nodes), x = t, y = clusters, group = factor(nodes))) +
          geom_segment(aes(x = lag_t, xend = t, y = lag_clusters, yend = clusters, linetype = linetype))  + 
          geom_point(data = cluster_frame2)+
          guides(linetype = FALSE) +
          facet_wrap(~factor(nodes), nrow = facet_control$nrow, ncol = facet_control$ncol) +
          labs(y = "cluster") + # Y axis label
          scale_color_manual(values = color_set)+
          scale_linetype_manual(values = ltype_set)+
          scale_x_continuous(limits = c(min(timestamp_vec), max(timestamp_vec)))+
          theme_bw() +  # 使用简洁主题
          theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
      }
    }else{
      if(is.null(facet_control)){
        cluster_frame$time = timemap$mappingvec[cluster_frame$t]
        cluster_frame$lag_time = timemap$mappingvec[cluster_frame$lag_t]
        cluster_frame2$time = timemap$mappingvec[cluster_frame2$t]
        cluster_changeplot = ggplot(cluster_frame, aes(color = factor(nodes), x = time, y = clusters)) +
          geom_segment(aes(x = lag_time, xend = time, y = lag_clusters, yend = clusters, linetype = linetype))  + 
          geom_point(data = cluster_frame2, aes(x=time, y = clusters))+
          facet_wrap(~factor(nodes)) + guides(linetype = FALSE) +
          scale_x_datetime(name = timemap$unit, date_labels = timemap$date_labels, limits = timemap$limits) +  # 格式化X轴时间显示
          scale_y_discrete(labels = function(x) ifelse(is.na(x), "", x)) +
          labs(y = "cluster") + # Y axis label
          scale_color_manual(values = color_set)+
          scale_linetype_manual(values = ltype_set)+
          theme_bw() +  # 使用简洁主题
          theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
      }else{
        cluster_frame$time = timemap$mappingvec[cluster_frame$t]
        cluster_frame$lag_time = timemap$mappingvec[cluster_frame$lag_t]
        cluster_frame2$time = timemap$mappingvec[cluster_frame2$t]
        cluster_changeplot = ggplot(cluster_frame, aes(color = factor(nodes), x = time, y = clusters, group = factor(nodes))) +
          geom_segment(aes(x = lag_time, xend = time, y = lag_clusters, yend = clusters, linetype = linetype))  + 
          geom_point(data = cluster_frame2)+
          guides(linetype = FALSE) + 
          facet_wrap(~factor(nodes), nrow = facet_control$nrow, ncol = facet_control$ncol) +
          scale_x_datetime(name = timemap$unit, date_labels = timemap$date_labels, limits = timemap$limits) +  # 格式化X轴时间显示
          labs(y = "cluster") + # Y axis label
          scale_color_manual(values = color_set)+
          scale_linetype_manual(values = ltype_set)+
          theme_bw() +  # 使用简洁主题
          theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
      }
    }
    cluster_changeplot
  }
}
