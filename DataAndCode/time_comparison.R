library(fda)
library(MFPCA)
library(tidygraph)
library(dplyr)
library(ggplot2)
library(igraph)
library(hawkes)
library(tidyr)
library(poweRlaw)
library(ggforce)
library(aricode)
library(Rtsne)
library(patchwork)
source("DynnetEDA\\functional_ASE.R")
source("DynnetEDA\\downstream_tasks.R")

# independent SBM, hyperparameter = (10, 4)
SBM_generation<-function(n,C,K,B){
  adjacency_matrix=matrix(0,nrow=n,ncol=n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      p = B[C[i],C[j]]
      adjacency_matrix[i,j]=sample(c(0,1),size=1,prob=c(1-p,p))
    }
    adjacency_matrix[i,i]=0.5
  }
  adjacency_matrix+t(adjacency_matrix)
}

N_list = 200
T_list = c(30, 60, 90)

result_df = c()
result_list = vector("list", 9)

for(i in 1:3){
  set.seed(10)
  n_nodes = N_list
  T = T_list[i]

  K = 2
  B = matrix(c(0.3,0.1,0.1,0.3),nrow=2)
  C = c(rep(1,times=n_nodes/2),rep(2,times=n_nodes/2))

  dynamic_adjacency_list = vector(mode="list",10)
  result_list[[i]] = vector("list", 10)
  for(j in 1:10){
    dynamic_adjacency = c()
    dynamic_network_adjacency = array(0, dim = c(n_nodes, n_nodes,T))
    for(t in 1:T){
      adj = SBM_generation(n_nodes,C,K,B)
      dynamic_adjacency = c(dynamic_adjacency,list(list(adj)))
      dynamic_network_adjacency[,,t] = adj
    }
    dynamic_adjacency_list[[j]] = dynamic_network_adjacency
  }

  delta_t1_record = 0
  delta_t2_record = 0
  delta_t3_record = 0
  delta_t1_times = 0
  delta_t2_times = 0
  delta_t3_times = 0
  for(j in 1:10){
    ## fASE method with stepsize = 1e-2
    t1 = Sys.time()
    result_SBM_originalmethod1 = tryCatch({
      fASE(dynamic_adjacency_list[[j]], 2, 10, batch_size = NULL, timestamp_vec = 1:T, scalable_method = FALSE, update_method = "GD", step_size = 1e-2, iteration_step = 500, epsilon=1e-6)},
      error = function(e){
        NA
      })
    t2 = Sys.time()
    delta_t1 = as.numeric(t2-t1, units = "mins")
    ## fASE method with stepsize = 1e-3
    t3 = Sys.time()
    result_SBM_originalmethod2 = tryCatch({
      fASE(dynamic_adjacency_list[[j]], 2, 10, batch_size = NULL, timestamp_vec = 1:T, scalable_method = FALSE, update_method = "GD", step_size = 1e-3, iteration_step = 500, epsilon=1e-6)},
      error = function(e){
        NA
      })
    t4 = Sys.time()
    delta_t2 = as.numeric(t4-t3, units = "mins")
    ## accelerated fASE method
    t5 = Sys.time()
    result_SBM_acceleratedmethod = fASE(dynamic_adjacency_list[[j]], 2, 10, batch_size = NULL, timestamp_vec = 1:T, scalable_method = TRUE, scalable_dim = 20, scalable_power = 6, iteration_step = 500, epsilon=1e-6)
    t6 = Sys.time()
    delta_t3 = as.numeric(t6-t5, units = "mins")

    if(any(is.na(result_SBM_originalmethod1))){
      delta_t1 = Inf
    }else if(result_SBM_originalmethod1[[5]]==500){
      delta_t1 = Inf
    }

    if(any(is.na(result_SBM_originalmethod2))){
      delta_t2 = Inf
    }else if(result_SBM_originalmethod2[[5]]==500){
      delta_t2 = Inf
    }


    ## record computational time data:
    result_list[[i]][[j]] = vector("list", 3)
    result_list[[i]][[j]][[1]] = result_SBM_originalmethod1
    result_list[[i]][[j]][[2]] = result_SBM_originalmethod2
    result_list[[i]][[j]][[3]] = result_SBM_acceleratedmethod

    if(!is.infinite(delta_t1)){
      delta_t1_record = delta_t1_record + delta_t1
    }else{
      delta_t1_times = delta_t1_times + 1
    }

    if(!is.infinite(delta_t2)){
      delta_t2_record = delta_t2_record + delta_t2
    }else{
      delta_t2_times = delta_t2_times + 1
    }

    if(!is.infinite(delta_t3)){
      delta_t3_record = delta_t3_record + delta_t3
    }else{
      delta_t3_times = delta_t3_times + 1
    }
  }

  res = c(delta_t1_record/(10-delta_t1_times), delta_t1_times, delta_t2_record/(10-delta_t2_times), delta_t2_times, delta_t3_record/(10-delta_t3_times), delta_t3_times, n_nodes, T)
  res[is.nan(res)] = Inf
  result_df = rbind(result_df, t(res))

  save(file="time_compare.RData",result_df, result_list)
}

