library(fda)
library(MFPCA)
library(latex2exp)
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
library(ggalt)
library(patchwork)
source("DataAndCode\\functional_ASE.R")
source("DataAndCode\\downstream_tasks.R")


# This R script is used to introduce the methods;
# 1.set the SBM + 10 hub nodes independent with time as an example.
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

set.seed(10)

n_nodes = 50
K = 2
T1 = 50
T2 = 100
T = T2
B = matrix(c(0.3,0.05,0.05,0.3),nrow=2)
C = c(rep(1,times=30),rep(2,times=20))
dynamic_adjacency = c()
dynamic_network_adjacency = array(0, dim = c(n_nodes, n_nodes,T2))

## snapshots:
for(T in 1:T1){
  adj = SBM_generation(n_nodes,C,K,B)
  dynamic_adjacency = c(dynamic_adjacency,list(list(adj)))
  dynamic_network_adjacency[,,T] = adj
}

K = 2
B = matrix(c(0.3,0.1,0.1,0.3),nrow=2)
C = c(rep(1,times=20),rep(2,times=30))
for(T in T1:T2){
  adj = SBM_generation(n_nodes,C,K,B)
  dynamic_adjacency = c(dynamic_adjacency,list(list(adj)))
  dynamic_network_adjacency[,,T] = adj
}

CC = c(rep(1,times=20),rep(3,times=10),rep(2,times=20))
adjacency_density = vector(mode = "list", T)
for(t in 1:T){
  adjacency_density[[t]] = adjacency_visualize(dynamic_network_adjacency[,,t], CC, permuted = TRUE, n=50)
}

plot_SBM1_snap = adjacency_density[[1]]+ labs(caption = "(a) T=1.") +
  adjacency_density[[50]]+ labs(caption = "(b) T=50.") +
  adjacency_density[[100]]+ labs(caption = "(c) T=100.") &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'), legend.position = "none",
        plot.margin = unit(c(0.3, 1.5, 0.1, 0.15), "cm"))
ggsave(plot_SBM1_snap ,file="SBM_example_fASE0_snapshots.pdf",width=24,height=7.5)

############ You may skip these costly procedures of embedding.
tuning_sim0 = fASE_tuning(10:20, 1:5, tuning_method = "heuristic", adjacency_tensor = dynamic_network_adjacency, spline_order = 4, batch_size = NULL, kernel_scale = FALSE , timestamp_vec = NULL, scalable_method = TRUE, scalable_dim = 20, scalable_power = 6, iteration_method = "sequential", iteration_step = 2500, epsilon=1e-6)
tuning_parameter0 = get_tuning(tuning_sim0)

result_SBM0 = fASE(dynamic_network_adjacency, tuning_parameter0[[2]], tuning_parameter0[[1]], batch_size = NULL, timestamp_vec = 1:50,scalable_method = TRUE, scalable_dim = 20, scalable_power = 6, iteration_step = 2500, epsilon=1e-6)
embedding_SBM0 = eval.fd(1:T, result_SBM0[[1]])
save(file="SBM_example_fASE0.RData",tuning_parameter0, result_SBM0, embedding_SBM0)
########################

## CD plot and functions:
load("vary_SBM_example_fASE0.RData")
PCA_result = fPCA(embedding_SBM0, 1:T, 40)
(l = which(cumsum(PCA_result[[3]])>0.8))
s = cluster_no_selection(PCA_result[[1]][,1:l[1]], 15,cluster_method="tk_means")
s_perm = s[[2]]$cluster
s_perm[s_perm==2]=4
s_perm[s_perm==3]=2
s_perm[s_perm==4]=3

PCA_score = fPCA(embedding_SBM0, 1:T, 2)
PCA_score_show = PCA_score[[1]]
PCA_score_show = as.data.frame(cbind(PCA_score_show,0))
names(PCA_score_show) = c("x1", "x2", "col")
PCA_score_show$label = factor(CC)
PCA_score_show$cluster = factor(s_perm)
score_plot0 = ggplot(PCA_score_show,aes(x=x1,y=x2,group=cluster,color=cluster))+geom_point(size=2)+
  geom_encircle(aes(group = cluster,fill=cluster),expand=0.02,spread=0.5,s_shape=1,size=3,linetype = 1,alpha=0.2)+
  scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+
  scale_fill_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+guides(fill='none')+
  theme_bw()
node_embedding_frame0 = function_sort(embedding_SBM0, 1:T)
node_embedding_frame0 = cbind(node_embedding_frame0, cluster = factor(rep(s_perm, each = T)))
node_embedding_frame0_star = filter(node_embedding_frame0, cluster==3)
func_plot_best_x1 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("x",1),group="i",color="cluster"))+geom_line(aes(alpha=cluster,size=cluster))+geom_point(size = 1)+scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+scale_size_manual(values= c(0.8,0.8,1))+scale_alpha_manual(values= c(0.8,0.8,1))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()
func_plot_best_x1 = func_plot_best_x1 +geom_point(size = 1)+ geom_line(data = node_embedding_frame0_star, alpha=1,size=1,color="#FC8677")
func_plot_best_x2 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("x",2),group="i",color="cluster"))+geom_line(aes(alpha=cluster,size=cluster))+geom_point(size = 1)+scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+scale_size_manual(values= c(0.8,0.8,1))+scale_alpha_manual(values= c(0.8,0.8,1))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()
func_plot_best_x2 = func_plot_best_x2 +geom_point(size = 1)+ geom_line(data = node_embedding_frame0_star, alpha=1,size=1,color="#FC8677")

layout <- "
AABBCC
AABBCC
"
plot_SBM0 = score_plot0+xlab(TeX("$x_1$"))+ylab(TeX("$x_2$"))+guides(color="none")+labs(caption = "(a) CD (community detection) plot.  ")+
  func_plot_best_x1+ylab(TeX("$Z_1$"))+labs(caption ="(b) first dimension of the \n multivariate functional data.")+
  func_plot_best_x2+ylab(TeX("$Z_2$"))+labs(caption ="(c) second dimension of the \n multivariate functional data.")+
  plot_layout(design = layout, guides = "collect") &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(plot_SBM0,file="SBM_example_fASE0.pdf",width=24,height=7.5)


# 2.set the 0-1 DCSBM + 10 hub nodes independent with time as an example.
theta_w_generation<-function(n,C,K,type,xmin=1,scaling=2.5,lambda=0.5){
  deg_vec=rpldis(n,xmin,scaling,discrete_max=n)
  theta=deg_vec
  kappa = rep(0, times=K)
  for(k in 1:K){
    theta[C==k]=deg_vec[C==k]/sum(deg_vec[C==k])
  }
  kappa=tapply(deg_vec,C,sum)
  m = sum(deg_vec)
  w_random=outer(kappa,kappa)/m
  if(type=="group"){
    w_planted = diag(kappa)
  }else if(type=="core-periphery") {
    w_planted = diag(kappa)
    K_max = which.max(kappa)
    K_min = which.min(kappa)
    w_planted = matrix(0, nrow=K, ncol=K)
    w_planted[K_min,K_min] = kappa[K_min]*0.3
    w_planted[K_max, K_min] = kappa[K_min]*0.7
    w_planted[K_min, K_max] = kappa[K_min]*0.7
    w_planted[K_max, K_max] = kappa[K_max]-kappa[K_min]*0.7
  }
  w = lambda*w_random + (1-lambda)*w_planted
  list(theta=theta,w=w)
}
DCBM_indeptot_generation <- function(n,nn,C,K,theta_w,p,T){
  theta = theta_w$theta
  w = theta_w$w
  adj = array(0, dim = c(n+nn,n+nn,T))
  for(t in 1:T){
    adj_init = matrix(0, nrow=n+nn,ncol=n+nn)
    for(i in 1:n){
      for(j in 1:i){
        adj_init[i,j] = rpois(1,theta[i]*w[C[i],C[j]]*theta[j])
        adj_init[j,i] = adj_init[i,j]
      }
    }

    for(i in (n+1):(n+nn)){
      for(j in 1:n){
        adj_init[i,j] = rpois(1,p*theta[j])
        adj_init[j,i] = adj_init[i,j]
      }
      for(j in (n+1):i){
        adj_init[i,j] = rpois(1,p^2/(mean(diag(w))))
        adj_init[j,i] = adj_init[i,j]
      }
    }

    adj[,,t] = adj_init
  }

  adj
}
set.seed(10)

n_nodes = 110
T = 100
K = 5
C = c(rep(1,times=25),rep(2,times=25),rep(3,times=25),rep(4,times=25))
theta_w = theta_w_generation(100,C,K,type="group",lambda=0.4,xmin=10,scaling=2.5)
p = 15
dynamic_network_adjacency = DCBM_indeptot_generation(100,10,C,K,theta_w,p,T)

C = c(rep(1,times=25),rep(2,times=25),rep(3,times=25),rep(4,times=25),rep(5,times=10))

## snapshots:
adjacency_density = vector(mode = "list", T)
for(t in 1:T){
  adjacency_density[[t]] = adjacency_visualize(dynamic_network_adjacency[,,t], C, permuted = TRUE, guide_break = c(5,10,15), n=110)
}

plot_DCBM_snap = adjacency_density[[1]]+ labs(caption = "(a) T=1.") +
  adjacency_density[[2]]+ labs(caption = "(b) T=2.") + guides(fill=FALSE)+
  adjacency_density[[3]]+ labs(caption = "(c) T=3.") + guides(fill=FALSE)+
  plot_layout(guides = "collect") &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'), legend.position = "none",
        plot.margin = unit(c(0.3, 1.5, 0.1, 0.15), "cm"))
ggsave(plot_DCBM_snap ,file="method_examples/DCBM_example_fASE0_snapshots.pdf",width=24,height=7.5)

############ You may skip these costly procedures of embedding.
tuning_sim0 = fASE_tuning(20:30, 1:5, tuning_method = "heuristic", adjacency_tensor = dynamic_network_adjacency, spline_order = 4, batch_size = NULL, kernel_scale = FALSE , timestamp_vec = NULL, scalable_method = TRUE, scalable_dim = 20, scalable_power = 6, iteration_method = "sequential", iteration_step = 2500, epsilon=1e-6)
tuning_parameter0 = get_tuning(tuning_sim0)

result_DCBM0 = fASE(dynamic_network_adjacency, tuning_parameter0[[2]], tuning_parameter0[[1]], batch_size = NULL, timestamp_vec = 1:T, scalable_dim = 20, scalable_power = 6, iteration_step = 2500, epsilon=1e-6)
embedding_DCBM0 = eval.fd(1:T, result_DCBM0[[1]])
save(file="DCBM_example_fASE0.RData",tuning_parameter0, result_DCBM0)
####################

load("DCBM_example_fASE0.RData")
embedding_DCBM0 = eval.fd(1:T, result_DCBM0[[1]])

## CD plot and functions:
PCA_result = fPCA(embedding_DCBM0, 1:T, 100)
(l = which(cumsum(PCA_result[[3]])>0.8))
s = cluster_no_selection(PCA_result[[1]][,1:max(l[1],2)], 15,cluster_method="tk_means")
s_perm = s[[2]]$cluster

PCA_score = fPCA(embedding_DCBM0, 1:T, 2)
PCA_score_show = PCA_score[[1]]
PCA_score_show = as.data.frame(cbind(PCA_score_show,0))
names(PCA_score_show) = c("x1", "x2", "col")
PCA_score_show$label = factor(C)
PCA_score_show$cluster_uncorrected = factor(s_perm)
score_plot0 = ggplot(PCA_score_show,aes(x=x1,y=x2,group=cluster_uncorrected,color=cluster_uncorrected))+geom_point()+geom_encircle(aes(group = cluster_uncorrected,fill=cluster_uncorrected),expand=0.05,spread=0.03,s_shape=1,size=3,linetype = 1,alpha=0.2)+
  scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+
  scale_fill_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+guides(fill='none')+
  theme_bw()

node_embedding_frame0 = function_sort(embedding_DCBM0, 1:T)
node_embedding_frame0 = cbind(node_embedding_frame0, cluster_uncorrected = factor(rep(s_perm, each = T)))
func_plot_best_x1 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("x",1),group="i",color="cluster_uncorrected"))+geom_line(aes(alpha=cluster_uncorrected,size=cluster_uncorrected))+geom_point(size = 1)+scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+scale_size_manual(values= c(0.5,0.5,0.5,0.5,0.8))+scale_alpha_manual(values= c(0.8,0.8,0.8,0.8,1))+guides(color=guide_legend(override.aes = list(size=3),title="cluster\nuncorrected"),size=guide_legend(title="cluster\nuncorrected"),alpha=guide_legend(title="cluster\nuncorrected"))+theme_bw()
func_plot_best_x2 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("x",2),group="i",color="cluster_uncorrected"))+geom_line(aes(alpha=cluster_uncorrected,size=cluster_uncorrected))+geom_point(size = 1)+scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+scale_size_manual(values= c(0.5,0.5,0.5,0.5,0.8))+scale_alpha_manual(values= c(0.8,0.8,0.8,0.8,1))+guides(color=guide_legend(override.aes = list(size=3),title="cluster\nuncorrected"),size=guide_legend(title="cluster\nuncorrected"),alpha=guide_legend(title="cluster\nuncorrected"))+theme_bw()

## DCCD plot and radius-corrected functions:
PCA_result = fPCA_angle(embedding_DCBM0, 1:T, 100)
(l = which(cumsum(PCA_result[[3]])>0.8))
s = cluster_no_selection(PCA_result[[1]][,1:max(l[1],2)], 15,cluster_method="tk_means")
s_perm = s[[2]]$cluster

PCA_score = fPCA_angle(embedding_DCBM0, 1:T, 2)
PCA_score_angle = as.data.frame(PCA_score[[1]])
names(PCA_score_angle) = c("theta1", "theta2")
PCA_score_angle$label = factor(C)
PCA_score_angle$cluster_corrected = factor(s_perm)
angle_plot0 = ggplot(PCA_score_angle,aes(x=theta1,y=theta2,group=cluster_corrected,color=cluster_corrected))+geom_point()+geom_encircle(aes(group = cluster_corrected,fill=cluster_corrected),expand=0.03,spread=0.5,s_shape=1,size=3,linetype = 1,alpha=0.2)+
  scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+
  scale_fill_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+guides(fill='none')+
  theme_bw()

node_embedding_frame0 = function_sort_angle(embedding_DCBM0, 1:T)
node_embedding_frame0 = cbind(node_embedding_frame0, cluster_corrected = factor(rep(s_perm, each = T)))
names(node_embedding_frame0)[1:4] = paste0("theta",1:4)
func_plot_best_a1 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("theta",1),group="i",color="cluster_corrected"))+geom_line(aes(alpha=cluster_corrected,size=cluster_corrected))+geom_point(size = 1)+scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+scale_size_manual(values= c(0.5,0.5,0.5,0.5,0.8))+scale_alpha_manual(values= c(0.8,0.8,0.8,0.8,1))+guides(color=guide_legend(override.aes = list(size=3),title="cluster\ncorrected"),size=guide_legend(title="cluster\ncorrected"),alpha=guide_legend(title="cluster\ncorrected"))+theme_bw()
func_plot_best_a2 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("theta",2),group="i",color="cluster_corrected"))+geom_line(aes(alpha=cluster_corrected,size=cluster_corrected))+geom_point(size = 1)+scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+scale_size_manual(values= c(0.5,0.5,0.5,0.5,0.8))+scale_alpha_manual(values= c(0.8,0.8,0.8,0.8,1))+guides(color=guide_legend(override.aes = list(size=3),title="cluster\ncorrected"),size=guide_legend(title="cluster\ncorrected"),alpha=guide_legend(title="cluster\ncorrected"))+theme_bw()

## DD plot and tSNE embeddings.
set.seed(10)
sim=cosine_similarity(eval.fd(seq(1,T,length=1000), result_DCBM0[[1]]), seq(1,T,length=1000))
tsne_result <- Rtsne(1-sim, dims = 2, perplexity = 30, is_distance=TRUE)
tsne_embeddings <- tsne_result$Y
tsne_embeddings=as.data.frame(tsne_embeddings)
colnames(tsne_embeddings) = c("x1","x2")
sim=cosine_similarity(eval.fd(1:100, result_DCBM0[[1]]), 1:100)
depth_stats = metric_depth(1-sim)
depth_reindex = order(depth_stats,decreasing =TRUE)
tsne_embeddings[depth_reindex,"depth_order"] = as.numeric(1:dim(tsne_embeddings)[1])
tsne_embeddings$depth = depth_stats
tsne_embeddings$cluster_corrected = factor(s_perm)
tsne_embeddings$degree = apply(dynamic_network_adjacency, 1,sum)
tsne_pos = ggplot(tsne_embeddings,aes(x=x1,y=x2))+
  geom_point(size = 2.5, aes(color=cluster_corrected,shape=cluster_corrected))+scale_shape_manual(values =  c(15, 15, 15, 15, 17))+
  scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+
  guides(color=FALSE,fill=FALSE,shape=FALSE)+theme_bw()
dd_plot = ggplot(tsne_embeddings,aes_string(x="depth",y="degree",color="cluster_corrected"))+geom_point(aes(shape = cluster_corrected),size=3)+scale_y_log10()+annotation_custom(
  grob = grid::linesGrob(y = unit(c(1, 0), "npc")), xmin = -Inf, ymin =Inf, xmax = Inf, ymax = -Inf)+
  geom_mark_rect(data=tsne_embeddings%>%filter(cluster_corrected==5),aes(group = cluster_corrected,fill=cluster_corrected),expand=0.02,size=0.8,linetype = 1,alpha=0.2)+
  scale_shape_manual(values =  c(15, 15, 15, 15, 17))+
  scale_color_manual(values=c("#93b69c","#FFBE7A","#FC8677","#82B0D2","#BEB8DC","#E7DAD2","gray"))+
  scale_fill_manual(values=c("violetred","#E7DAD2","gray"))+
  guides(color=guide_legend(title="cluster\ncorrected"),size=guide_legend(title="cluster\ncorrected"),shape=guide_legend(title="cluster\ncorrected"),fill="none")+
  theme_bw()+theme(text = element_text(size = 28),axis.title = element_text(size = 24),legend.text = element_text(size = 24))

layout <- "
CCAABB
CCAABB
DDEEFF
DDEEFF
"
plot_DCBM0 = func_plot_best_x1+ylab(TeX("$Z_1$"))+labs(caption ="(b) first dimension of the \n multivariate functional data.")+
  func_plot_best_x2+ylab(TeX("$Z_2$"))+labs(caption ="(c) second dimension of the \n multivariate functional data.")+
  score_plot0+xlab(TeX("$x_1$"))+ylab(TeX("$x_2$"))+guides(color="none")+labs(caption = "(a) CD (community detection) plot.\n")+
  angle_plot0+xlab(TeX("$x_1^*$"))+ylab(TeX("$x_2^*$"))+guides(color="none")+labs(caption = "(d) DCCD (degree-corrected \n community detection) plot.")+
  func_plot_best_a1+ylab(TeX("$Z_1^*$"))+labs(caption ="(e) first dimension of the multivariate \n radius-corrected functional data.")+
  func_plot_best_a2+ylab(TeX("$Z_2^*$"))+labs(caption ="(f) second dimension of the multivariate \n radius-corrected functional data.")+
  plot_layout(design = layout, guides = "collect") &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(plot_DCBM0,file="DCBM_example_fASE0.pdf",width=26,height=16)

layout <- "
EAACCBBD
EAACCBBD
"
plot_DCBM00 = dd_plot+labs(caption = "(a) DD(degree-depth) plot for nodes.")+
  tsne_pos+xlab(TeX("$y_1$"))+ylab(TeX("$y_2$"))+labs(caption = "(b) the tsne embedding of nodes using sine values.")+
  plot_spacer()+plot_spacer()+plot_spacer()+
  plot_layout(guides = "collect", design = layout, heights = c(4,4),width = c(4.5,7,7,3,3,7,7,4.5)) &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(plot_DCBM00,file="DCBM_example_fASE0_depth.pdf",width=24,height=7.5)


# 3.set the time-varying 0-1 SBM.
SBM_varynodes <- function(n,C,K,B,t){
  adjacency_matrix=matrix(0,nrow=n,ncol=n)

  if(t<=20){
    active_nodes = unique(c(1:30, 61:90))
  }else if(t>60){
    active_nodes = unique(c(31:60, 91:120))
  }else if(t>20&&t<=60){
    active_nodes = c(16:45,76:105)
  }
  for(i in active_nodes){
    for(j in active_nodes){
      if(i>j){
        p = B[C[i],C[j]]
        adjacency_matrix[i,j]=sample(c(0,1),size=1,prob=c(1-p,p))
      }
    }
    adjacency_matrix[i,i]=0.5
  }
  adjacency_matrix+t(adjacency_matrix)
}

set.seed(10)
n_nodes = 120
T = 80
K = 2
B = matrix(c(0.3,0.1,0.1,0.3),nrow=2)
C = c(rep(1,times=60),rep(2,times=60))
dynamic_adjacency = c()
dynamic_network_adjacency = array(0, dim = c(n_nodes, n_nodes,T))
for(t in 1:T){
  adj = SBM_varynodes(n_nodes,C,K,B,t)
  dynamic_adjacency = c(dynamic_adjacency,list(list(adj)))
  dynamic_network_adjacency[,,t] = adj
}

## snapshot:
adjacency_density = vector(mode = "list", T)
for(t in 1:T){
  adjacency_density[[t]] = adjacency_visualize(dynamic_network_adjacency[,,t],C,n=n_nodes)
}
plot_vary1_snap = adjacency_density[[1]]+ labs(caption = "(a) T=1.") + guides(fill="none")+
  adjacency_density[[40]]+ labs(caption = "(e) T=40.") +  guides(fill="none")+
  adjacency_density[[80]]+ labs(caption = "(f) T=80.") +  guides(fill="none")+
  plot_layout(guides = "collect") &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'), legend.position = "none",
        plot.margin = unit(c(0.3, 1.5, 0.1, 0.15), "cm"))
ggsave(plot_vary1_snap ,file="varynodes_example_fASE0_snapshots.pdf",width=24,height=7.5)

############ You may skip these costly procedures of embedding.
tuning_sim0 = fASE_tuning(5:15, 1:10, tuning_method = "heuristic", adjacency_tensor = dynamic_network_adjacency, spline_order = 4,  batch_size = 20,  kernel_scale = FALSE , timestamp_vec = NULL, scalable_method = TRUE, scalable_dim = 30, scalable_power = 6, iteration_method = "sequential", iteration_step = 2500, epsilon=1e-6, parallel = TRUE)
tuning_parameter0 = get_tuning(tuning_sim0)

result_varySBM0 = fASE(dynamic_network_adjacency, tuning_parameter10[[2]], tuning_parameter10[[1]],  batch_size = 20,  timestamp_vec = 1:T, scalable_dim = 30, scalable_power = 6,  iteration_step = 2500, epsilon=1e-6)
embedding_varySBM0 = eval.fd(1:T, result_varySBM0[[1]])
save(file="vary_SBM_example_fASE0.RData",tuning_parameter0, result_varySBM0, embedding_varySBM0)
#############################

load("vary_SBM_example_fASE0.RData")
active_points = lapply(1:120, active_calculation, dynamic_network_adjacency)
mean_time = sapply(active_points, function(x){mean(x)},simplify=TRUE)
active_duration = sapply(active_points, function(x){length(x)},simplify=TRUE)

## generalized CD plot:
PCA_result = fPCA(embedding_varySBM0,  1:80, 100, active_points)
(l = which(cumsum(PCA_result[[3]])>0.8))
s = cluster_no_selection(PCA_result[[1]][,1:l[1]], 15,cluster_method="tk_means")
s_perm = s[[2]]$cluster
PCA_result = fPCA(embedding_varySBM0, 1:80, 2, active_points)
PCA_score_show = PCA_result[[1]]
PCA_score_show = as.data.frame(cbind(PCA_score_show,0))
names(PCA_score_show) = c("x1", "x2", "col")
PCA_score_show$label = factor(C)
PCA_score_show$cluster = factor(s_perm)
PCA_score_show$mean_time = mean_time
PCA_score_show$active_duration = active_duration
score_plot0 = ggplot(PCA_score_show,aes(x=x1,y=x2,group=cluster,color=cluster,shape=label))+geom_point(size = 3)+geom_encircle(aes(group = cluster,fill=cluster),expand=0.03,spread=0.5,s_shape=1,size=3,linetype = 1,alpha=0.2)+scale_color_manual(values=c("#FFBE7A","#698aa2","#FC8677","#5383c3","#93b69c","#a59aca"))+scale_fill_manual(values=c("#FFBE7A","#698aa2","#f4a582","#5383c3","#93b69c","#a59aca"))+guides(fill="none")+theme_bw()
score_plot0_2 = ggplot(PCA_score_show,aes(x=x1,y=x2,group=cluster,color=mean_time,shape=label,size=active_duration))+geom_point()+scale_color_gradientn(values = scales::rescale(c(0,20,60,80), to = c(0, 1)), colors = c("#CCE5FF","lightskyblue","#5383c3","midnightblue"),breaks = c(20,40,60))+scale_size_continuous(range = c(2.5,4.5))+guides(size=guide_legend(override.aes = list(size=c(1.5,2,2.5,3,3.5))))+theme_bw()

## generalized DCCD plot:
PCA_result = fPCA_angle(embedding_varySBM0, 1:80, 100, active_points)
(l = which(cumsum(PCA_result[[3]])>0.8))
s = cluster_no_selection(PCA_result[[1]][,1:l[1]], 15,cluster_method="tk_means")
s_perm = s[[2]]$cluster
s_perm[s[[2]]$cluster==1] = 1
s_perm[s[[2]]$cluster==2] = 2
s_perm[s[[2]]$cluster==3] = 5
s_perm[s[[2]]$cluster==4] = 4
s_perm[s[[2]]$cluster==5] = 6
s_perm[s[[2]]$cluster==6] = 3
PCA_score_angle= fPCA_angle(embedding_varySBM0, 1:80, 2, active_points)
PCA_score_angle = as.data.frame(PCA_score_angle[[1]])
names(PCA_score_angle) = c("theta1", "theta2")
PCA_score_angle$label = factor(C)
PCA_score_angle$cluster = factor(s_perm)
PCA_score_angle$mean_time = mean_time
PCA_score_angle$active_duration = active_duration
angle_plot0 = ggplot(PCA_score_angle,aes(x=theta1,y=theta2,group=cluster,color=cluster,shape=label))+geom_point(size = 3)+geom_encircle(aes(group = cluster,fill=cluster),expand=0.03,spread=0.5,s_shape=1,size=3,linetype = 1,alpha=0.2)+scale_color_manual(values=c("#FC8677","#FFBE7A","#93b69c","#698aa2","#5383c3","#a59aca","#ea5506","#007b43","midnightblue","#5383c3","#7058a3"))+scale_fill_manual(values=c("#f4a582","#FFBE7A","#93b69c","#698aa2","#5383c3","#a59aca","#ea5506","#007b43","midnightblue","#5383c3","#7058a3"))+guides(fill="none")+theme_bw()
angle_plot0_2 = ggplot(PCA_score_angle,aes(x=theta1,y=theta2,group=cluster,color=mean_time,shape=label,size=active_duration))+geom_point()+scale_color_gradientn(values = scales::rescale(c(0,20,60,80), to = c(0, 1)), colors = c("#CCE5FF","lightskyblue","#5383c3","midnightblue"),breaks = c(20,40,60))+scale_size_continuous(range = c(2.5,4.5))+guides(size=guide_legend(override.aes = list(size=c(1.5,2,2.5,3,3.5))))+theme_bw()

s_perm2 = s_perm
s_perm2[16:30] = "2-1"
s_perm2[31:45] = "2-2"
s_perm2[76:90] = "5-1"
s_perm2[91:105] = "5-2"

## functions:
node_embedding_frame0 = function_sort(embedding_varySBM0, 1:T)
node_embedding_frame0 = cbind(node_embedding_frame0, label = factor(rep(C, each = T)), cluster = factor(rep(s_perm, each = T)), subcluster = factor(rep(s_perm2, each = T)))
func_plot_best_x1 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("x",1),group="i",color="cluster"))+facet_wrap("subcluster", nrow = 2)+geom_line(aes(alpha=cluster,size=cluster,linetype=label))+scale_color_manual(values=c("#FC8677","#FFBE7A","#93b69c","#698aa2","#5383c3","#a59aca","#ea5506","#007b43","midnightblue","#5383c3","#7058a3"))+scale_size_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+scale_alpha_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()+theme(strip.background = element_rect(fill = "white"))
func_plot_best_x2 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("x",2),group="i",color="cluster"))+facet_wrap("subcluster", nrow = 2)+geom_line(aes(alpha=cluster,size=cluster,linetype=label))+scale_color_manual(values=c("#FC8677","#FFBE7A","#93b69c","#698aa2","#5383c3","#a59aca","#ea5506","#007b43","midnightblue","#5383c3","#7058a3"))+scale_size_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+scale_alpha_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()+theme(strip.background = element_rect(fill = "white"))
func_plot_best_x3 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("x",3),group="i",color="cluster"))+facet_wrap("subcluster", nrow = 2)+geom_line(aes(alpha=cluster,size=cluster,linetype=label))+scale_color_manual(values=c("#FC8677","#FFBE7A","#93b69c","#698aa2","#5383c3","#a59aca","#ea5506","#007b43","midnightblue","#5383c3","#7058a3"))+scale_size_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+scale_alpha_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()+theme(strip.background = element_rect(fill = "white"))

layout <- "
AA
AA
BB
BB
"
plot_varySBM0_11 = (score_plot0+guides(color="none")+xlab(TeX("$x_1$"))+ylab(TeX("$x_2$"))+labs(caption = "(a) CD plot with estimated\n cluster for each node.")+
                    angle_plot0+xlab(TeX("$x_1^*$"))+ylab(TeX("$x_2^*$"))+labs(caption = "(c) DCCD plot with estimated\n cluster for each node.")+
                    plot_layout(guides = "collect",  design = layout, heights = c(4,4,4,4),width = c(4.5,7,7,3,3,7,7,4.5)))
plot_varySBM0_12 = (score_plot0_2+xlab(TeX("$x_1$"))+ylab(TeX("$x_2$"))+guides(shape=guide_legend(override.aes = list(size=3)))+labs(caption = "(b) CD plot with active\n duration for each node.")+
                    angle_plot0_2+xlab(TeX("$x_1^*$"))+ylab(TeX("$x_2^*$"))+guides(shape = "none")+labs(caption = "(d) DCCD plot with active\n duration for each node.")+
                    plot_layout(guides = "collect", design = layout, heights = c(4,4,4,4),width = c(4.5,7,7,3,3,7,7,4.5)))
layout <- "
ABCDE
"
plot_varySBM0_1 = (plot_spacer()+plot_varySBM0_11+plot_spacer()+plot_varySBM0_12+plot_spacer()) + plot_layout(design = layout, width = c(4.5,14,4,14,4.5)) & theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))


plot_varySBM0_3 = func_plot_best_x1+ylab(TeX("$Z_1$")) +labs(caption = "(a) the first dimension of the multivariate functional data.") + theme(axis.text=element_text(size=18), axis.title =element_text(size=20), strip.text =element_text(size=20), plot.caption =element_text(size=24, hjust=0.5), legend.text=element_text(size=22), legend.title =element_text(size=22), legend.background = element_rect(fill = 'white', colour = 'black'))
plot_varySBM0_4 = func_plot_best_x2+ylab(TeX("$Z_2$")) +labs(caption = "(b) the second dimension of the multivariate functional data.") + theme(axis.text=element_text(size=18), axis.title =element_text(size=20), strip.text =element_text(size=20), plot.caption =element_text(size=24, hjust=0.5), legend.text=element_text(size=22), legend.title =element_text(size=22), legend.background = element_rect(fill = 'white', colour = 'black'))
plot_varySBM0_2 = plot_varySBM0_3 / plot_varySBM0_4 +
  plot_layout(guides = "collect") & theme(axis.text=element_text(size=18), axis.title =element_text(size=20), plot.caption =element_text(size=24, hjust=0.5), legend.text=element_text(size=22), legend.title =element_text(size=22), legend.background = element_rect(fill = 'white', colour = 'black'))

ggsave(plot_varySBM0_1,file="varynodes_example_fASE0.pdf",width=27,height=15)
ggsave(plot_varySBM0_2,file="varynodes_example_fASE0_func.pdf",width=20,height=18)
