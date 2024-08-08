library(dplyr)
library(tidyr)
library(readr)
library(cluster)
library(doParallel)
library(RColorBrewer)
library(distr)
library(igraph)
library(ggplot2)
library(reshape2)
library(patchwork)
library(cowplot)
library(ggraph)
library(tidygraph)
library(flexclust)
library(PRROC)
library(mgsub)
library(scales)
library(lubridate)
library(ggpubr)
library(ggforce)

############  You may skip this procedure of data preparation.
load("primary_school.RData")
X = primaryschool%>%dplyr::select(i=2,j=3,time=1)%>%mutate(time=(time-31220))
X = as.data.frame(X)
X = mutate(X, time = as.numeric(time))
nodes = unique(c(X[,1],X[,2],recursive=TRUE))
node_num = length(nodes)
label = c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","Teachers")

index = rep(0, times = node_num)
for(q in 1:length(nodes)){
  node = nodes[q]
  index1 = which(primaryschool[,2]==node)
  index2 = which(primaryschool[,3]==node)
  if(length(index1)!=0){
    index1 = index1[1]
    index[q] = match(primaryschool[index1,4],label)
  }else{
    index2 = index2[1]
    index[q] = match(primaryschool[index2,5],label)
  }
}
XX = X
for(s in 1:length(nodes)){
  XX[X[,1]==nodes[s],1]=s
  XX[X[,2]==nodes[s],2]=s
}
X = XX

temp<-filter(X,i<j)
temp2<-filter(X,i>j)
temp3<-temp$i
temp$i<-temp$j
temp$j<-temp3
X<-rbind(temp,temp2)
T = max(X[,3])
T1 = 31080
T2 = 86020
T3 = T
T = T3-T2
X = filter(X, time<T1)

rm(index1, index2, q, s)

dynamic_network_adjacency = array(0, dim = c(242,242,26))
for(a in 1:26){
  temp = filter(X, (time>=(a-1)*1200) & (time<a*1200))
  temp2 = group_by(temp, i, j)%>%count()
  snapshot = igraph::graph_from_edgelist(as.matrix(temp2[,1:2]), directed = FALSE)
  E(snapshot)$weight = c(temp2[,3],recursive=TRUE)
  adjmat = as.matrix(get.adjacency(snapshot,type="both",attr="weight",sparse=FALSE))
  if(dim(adjmat)[1]<242){
    adjmat = cbind(rbind(adjmat,matrix(0, nrow = 242-dim(adjmat)[1], ncol = dim(adjmat)[1])),matrix(0, ncol = 242 - dim(adjmat)[1], nrow = 242))
  }
  dynamic_network_adjacency[,,a] = adjmat
}
## remove nodes that degree=0
equalzero = which(apply(dynamic_network_adjacency,1,sum)==0)
index = index[-equalzero]
dynamic_network_adjacency = dynamic_network_adjacency[-equalzero,-equalzero,]

adjacency_density = vector(mode = "list", 26)
for(t in 1:26){
  adjacency_density[[t]] = adjacency_visualize(dynamic_network_adjacency[,,t], index, permuted = TRUE, guide_break = c(6,12,18), n=238)
}
save(dynamic_network_adjacency, index, file = "primary_school_train.RData")

plot_primary_snap = adjacency_density[[6]]+ labs(caption = "(a) T=6, 10:05-10:25.") + guides(fill=FALSE)+
  adjacency_density[[14]]+ labs(caption = "(b) T=14, 13:05-13:25.") +
  adjacency_density[[22]]+ labs(caption = "(c) T=22, 16:05-16:25.") +  guides(fill=FALSE)+
  plot_layout(guides = "collect") &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'), plot.margin = unit(c(0.3, 1.5, 0.1, 0.15), "cm"))
ggsave(plot_primary_snap ,file="primary_snapshots.pdf",width=24,height=7.5)
####################

############ You may skip these costly procedures of embedding.
tuning_primary = fASE_tuning(5:10, 1:5, tuning_method = "heuristic", adjacency_tensor = dynamic_network_adjacency, spline_order = 4, batch_size = NULL, kernel_scale = FALSE , timestamp_vec = NULL, scalable_method = TRUE, scalable_dim = 20, scalable_power = 6, iteration_method = "sequential", iteration_step = 2500, epsilon=1e-6)
tuning_primary = get_tuning(tuning_primary)

result_primary = fASE(dynamic_network_adjacency, tuning_primary[[2]], tuning_primary[[1]], batch_size = NULL, timestamp_vec = 1:26, scalable_dim = 20, scalable_power = 6, iteration_step = 2500, epsilon=1e-6)
save(file="primary_fASE.RData",tuning_primary, result_primary)
####################

source("DataAndCode\\functional_ASE.R")
source("DataAndCode\\downstream_tasks.R")
load("primary_fASE.RData")
load("primary_school_train.RData")
T = 26
embedding_primary = eval.fd(1:T, result_primary[[1]])

set.seed(1)

## clusters obtained by radius-corrected MVFPCA.
PCA_angle = fPCA_angle(embedding_primary,  1:26, 130)
(l = which(cumsum(PCA_angle[[3]])>0.8))
s = cluster_no_selection(PCA_angle[[1]][,1:l[1]], 15,cluster_method = "tk_means",K_selection = "Silhouette")
library(aricode)
(nmi_primary = NMI(s[[2]]$cluster, index))
# 0.8321506

## clusters obtained by MVFPCA.
PCA_score = fPCA(embedding_primary,  1:26, 130)
(l = which(cumsum(PCA_score[[3]])>0.8))
u = cluster_no_selection(PCA_score[[1]][,1:l[1]], 15,cluster_method = "tk_means",K_selection = "Silhouette")

## visualization of clusters obtained by radius-corrected MVFPCA.
PCA_cluster = data.frame(nodes = 1:dim(dynamic_network_adjacency)[1], group = index, cluster_tmp = s[[2]]$cluster, cluster = 0)
PCA_label = PCA_cluster%>%group_by(cluster_tmp)%>%count(group)%>%as.data.frame()
PCA_sub = rep(0,times = 10)
for(h in 1:10){
  PCA_sub[h] = PCA_label[PCA_label$group==h,"cluster_tmp"][which.max(PCA_label[PCA_label$group==h,"n"])]
}
PCA_sub2 = rep(0,times = max(s[[2]]$cluster))
PCA_sub2[PCA_sub] = 1:10
for(h in 1:max(s[[2]]$cluster)){
  if(PCA_sub2[h]==0){
    PCA_sub2[h] = PCA_label[PCA_label$cluster_tmp==h,"group"][which.max(PCA_label[PCA_label$cluster_tmp==h,"n"])]
  }
}
PCA_index = order(PCA_sub2,decreasing=FALSE)
for(h in 1:max(s[[2]]$cluster)){
  PCA_cluster[PCA_cluster$cluster_tmp==PCA_index[h],"cluster"] = h
}
PCA_cluster$cluster = factor(PCA_cluster$cluster)
PCA_cluster$group = factor(PCA_cluster$group)
img_cluster = ggplot(PCA_cluster,y=0:1,aes(x=cluster)) +geom_bar(stat = 'count',aes(fill=group,width = 0.5, position = "fill"))+
  scale_x_discrete()+
  scale_fill_manual(values = c("#C82423","#F89100","#FFC94C", "#B7CD4D","#A5CE9D","forestgreen","lightskyblue","#5383c3","#82B0D2","#ccccff","#D3D3D3","#33001A"))+
  theme_bw()+ theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black')) +
  scale_y_continuous(position = "left", expand=c(0,0.5))+
  labs(y=NULL)

## CD plot and DCCD plot.
PCA_nonangle = fPCA(embedding_primary,  1:26, 2)
PCA_score_nonangle=PCA_nonangle[[1]]
PCA_score_nonangle = as.data.frame(PCA_score_nonangle)
names(PCA_score_nonangle) = c("x1", "x2")
PCA_score_nonangle$cluster = factor(u[[2]]$cluster)
nonangle_plot = ggplot(PCA_score_nonangle,aes(x=x1,y=x2,group=cluster,color=cluster))+geom_point(size=2)+theme_bw()+
  geom_encircle(aes(group = cluster,fill=cluster),expand=0,spread=0.5,s_shape=1,size=3,linetype = 1,alpha=0.4)+
  scale_color_manual(values = c("#698aa2","#FFBE7A","#FC8677","#82B0D2", "#B9D993", "#A5CE9D", "#B9BED7", "#5383c3",  "#D8A39C", "#D78C9F",  "#FFD9E6", "#ccccff", "#D3D3D3", "#D6ecf0", "#d3b17d" ))+
  scale_fill_manual(values = c("#698aa2","#FFBE7A","#FC8677","#82B0D2", "#B9D993", "#A5CE9D", "#B9BED7", "#5383c3",  "#D8A39C", "#D78C9F",  "#FFD9E6", "#ccccff", "#D3D3D3", "#D6ecf0", "#d3b17d"))+
  guides(fill="none",color=guide_legend(ncol=2))+
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
PCA_angle = fPCA_angle(embedding_primary,  1:26, 2)
PCA_score_angle=PCA_angle[[1]]
PCA_score_angle = as.data.frame(PCA_score_angle)
names(PCA_score_angle) = c("theta1", "theta2")
PCA_score_angle$class = factor(index)
PCA_score_angle$cluster = factor(s[[2]]$cluster)
library(ggalt)
angle_plot2 = ggplot(PCA_score_angle,aes(x=theta1,y=theta2,group=cluster,color=cluster))+geom_point(size=2)+theme_bw()+
  geom_encircle(aes(group = cluster,fill=cluster),expand=0,spread=0.5,s_shape=1,size=3,linetype = 1,alpha=0.4)+
  scale_color_manual(values = c("#698aa2","#FFBE7A","#FC8677","#82B0D2", "#B9D993",  "#B9BED7",   "#D8A39C", "#D78C9F",  "#FFD9E6"))+
  scale_fill_manual(values = c("#698aa2","#FFBE7A","#FC8677","#82B0D2", "#B9D993",  "#B9BED7",   "#D8A39C", "#D78C9F",  "#FFD9E6"))+
  guides(fill="none")+
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
angle_plot = nonangle_plot +xlab(TeX("$x_1$"))+ylab(TeX("$x_2$"))+ labs(caption = "(a) CD plot with estimated community.") + angle_plot2 +xlab(TeX("$x_1^*$"))+ylab(TeX("$x_2^*$")) + labs(caption = "(b) DCCD plot with estimated community.") + img_cluster + labs(caption  = "(c) the components of estimated\n community in DCCD plot.")
ggsave(angle_plot ,file="fPCA_plot.pdf",width=27,height=7.5)

