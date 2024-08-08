library(fda)
library(mrfDepth)
library(MFPCA)
library(tidygraph)
library(dplyr)
library(ggplot2)
library(igraph)
library(hawkes)
library(tidyr)
library(poweRlaw)
library(aricode)
library(mrfDepth)
library(DepthProc)
library(Rtsne)
library(viridis)
library(ClusterR)
library(patchwork)
source("DataAndCode\\functional_ASE.R")
source("DataAndCode\\downstream_tasks.R")

############  You may skip this costly procedure of data preparation.
library(R.matlab)
CiteeDynamicFinal <- readMat("statistical_coauthorship_citation network/Ready-to-use data matrice/Co-citation networks/CiteeDynamicFinal.mat")
authorname_all <- read_delim("statistical_coauthorship_citation network/Ready-to-use data matrice/Co-citation networks/author_name.txt",
                             delim = "\t", escape_double = FALSE,
                             col_names = FALSE, trim_ws = TRUE)

papers = filter(AuPapMat, journal%in%c("Bka","JRSSB", "JASA", "AoS"), year>=1976)
paper_included = papers$idxPap
authors_included = unique(papers$idxAu)
authors_included2 = authors_included
authorname = c(authorname_all,recursive=TRUE)[authors_included]
citation_included = filter(PapPapMat, ToPap%in%paper_included, FromYear>=1976)
citation_included = citation_included%>%arrange(FromYear, FromPap)

acc_network_list = c()
h = 0
for(paper_index in unique(citation_included$FromPap)){
  h  = h+1
  all_cited_paper = citation_included[citation_included$FromPap==paper_index,"ToPap"]

  fromauthor = match(papers[papers$idxPap==paper_index,"idxAu"],authors_included)
  toauthor = na.omit(match(papers[papers$idxPap%in%c(all_cited_paper),"idxAu"],authors_included))
  if(length(toauthor)!=0){
    citation_Dframe = cbind(rep(fromauthor,times=length(toauthor)), rep(toauthor,each=length(fromauthor)))
    acc_network_list = rbind(acc_network_list, citation_Dframe)
  }
}
save(acc_network_list, file = "citation_network_temp.RData")

colnames(acc_network_list) = c("out_nodes", "in_nodes")
acc_network_list = as.data.frame(acc_network_list)%>%group_by(out_nodes,in_nodes)%>%count()
accnetwork = igraph::graph_from_edgelist(as.matrix(acc_network_list[,1:2]), directed = FALSE)
E(accnetwork)$weight = acc_network_list$n
accnetwork_adj = igraph::as_adjacency_matrix(accnetwork,attr="weight")
deg_accnetwork =igraph::degree(accnetwork)
h = 0
while(any(deg_accnetwork<30)){
  h = h+1
  authors_included = authors_included[deg_accnetwork>=30]
  accnetwork_adj = accnetwork_adj[deg_accnetwork>=30,deg_accnetwork>=30]
  deg_accnetwork = apply(accnetwork_adj,1,sum)
}
authorname = c(authorname_all,recursive=TRUE)[authors_included]

# 2 years a snapshot
acc_dynnetwork_list = c()
h = 0
for(paper_index in unique(citation_included$FromPap)){
  cited_year = citation_included[citation_included$FromPap==paper_index,"FromYear"][1]
  t = ceiling((cited_year-1975)/2)
  h  = h+1
  all_cited_paper = citation_included[citation_included$FromPap==paper_index,"ToPap"]

  fromauthor = na.omit(match(papers[papers$idxPap==paper_index,"idxAu"],authors_included))
  toauthor = na.omit(match(papers[papers$idxPap%in%c(all_cited_paper),"idxAu"],authors_included))

  if(length(toauthor)!=0 && length(fromauthor)!=0){
    citation_Dframe = cbind(rep(fromauthor,times=length(toauthor)), rep(toauthor,each=length(fromauthor)), t)
    acc_dynnetwork_list = rbind(acc_dynnetwork_list, citation_Dframe)
  }
}

rm(list = ls())
gc()
load("citation_network_temp.RData")
dynamic_network = array(0, dim = c(length(authorname),length(authorname),20) )
acc_dynnetwork_list = as.data.frame(acc_dynnetwork_list)
for(T in 1:20){
  acc_snap_list = acc_dynnetwork_list%>%dplyr::filter(t==T)%>%rename("out_nodes"=1, "in_nodes"=2, "t"=3)
  acc_snap_list = as.data.frame(acc_snap_list)%>%group_by(out_nodes,in_nodes)%>%count()
  accnetwork = igraph::graph_from_edgelist(as.matrix(acc_snap_list[,1:2]), directed = FALSE)
  E(accnetwork)$weight = acc_snap_list$n
  adj_snap = igraph::as_adjacency_matrix(accnetwork,attr="weight")
  if(dim(adj_snap)[1]<length(authorname)){
    adj_snap = cbind(rbind(adj_snap, matrix(0,nrow=length(authorname)-dim(adj_snap)[1],ncol=dim(adj_snap)[1])),
                     matrix(0,nrow=length(authorname),ncol=length(authorname)-dim(adj_snap)[1]))
  }
  dynamic_network[,,T] = as.matrix(adj_snap)
}
save(acc_dynnetwork_list, dynamic_network, authorname, file = "citation_data.RData")
############################

############  You may skip this costly procedure of embedding.
set.seed(1)
t1 = Sys.time()
result = fASE(dynamic_network_adjacency, 10, 10, timestamp_vec = 1:20, scalable_power = 5, scalable_dim = 30, epsilon = 1e-6, iteration_step = 10000, step_size = 0.1)
t2 = Sys.time()
t = as.numeric(t2-t1, units = "mins")
save(file = "citation_result.RData",result,t)
############################

load("citation_data.RData")
load("citation_result.RData")

embedding_result_fASE = eval.fd(1:20, result[[1]])

## Calculate the depth of the whole network:
sim=cosine_similarity(eval.fd(1:20, result[[1]]), 1:20)
set.seed(1)
gc()
tsne_result <- Rtsne(1-sim, dims = 2, perplexity = round(dim(dynamic_network)[1]/4), is_distance=TRUE)
tsne_embeddings <- tsne_result$Y
tsne_embeddings=as.data.frame(tsne_embeddings)
colnames(tsne_embeddings) = c("x1","x2")
rownames(tsne_embeddings) = 1:dim(dynamic_network)[1]
tsne_pos = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="label"))+geom_point()
depth_stats = metric_depth(1-sim)
depth_reindex = order(depth_stats,decreasing =TRUE)
tsne_embeddings$degree = apply(dynamic_network,1,sum)
tsne_embeddings$depth_order = 0
tsne_embeddings$name = 0
tsne_embeddings[,"depth"] = depth_stats
tsne_embeddings[depth_reindex,"depth_order"] = 1:dim(tsne_embeddings)[1]
tsne_embeddings[,"name"] = authorname
tsne_embeddings[,"name_label"] = paste0(tsne_embeddings[,"name"] , ",", tsne_embeddings[,"depth_order"])

set.seed(1)
tsne_pos = ggplot(tsne_embeddings,aes(x=depth,y=degree,color=depth*log(degree)))+scale_y_log10()+geom_point()+
  geom_hline(yintercept = 1100)+
  annotation_custom(
    grob = grid::linesGrob(),
    xmin = 0.463, xmax = 0.463, ymin = -Inf, ymax = Inf)+
  ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:60)|degree>1100), aes(label = name),size=7,force=60,max.iter=1000000,)+
  scale_color_gradientn(values = scales::rescale(c(1,3,4), to = c(0, 1)),colors = c("lightskyblue","dodgerblue","midnightblue"))+
  theme_bw()+theme(text = element_text(size = 28),axis.title = element_text(size = 24),legend.text = element_text(size = 24))+labs(color = "depth*\nlog(degree)")+
  theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'), legend.key.height = unit(2, "cm"))
ggsave(file=paste0("citation_networks/all_cendeg.pdf"), tsne_pos, height = 18, width = 24)
save(sim, depth_stats, tsne_embeddings, file = "citation_networks/all_cendeg.RData")

## cluster with radius-corrected MVFPCA:
############  You may skip this costly procedure of clustering.
set.seed(1)
active_points = lapply(1:length(authorname), active_calculation, dynamic_network)
active_points = lapply(active_points, function(x){x[1]:20})
PCA_angle = fPCA_angle(embedding_result_fASE,  1:20, 200, active_points)
(l = which(cumsum(PCA_angle[[3]])>0.8))
s = cluster_no_selection(PCA_angle[[1]][,1:l[1]], 50,cluster_method = "tk_means",K_selection = "Silhouette")
save(file = "citation_result.RData",result,t,s)
#########################

## generalized DCCD plot:
set.seed(1)
basic_colors = c("#992f87", "#c9a2c6", "pink","#ef7a82","#FF0000","#F89100","#FFDA9F","forestgreen","#bee6ff","#82B0D2","#126085")
color_palette <- c()
for(i in 1:11){
  color_palette <- c(color_palette,colorRampPalette(c(basic_colors[i], basic_colors[i+1]))(4)[1:3])
}
active_points = lapply(1:length(authorname), active_calculation, dynamic_network)
active_points = lapply(active_points, function(x){x[1]:20})
PCA_angle = fPCA_angle(embedding_result_fASE,  1:20, 2, active_points)
PCA_score_angle=PCA_angle[[1]]
PCA_score_angle = as.data.frame(PCA_score_angle)
names(PCA_score_angle) = c("theta1", "theta2")
PCA_score_angle$theta1=-PCA_score_angle$theta1
PCA_score_angle$cluster = factor(s[[2]]$cluster)
mass_center1 = tapply(PCA_score_angle$theta1, PCA_score_angle$cluster, mean)
mass_center2 = tapply(PCA_score_angle$theta2, PCA_score_angle$cluster, mean)
mass_center = as.data.frame(cbind(mass_center1,mass_center2))
mass_distance = as.matrix(dist(mass_center))
col_indices = rep(0, times = dim(mass_center)[1])
col_indices[21] <- 1
cluster_now = 21
cluster_left = c(1:20,22:dim(mass_center)[1])
h = 1
while(length(cluster_left)>0) {
  h = h+1
  closest_cluster = cluster_left[which.min(mass_distance[cluster_now,cluster_left])]
  col_indices[closest_cluster] <- h
  cluster_now = closest_cluster
  cluster_left = setdiff(cluster_left, closest_cluster)
}
color_palette = append(setdiff(color_palette,"#AAD4F0"),"#663D74",after = 0)
col_indices[col_indices==26] = 0  ## change the order of clusters
col_indices[col_indices<26] = col_indices[col_indices<26]+1
PCA_score_angle$cluster = factor(col_indices[s[[2]]$cluster])
PCA_score_angle$name = authorname
PCA_score_angle$degree = apply(dynamic_network,1,sum)
PCA_score_angle$degree_order = 0
PCA_score_angle$degree_order[order(PCA_score_angle$degree, decreasing=TRUE)] = 1:dim(PCA_score_angle)[1]
angle_plot = ggplot(PCA_score_angle,aes(x=theta1,y=theta2,group=cluster,color=cluster))+geom_point()+xlab(TeX("$x_1^*$"))+ylab(TeX("$x_2^*$"))+
  ggrepel::geom_label_repel(data = subset(PCA_score_angle, degree_order%in%c(1:60)), aes(label = name),size=5,force=50,show.legend=FALSE)+
  scale_color_manual(values=color_palette)+theme_bw()+
  theme(legend.position = "none", axis.text=element_text(size=28), axis.title =element_text(size=28),  plot.caption =element_text(size=30, hjust=0.5))
ggsave(file=paste0("fPCA_results_clusters.pdf"), angle_plot, height = 12, width = 18)

PCA_score_angle$first_active_year = seq(1976,2015,by=2)[sapply(active_points, function(x){x[[1]]}, simplify = TRUE)]
PCA_score_angle$mean_time = seq(1976,2015,by=2)[sapply(active_points, function(x){mean(x)}, simplify = TRUE)]
start_time_plot = ggplot(PCA_score_angle,aes(x=theta1,y=theta2,color=mean_time))+geom_point()+
  ggrepel::geom_label_repel(data = subset(PCA_score_angle, degree_order%in%c(1:60)), aes(label = name),size=5,force=50)+
  scale_color_gradientn(values = scales::rescale(c(1976, 1985, 1995, 2000, 2015), to = c(0, 1)),colors = c("midnightblue","royalblue","#5383c3","lightskyblue","#CCE5FF"))+theme_bw()+
  theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(file=paste0("fPCA_results_withT.pdf"), start_time_plot, height = 12, width = 18)

## tSNE embedding:
tsne_embeddings[,"cluster"] = PCA_score_angle$cluster
tsne_embeddings[,"first_active_year"] = PCA_score_angle$first_active_year
tsne_embeddings$cluster = factor(tsne_embeddings$cluster)
tsne_pos2 = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="cluster"))+geom_point()+
  ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:30)), aes(label = name_label),size=7,force=50)+
  scale_color_manual(values = color_palette)+
  theme_bw()+
  theme(axis.text=element_text(size=28), axis.title =element_text(size=28),  plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(file=paste0("allwithcluster.pdf"), tsne_pos2, height = 18, width = 24)
tsne_pos3 = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="first_active_year"))+geom_point()+
  ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:30)), aes(label = name_label),size=7,force=50)+
  scale_color_gradientn(values = scales::rescale(c(1,round(dim(tsne_embeddings)[1]/2),dim(tsne_embeddings)[1]), to = c(0, 1)),colors = c("lightskyblue","dodgerblue","midnightblue"))+
  theme_bw()+
  theme(axis.text=element_text(size=28), axis.title =element_text(size=28), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(file=paste0("allwithactiveyear.pdf"), tsne_pos3, height = 18, width = 24)


## time-cluster plot:
tsne_embeddings$depth = depth_stats
tsne_embeddings_group = group_by(tsne_embeddings, cluster)%>%summarise(med_first_active_year = median(first_active_year), med_depth = median(depth))
tsne_embeddings_group$med_first_active_time = tsne_embeddings_group$med_first_active_year
tsne_embeddings_group$med_first_active_year = as.POSIXct(paste0(tsne_embeddings_group$med_first_active_year, "-01-01"), format = "%Y-%m-%d")
tsne_embeddings_group$med_first_active_width = as.POSIXct("2016-01-01", format = "%Y-%m-%d")-tsne_embeddings_group$med_first_active_year
tsne_embeddings_group$group_name = c("Likelihood Theory","Nonparametric I", "Nonparametric II", "Nonparametric III", "Wavelet",
                                     "Measurement Error", "High Dimension I", "High Dimension II", "Multiple Tests",
                                     "Time Series", "High Dimension III", "Functional Data", "Stochastic Ordering",
                                     "Bayesian I", "Computer Experiments", "Robust Statistics", "Nonparametric IV",
                                     "Nonparametric V", "High Dimension IV", "Statistical Learning", "Bayesian II",
                                     "Bayesian III", "Bayesian IV", "Bayesian V", "Bayesian VI",
                                     "Biostatistics I", "Biostatistics II", "Biostatistics III",
                                     "Biostatistics IV", "Biostatistics V", "Biostatistics VI", "Biostatistics VII")
tsne_embeddings_group$color = color_palette[1:32]
tsne_embeddings_group$cluster =  as.numeric(as.character(tsne_embeddings_group$cluster))
tsne_embeddings_group = tsne_embeddings_group[c(16,1,2,14,18,9,13,5,6,15,17,3,20,12,11,10,4,19,7,8,23,25,24,22,21,29,31,27,32,26,28,30),]
tsne_embeddings_group$cluster = 1:32
tsne_embeddings_group$biggroup_name = c("Non/Semi & HD I", "Non/Semi & HD II", "Non/Semi & HD III", "Non/Semi & HD IV",
                                        "Non/Semi & HD V", "Non/Semi & HD VI", "Non/Semi & HD VII", "Non/Semi & HD VIII",
                                        "Non/Semi & HD IX", "Non/Semi & HD X", "Non/Semi & HD XI", "Non/Semi& HD XII",
                                        "Non/Semi & HD XIII", "Non/Semi & HD XVI", "Non/Semi & HD XV", "Non/Semi & HD XVI",
                                        "Non/Semi & HD XVII", "Non/Semi & HD XVIII", "Non/Semi & HD XIX", "Non/Semi & HD XX", "Bayesian I",
                                        "Bayesian II", "Bayesian III", "Bayesian IV", "Bayesian V",
                                        "Biostatistics I", "Biostatistics II", "Biostatistics III",
                                        "Biostatistics IV", "Biostatistics V", "Biostatistics VI", "Biostatistics VII")
cluster_timeplot1 = ggplot(tsne_embeddings_group, aes(fill = factor(cluster), color = factor(cluster))) +
  geom_rect(aes(xmin = med_first_active_year, xmax = med_first_active_year+med_first_active_width, ymin = cluster-0.3, ymax = cluster+0.3)) +  # 使用geom_tile绘制色块
  scale_x_datetime(name = "Year", date_labels = "%Y", limits = c(as.POSIXct("1974-01-01", format = "%Y-%m-%d"),as.POSIXct("2016-01-01", format = "%Y-%m-%d"))) +  # 格式化X轴时间显示
  geom_text(aes(x = med_first_active_year - 3600*24*1800, y = cluster, label = biggroup_name), size = 12) +
  labs(y = "cluster") + guides(fill = FALSE, color = FALSE) + # Y axis label
  scale_fill_manual(values = tsne_embeddings_group$color)+
  scale_color_manual(values = tsne_embeddings_group$color)+
  theme_bw() +
  theme(axis.text=element_text(size=28), axis.title =element_text(size=28), axis.text.y = element_blank(), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(file=paste0("clusteractiveyear.pdf"), cluster_timeplot1, height = 18, width = 24)


## Dynamic performance of clusters
set.seed(1)
active_points = lapply(1:length(authorname), active_calculation, dynamic_network)
active_nodes = which(!is.na(lapply(active_points, function(x){if(x[1]>5){NA}else{x[1]}})))
active_points = lapply(active_points, function(x){x[1]:5})[active_nodes]
tsne_embeddings_first10Y = embedding_result_fASE[1:5, active_nodes,]
authorname_first10Y = authorname[active_nodes]
PCA_angle = fPCA_angle(tsne_embeddings_first10Y,  1:5, 200, active_points)
which(cumsum(PCA_angle[[3]])>0.8)
l = which(cumsum(PCA_angle[[3]])>0.8)[1]
s_first10Y = cluster_no_selection(PCA_angle[[1]][,1:l], 20,cluster_method = "tk_means")
save(s_first10Y, authorname_first10Y, tsne_embeddings_first10Y, file = "citation_networks/dynamic_clustering_first10Y.RData")
edges_first10Y = c()
for(i in 1:max(s_first10Y[[2]]$cluster)){
  edges_first10Y = c(edges_first10Y, median(apply(dynamic_network[active_nodes, ,1:5][s_first10Y[[2]]$cluster==i,,],1,sum)))
}
set.seed(1)
for(i in 1:max(s_first10Y[[2]]$cluster)){
  nodes_i = which(s_first10Y[[2]]$cluster%in%c(i))
  sim=cosine_similarity(eval.fd(seq(1,5,length=100), result[[1]][active_nodes][nodes_i]), seq(1,5,length=100))
  tsne_result <- Rtsne(1-sim, dims = 2, perplexity = round(length(nodes_i)/4), is_distance=TRUE)
  tsne_embeddings <- tsne_result$Y
  tsne_embeddings=as.data.frame(tsne_embeddings)
  colnames(tsne_embeddings) = c("x1","x2")
  rownames(tsne_embeddings) = nodes_i
  tsne_pos = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="label"))+geom_point()
  depth_stats = metric_depth(1-sim)
  # depth_stats = depth(tsne_embeddings,tsne_embeddings)
  depth_reindex = order(depth_stats,decreasing =TRUE)
  tsne_embeddings$depth_order = 0
  tsne_embeddings$name = 0
  tsne_embeddings[depth_reindex,"depth_order"] = 1:dim(tsne_embeddings)[1]
  tsne_embeddings[,"name"] = authorname_first10Y[nodes_i]
  tsne_embeddings[,"name_label"] = paste0(tsne_embeddings[,"name"] , ",", tsne_embeddings[,"depth_order"])
  if(i==1){
    tsne_pos1 = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="depth_order"))+geom_point()+
      ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1,3:20)), aes(label = name_label),size=2.5,force=50)+
      scale_color_gradientn(values = scales::rescale(c(1,round(dim(tsne_embeddings)[1]/2),dim(tsne_embeddings)[1]), to = c(0, 1)),colors = c("lightskyblue","dodgerblue","midnightblue"))+
      theme_bw()
  }else{
    tsne_pos1 = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="depth_order"))+geom_point()+
      ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:20)), aes(label = name_label),size=2.5,force=50)+
      scale_color_gradientn(values = scales::rescale(c(1,round(dim(tsne_embeddings)[1]/2),dim(tsne_embeddings)[1]), to = c(0, 1)),colors = c("lightskyblue","dodgerblue","midnightblue"))+
      theme_bw()
  }
  # regression? dimension reduction? biostatistics?
  ggsave(file=paste0("citation_networks/first10Y/group",i,".pdf"), tsne_pos1, height = 9, width = 12)
}

set.seed(1)
active_points = lapply(1:length(authorname), active_calculation, dynamic_network)
active_nodes = which(!is.na(lapply(active_points, function(x){if(x[1]>10){NA}else{x[1]}})))
active_points = lapply(active_points, function(x){max(x[1],6):10})[active_nodes]
tsne_embeddings_mid10Y = embedding_result_fASE[6:10, active_nodes,]
authorname_mid10Y = authorname[active_nodes]
PCA_angle = fPCA_angle(tsne_embeddings_mid10Y,  6:10, 200, active_points)
which(cumsum(PCA_angle[[3]])>0.8)
l = which(cumsum(PCA_angle[[3]])>0.8)[1]
s_mid10Y = cluster_no_selection(PCA_angle[[1]][,1:l], 20,cluster_method = "tk_means")
save(s_mid10Y, file = "citation_networks/dynamic_clustering_mid10Y.RData")
edges_mid10Y = c()
for(i in 1:max(s_mid10Y[[2]]$cluster)){
  edges_mid10Y = c(edges_mid10Y, median(apply(dynamic_network[active_nodes, ,6:10][s_mid10Y[[2]]$cluster==i,,],1,sum)))
}
set.seed(1)
for(i in 1:max(s_mid10Y[[2]]$cluster)){
  nodes_i = which(s_mid10Y[[2]]$cluster%in%c(i))
  sim=cosine_similarity(eval.fd(seq(6,10,length=100), result[[1]][active_nodes][nodes_i]), seq(6,10,length=100))
  tsne_result <- Rtsne(1-sim, dims = 2, perplexity = round(length(nodes_i)/4), is_distance=TRUE)
  tsne_embeddings <- tsne_result$Y
  tsne_embeddings=as.data.frame(tsne_embeddings)
  colnames(tsne_embeddings) = c("x1","x2")
  rownames(tsne_embeddings) = nodes_i
  tsne_pos = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="label"))+geom_point()
  depth_stats = metric_depth(1-sim)
  # depth_stats = depth(tsne_embeddings,tsne_embeddings)
  depth_reindex = order(depth_stats,decreasing =TRUE)
  tsne_embeddings$depth_order = 0
  tsne_embeddings$name = 0
  tsne_embeddings[depth_reindex,"depth_order"] = 1:dim(tsne_embeddings)[1]
  tsne_embeddings[,"name"] = authorname_mid10Y[nodes_i]
  tsne_embeddings[,"name_label"] = paste0(tsne_embeddings[,"name"] , ",", tsne_embeddings[,"depth_order"])
  tsne_pos1 = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="depth_order"))+geom_point()+
    ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:20)), aes(label = name_label),size=2.5,force=50)+
    scale_color_gradientn(values = scales::rescale(c(1,round(dim(tsne_embeddings)[1]/2),dim(tsne_embeddings)[1]), to = c(0, 1)),colors = c("lightskyblue","dodgerblue","midnightblue"))+
    theme_bw()
  # regression? dimension reduction? biostatistics?
  ggsave(file=paste0("citation_networks/mid10Y/group",i,".pdf"), tsne_pos1, height = 9, width = 12)
}

set.seed(1)
active_points = lapply(1:length(authorname), active_calculation, dynamic_network)
active_nodes = which(!is.na(lapply(active_points, function(x){if(x[1]>15){NA}else{x[1]}})))
active_points = lapply(active_points, function(x){max(x[1],11):15})[active_nodes]
tsne_embeddings_midd10Y = embedding_result_fASE[11:15, active_nodes,]
authorname_midd10Y = authorname[active_nodes]
PCA_angle = fPCA_angle(tsne_embeddings_midd10Y,  11:15, 200, active_points)
which(cumsum(PCA_angle[[3]])>0.8)
l = which(cumsum(PCA_angle[[3]])>0.8)[1]
s_midd10Y = cluster_no_selection(PCA_angle[[1]][,1:l], 20,cluster_method = "tk_means")
save(s_midd10Y, file = "citation_networks/dynamic_clustering_midd10Y.RData")
edges_midd10Y = c()
for(i in 1:max(s_midd10Y[[2]]$cluster)){
  edges_midd10Y = c(edges_midd10Y, median(apply(dynamic_network[active_nodes, ,11:15][s_midd10Y[[2]]$cluster==i,,],1,sum)))
}
set.seed(1)
for(i in 1:max(s_midd10Y[[2]]$cluster)){
  nodes_i = which(s_midd10Y[[2]]$cluster%in%c(i))
  sim=cosine_similarity(eval.fd(seq(6,10,length=100), result[[1]][active_nodes][nodes_i]), seq(6,10,length=100))
  tsne_result <- Rtsne(1-sim, dims = 2, perplexity = round(length(nodes_i)/4), is_distance=TRUE)
  tsne_embeddings <- tsne_result$Y
  tsne_embeddings=as.data.frame(tsne_embeddings)
  colnames(tsne_embeddings) = c("x1","x2")
  rownames(tsne_embeddings) = nodes_i
  tsne_pos = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="label"))+geom_point()
  depth_stats = metric_depth(1-sim)
  # depth_stats = depth(tsne_embeddings,tsne_embeddings)
  depth_reindex = order(depth_stats,decreasing =TRUE)
  tsne_embeddings$depth_order = 0
  tsne_embeddings$name = 0
  tsne_embeddings[depth_reindex,"depth_order"] = 1:dim(tsne_embeddings)[1]
  tsne_embeddings[,"name"] = authorname_midd10Y[nodes_i]
  tsne_embeddings[,"name_label"] = paste0(tsne_embeddings[,"name"] , ",", tsne_embeddings[,"depth_order"])
  tsne_pos1 = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="depth_order"))+geom_point()+
    ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:20)), aes(label = name_label),size=2.5,force=50)+
    scale_color_gradientn(values = scales::rescale(c(1,round(dim(tsne_embeddings)[1]/2),dim(tsne_embeddings)[1]), to = c(0, 1)),colors = c("lightskyblue","dodgerblue","midnightblue"))+
    theme_bw()
  # regression? dimension reduction? biostatistics?
  ggsave(file=paste0("citation_networks/midd10Y/group",i,".pdf"), tsne_pos1, height = 9, width = 12)
}

set.seed(1)
active_points = lapply(1:length(authorname), active_calculation, dynamic_network)
active_points = lapply(active_points, function(x){max(x[1],16):20})
tsne_embeddings_final10Y = embedding_result_fASE[16:20, ,]
authorname_final10Y = authorname
PCA_angle = fPCA_angle(tsne_embeddings_final10Y,  16:20, 200, active_points)
which(cumsum(PCA_angle[[3]])>0.8)
l = which(cumsum(PCA_angle[[3]])>0.8)[1]
s_final10Y = cluster_no_selection(PCA_angle[[1]][,1:l], 20,cluster_method = "tk_means")
save(s_final10Y, file = "citation_networks/dynamic_clustering_final10Y.RData")
sim=cosine_similarity(eval.fd(16:20, result[[1]]), 16:20)
depth_stats = metric_depth(1-sim)
depth_reindex = order(depth_stats,decreasing =TRUE)
depth_final10Y = tapply(depth_stats,s_final10Y[[2]]$cluster,median)
edges_final10Y = c()
for(i in 1:max(s_final10Y[[2]]$cluster)){
  edges_final10Y = c(edges_final10Y, median(apply(dynamic_network[, ,16:20][s_final10Y[[2]]$cluster==i,,],1,sum)))
}
set.seed(1)
for(i in 1:max(s_final10Y[[2]]$cluster)){
  nodes_i = which(s_final10Y[[2]]$cluster%in%c(i))
  sim=cosine_similarity(eval.fd(seq(6,10,length=100), result[[1]][nodes_i]), seq(6,10,length=100))
  tsne_result <- Rtsne(1-sim, dims = 2, perplexity = round(length(nodes_i)/4), is_distance=TRUE)
  tsne_embeddings <- tsne_result$Y
  tsne_embeddings=as.data.frame(tsne_embeddings)
  colnames(tsne_embeddings) = c("x1","x2")
  rownames(tsne_embeddings) = nodes_i
  tsne_pos = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="label"))+geom_point()
  depth_stats = metric_depth(1-sim)
  # depth_stats = depth(tsne_embeddings,tsne_embeddings)
  depth_reindex = order(depth_stats,decreasing =TRUE)
  tsne_embeddings$depth_order = 0
  tsne_embeddings$name = 0
  tsne_embeddings[depth_reindex,"depth_order"] = 1:dim(tsne_embeddings)[1]
  tsne_embeddings[,"name"] = authorname_final10Y[nodes_i]
  tsne_embeddings[,"name_label"] = paste0(tsne_embeddings[,"name"] , ",", tsne_embeddings[,"depth_order"])
  tsne_pos1 = ggplot(tsne_embeddings,aes_string(x="x1",y="x2",color="depth_order"))+geom_point()+
    ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:20)), aes(label = name_label),size=2.5,force=50)+
    scale_color_gradientn(values = scales::rescale(c(1,round(dim(tsne_embeddings)[1]/2),dim(tsne_embeddings)[1]), to = c(0, 1)),colors = c("lightskyblue","dodgerblue","midnightblue"))+
    theme_bw()
  # regression? dimension reduction? biostatistics?
  ggsave(file=paste0("citation_networks/final10Y/group",i,".pdf"), tsne_pos1, height = 9, width = 12)
}


## change membership of individuals.
# combine clusters:
# 2,3,4,17,18 -> 1, nonparametric
# 13 -> 2, stochastic ordering
# 5 -> 3, wavelet
# 12 -> 4, functional data
# 7,8,11,19 -> 5, high dimension
# 10 -> 6, time series
# 14,21,22,23,24,25 -> 7, Bayesian
# 6 -> 8, measurement error
# 15 -> 9, computer experiment
# 26,27,28,29,30,31,32 -> 10, biostatistics
# 20 -> 11, statistical learning
# 16 -> 12, robust statistics
# 9 -> 13, multiple tests
# 1 -> 14, likelihood theory
cluster_combine <- function(clus_vec){
  clus_vec2 = clus_vec
  clus_vec[clus_vec2%in%c(2,3,4,17,18)] = 1
  clus_vec[clus_vec2==13] = 2
  clus_vec[clus_vec2==5] = 3
  clus_vec[clus_vec2==12] = 4
  clus_vec[clus_vec2%in%c(7,8,11,19)] = 5
  clus_vec[clus_vec2==10] = 6
  clus_vec[clus_vec2%in%c(14,21,22,23,24,25)] = 7
  clus_vec[clus_vec2==6] = 8
  clus_vec[clus_vec2==15] = 9
  clus_vec[clus_vec2%in%c(26,27,28,29,30,31,32)] = 10
  clus_vec[clus_vec2==20] = 11
  clus_vec[clus_vec2==16] = 12
  clus_vec[clus_vec2==9] = 13
  clus_vec[clus_vec2==1] = 14
  clus_vec
}
combined_clusname <- c("NP", "SO", "WV",
                       "FD", "HD", "TS",
                       "BY", "ME", "CE",
                       "BS", "SL", "RS",
                       "MT", "LT")
tsne_embeddings$cluster = PCA_score_angle$cluster
active_cluster = group_by(tsne_embeddings, cluster)%>%summarise(median(first_active_year))
active_cluster = c((active_cluster[,2]-1974)/2,recursive = TRUE)
names(active_cluster) = 1:32

node_Ramsay =  which(authorname =="James O. Ramsay")
name_Ramsay = authorname[node_Ramsay]
sim_Ramsay = cosine_similarity_partial(embedding_result_fASE, node_Ramsay)
cluster_Ramsay = minimal_clusters(sim_Ramsay, node_Ramsay, tsne_embeddings$cluster, active_points = active_points, active_cluster = active_cluster)
active_points_partial = active_points[node_Ramsay]
active_points_partial = lapply(active_points_partial, function(x){x[x%in%4:18]})
cluster_Ramsay[[1]][[1]] = cluster_combine(cluster_Ramsay[[1]][[1]])
cluster_Ramsay[[1]][[1]] = cluster_Ramsay[[1]][[1]][names(cluster_Ramsay[[1]][[1]])%in%4:18]
timemap = list(mappingvec = as.POSIXct(paste0(seq(1976,2015, by=2), "-01-01"), format = "%Y-%m-%d"), unit = "Year", date_labels = "%Y", limits = c(as.POSIXct("1980-01-01", format = "%Y-%m-%d"),as.POSIXct("2010-01-01", format = "%Y-%m-%d")))
membership_changeplot_Ramsay = membership_visualize(cluster_Ramsay, name_Ramsay, active_points_partial, 1:20, timemap = timemap, cluster_no = combined_clusname)
membership_changeplot_Ramsay = membership_changeplot_Ramsay + guides(color = FALSE, linetype = FALSE) + theme(strip.text = element_text(size = 24),strip.background = element_rect(fill = "white"))

node_yu = which(authorname == "Bin Yu")
name_yu = authorname[node_yu]
sim_yu = cosine_similarity_partial(embedding_result_fASE, node_yu)
cluster_yu = minimal_clusters(sim_yu, node_yu, tsne_embeddings$cluster, active_points = active_points, active_cluster = active_cluster)
active_points_partial = active_points[node_yu]
active_points_partial = lapply(active_points_partial, function(x){x[x%in%4:18]})
cluster_yu[[1]][[1]] = cluster_combine(cluster_yu[[1]][[1]])
cluster_yu[[1]][[1]] = cluster_yu[[1]][[1]][names(cluster_yu[[1]][[1]])%in%4:18]
timemap = list(mappingvec = as.POSIXct(paste0(seq(1976,2015, by=2), "-01-01"), format = "%Y-%m-%d"), unit = "Year", date_labels = "%Y", limits = c(as.POSIXct("1982-01-01", format = "%Y-%m-%d"),as.POSIXct("2010-01-01", format = "%Y-%m-%d")))
membership_changeplot_yu = membership_visualize(cluster_yu, name_yu, active_points_partial, 1:20, timemap = timemap, cluster_no = combined_clusname)
membership_changeplot_yu = membership_changeplot_yu + guides(color = FALSE, linetype = FALSE) + theme(strip.text = element_text(size = 24),strip.background = element_rect(fill = "white"))

node_larry = which(authorname == "Larry Wasserman")
name_larry = authorname[node_larry]
sim_larry = cosine_similarity_partial(embedding_result_fASE, node_larry)
cluster_larry = minimal_clusters(sim_larry, node_larry, tsne_embeddings$cluster, active_points = active_points, active_cluster = active_cluster)
active_points_partial = active_points[node_larry]
active_points_partial = lapply(active_points_partial, function(x){x[x%in%4:18]})
cluster_larry[[1]][[1]] = cluster_combine(cluster_larry[[1]][[1]])
cluster_larry[[1]][[1]] = cluster_larry[[1]][[1]][names(cluster_larry[[1]][[1]])%in%4:18]
timemap = list(mappingvec = as.POSIXct(paste0(seq(1976,2015, by=2), "-01-01"), format = "%Y-%m-%d"), unit = "Year", date_labels = "%Y", limits = c(as.POSIXct("1982-01-01", format = "%Y-%m-%d"),as.POSIXct("2010-01-01", format = "%Y-%m-%d")))
membership_changeplot_larry = membership_visualize(cluster_larry, name_larry, active_points_partial, 1:20, timemap = timemap, cluster_no = combined_clusname)
membership_changeplot_larry = membership_changeplot_larry + guides(color = FALSE, linetype = FALSE) + theme(strip.text = element_text(size = 24),strip.background = element_rect(fill = "white"))

layout <- "
AABBCC
AABBCC
"
membership_changeplot = membership_changeplot_Ramsay+membership_changeplot_yu +
  membership_changeplot_larry +
  plot_layout(design = layout) &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=26), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=28), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(membership_changeplot, file = paste0("membership_change_all.pdf"), width = 27, height = 9)

