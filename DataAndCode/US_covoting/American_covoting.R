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
library(ClusterR)
library(aricode)
library(mrfDepth)
library(Rtsne)
library(patchwork)
library(ggalt)

############  You may skip this procedure of data preparation.
## Due to the upload limit of Github, you can download the raw files in: https://voteview.com/data
HSall_votes = read_csv("HSall_votes.csv")
Sall_members = read_csv("Sall_members.csv")

HSall_data = dplyr::filter(HSall_votes, congress>=40&(chamber=="Senate"|chamber=="President"))

American_covoting = list()
for (congress_now in 40:117){
  covoting_thisyear = c()
  voting_thisyear = filter(HSall_data, congress == congress_now)
  active_id = unique(voting_thisyear$icpsr)

  for(i in active_id){
    for(j in active_id){
      if(i!=j){
        i_voting_thisyear = filter(voting_thisyear, icpsr == i)%>%dplyr::select(rollnumber, cast_code)%>%rename(cast_code1 = 2)
        j_voting_thisyear = filter(voting_thisyear, icpsr == j)%>%dplyr::select(rollnumber, cast_code)%>%rename(cast_code2 = 2)

        ij_voting_thisyear = inner_join(i_voting_thisyear, j_voting_thisyear)

        ij_voting_thisyear = ij_voting_thisyear%>%mutate(covoting = (cast_code1 * cast_code2 == 1|cast_code1 * cast_code2 == 36))
        if(!is.null(ij_voting_thisyear)&&dim(ij_voting_thisyear)[1]!=0){
          if(sum(ij_voting_thisyear$covoting)/length(ij_voting_thisyear$covoting)>=0.75){
            covoting_thisyear = rbind(covoting_thisyear, t(c(i,j)))
          }
        }
      }
      cat(congress_now, "\n")
    }
  }
  American_covoting = c(American_covoting, list(list(covoting_thisyear)))
}
index = pmatch(id, Sall_members$icpsr)
party = c(Sall_members[index, "party_code"], recursive=TRUE)
save(Sall_members, American_covoting, id, party, file="American_covoting.RData")
########################

source("DataAndCode\\functional_ASE.R")
source("DataAndCode\\downstream_tasks.R")
load("American_covoting.RData")

nodes = id
n_nodes = length(nodes)
m_layers = length(79:117)
dynamic_network_adjacency = array(0, dim = c(n_nodes, n_nodes, m_layers), dimnames = list(nodes,nodes,79:117))
for(congress_now in 79:117){
  if(!is.null(American_covoting[[congress_now-39]][[1]])){
    list_now = (American_covoting[[congress_now-39]][[1]])
    list_now = list_now[list_now[,1]<list_now[,2],]
    if(is.null(dim(list_now)[1])){
      list_now = t(list_now)
    }


    for(d in 1:dim(list_now)[1]){
      list_now[d,1] = which(list_now[d,1]==id)
      list_now[d,2] = which(list_now[d,2]==id)
    }
    i = congress_now - 78
    max_active = max(list_now)

    network_slice_graph = graph.edgelist(list_now, directed = FALSE)
    network_now = as.matrix( get.adjacency(network_slice_graph))

    dynamic_network_adjacency[1:max_active,1:max_active,i] = network_now
  }
}
dynamic_network_adjacency_spare = dynamic_network_adjacency

id_out = c()
h = 0
while(length(id_out)==0){
  h = h+1
  for(j in 1:n_nodes){
    if(all(dynamic_network_adjacency[j,,]==0)&&all(dynamic_network_adjacency[,j,]==0)){
      id_out = c(id_out, j)
    }else{
      degs = apply(dynamic_network_adjacency[j,,],2,sum)
      if(length(which(degs!=0))<3){
        id_out = c(id_out, j)
      }
    }
  }
  if(length(id_out)==0) break
  nodes = nodes[-id_out]
  party = party[-id_out]
  dynamic_network_adjacency = dynamic_network_adjacency[-id_out,-id_out,]
  n_nodes = length(nodes)
  if(length(id_out)!=0) id_out = c()
  else id_out = TRUE
}

party2 = party
party2[party==100] = "Democratic"
party2[party==200] = "Republican"
party2[party==328] = "Independent"
party = party2

m_layers = 117-78

############  You may skip this procedure of embedding.
tuning_sim_US = fASE_tuning(10:15, 5:10, tuning_method = "heuristic", adjacency_tensor = dynamic_network_adjacency, spline_order = 4, batch_size = NULL, kernel_scale = FALSE , timestamp_vec = NULL, scalable_method = TRUE, scalable_dim = 30, scalable_power = 6, iteration_method = "sequential", iteration_step = 2500, epsilon=1e-6)
tuning_parameter_US = get_tuning(tuning_sim_US)

t1 = Sys.time()
result = fASE(dynamic_network_adjacency, tuning_parameter_US[[2]], tuning_parameter_US[[1]], batch_size = NULL, timestamp_vec = 79:117, scalable_power = 10, scalable_dim = 30, epsilon = 1e-6, iteration_step = 2500)
t2 = Sys.time()
embedding_result_fASE = eval.fd(79:117, result[[1]])

save(file = "US_covoting.RData", tuning_sim_US, result, embedding_result_fASE, dynamic_network_adjacency, nodes, party)
#######################

load("US_covoting.RData")
n_nodes = dim(dynamic_network_adjacency)[2]
timestamp_vec = 79:117
active_points = lapply(1:n_nodes, active_calculation, dynamic_network_adjacency, timestamp_vec = timestamp_vec)
active_time = sapply(active_points, function(x){mean(x)}, simplify=TRUE)
active_duration = sapply(active_points, function(x){length(x)},simplify=TRUE)

## generalized DCCD plot.
set.seed(1)
PCA_angle = fPCA_angle(embedding_result_fASE, 79:117, 200, active_points)
which(cumsum(PCA_angle[[3]])>0.8)
s = cluster_no_selection(PCA_angle[[1]][,1:32], 10,cluster_method = "tk_means", K_selection = "Silhouette")
s[[2]]$cluster = pmatch(s[[2]]$cluster, c(2,5,1,3,6,4,7),duplicates.ok = TRUE)
PCA_angle = fPCA_angle(embedding_result_fASE, 79:117, 2, active_points)
PCA_score_angle=PCA_angle[[1]]
PCA_score_angle = as.data.frame(PCA_score_angle)
names(PCA_score_angle) = c("theta1", "theta2")
PCA_score_angle$label = factor(party, levels = c("Democratic", "Republican", "Independent"))
PCA_score_angle$cluster = factor(s[[2]]$cluster)
PCA_score_angle$mean_time = active_time
PCA_score_angle$active_duration = active_duration
angle_plot1 = ggplot(PCA_score_angle,aes(x=theta1,y=theta2,group=cluster,shape=label,color=cluster,fill=cluster))+guides(fill="none")+geom_point(size=3)+scale_color_manual(values=c("#CFEAF1","#82B0D2","#386CB0","#FDC086","#F0988C","#C82423","#E7DAD2"))+scale_fill_manual(values=c("#CFEAF1","#82B0D2","#386CB0","#FDC086","#F0988C","#C82423","#E7DAD2"))+geom_encircle(aes(group = cluster,fill=cluster),expand=0.05,spread=0.5,s_shape=1,size=3,linetype = 1,alpha=0.2)+theme_bw()
angle_plot2 = ggplot(PCA_score_angle,aes(x=theta1,y=theta2,group=cluster,shape=label,color=mean_time,size=active_duration))+geom_point()+scale_color_gradientn(values = scales::rescale(c(80,100,110,120), to = c(0, 1)), colors = c("#CCE5FF","lightskyblue","#5383c3","midnightblue"),breaks = c(90,100,110))+scale_size_continuous(range = c(2.5,4.5))+guides(shape=guide_legend(override.aes = list(size=c(3,3,3))))+theme_bw()

layout <- "ABCDE"
fPCA_score_plot = (plot_spacer()+angle_plot1+xlab(TeX("$x_1^*$"))+ylab(TeX("$x_2^*$"))+labs(caption =  "(a) DCCD plot with estimated\n cluster for each node.")+plot_spacer()+angle_plot2 +xlab(TeX("$x_1^*$"))+ylab(TeX("$x_2^*$"))+ labs(caption =  "(b) DCCD plot with active\n duration for each node.") +plot_spacer())+
   plot_layout(design=layout, widths = c(4.5,14,4,14,4.5)) &
   theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))

ggsave(fPCA_score_plot, file = "fPCA_score.pdf", width = 27, height = 7.8)

## radius-corrected functions.
node_embedding_frame = function_sort_angle(embedding_result_fASE, 79:117, active_points)
node_embedding_frame = cbind(node_embedding_frame, cluster = rep(factor(s[[2]]$cluster), each = 117-79+1))
node_embedding_frame = cbind(node_embedding_frame, party = rep(factor(party, levels = c("Democratic", "Republican", "Independent")), each = 117-79+1))
func_plot_best_x1_angle = ggplot(node_embedding_frame,aes_string(x="t",y=paste0("x",1),group="i",color="cluster"))+facet_wrap(~cluster)+geom_line(aes(linetype=party))+scale_color_manual(values=c("#CFEAF1","#82B0D2","#386CB0","#FDC086","#F0988C","#C82423","#E7DAD2"))+scale_size_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+scale_alpha_manual(values= rep(1,times=5))+scale_linetype_manual(values = c("solid","dotted","longdash"))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()+theme(strip.text = element_text(size = 18),strip.background = element_rect(fill = "white"),plot.margin = unit(c(1, 9, 1, 9), "lines"))
func_plot_best_x2_angle = ggplot(node_embedding_frame,aes_string(x="t",y=paste0("x",2),group="i",color="cluster"))+facet_wrap(~cluster)+geom_line(aes(linetype=party))+scale_color_manual(values=c("#CFEAF1","#82B0D2","#386CB0","#FDC086","#F0988C","#C82423","#E7DAD2"))+scale_size_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+scale_alpha_manual(values= rep(1,times=5))+scale_linetype_manual(values = c("solid","dotted","longdash"))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()+theme(strip.text = element_text(size = 18),strip.background = element_rect(fill = "white"),plot.margin = unit(c(1, 9, 1, 9), "lines"))
func_plot_best_x3_angle = ggplot(node_embedding_frame,aes_string(x="t",y=paste0("x",3),group="i",color="cluster"))+facet_wrap(~cluster)+geom_line(aes(linetype=party))+scale_color_manual(values=c("#CFEAF1","#82B0D2","#386CB0","#FDC086","#F0988C","#C82423","#E7DAD2"))+scale_size_manual(values= c(0.8,0.8,0.8,0.8,0.8,0.8,0.8))+scale_alpha_manual(values= rep(1,times=5))+scale_linetype_manual(values = c("solid","dotted","longdash"))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()+theme(strip.text = element_text(size = 18),strip.background = element_rect(fill = "white"),plot.margin = unit(c(1, 9, 1, 9), "lines"))

function_plot_best_angle = (func_plot_best_x1_angle +ylab(TeX("$Z_1^*$"))+guides(color="none", linetype="none") +labs(caption = "(a) the first dimension of the (radius corrected) multivariate functional data.") )/
       (func_plot_best_x2_angle +ylab(TeX("$Z_2^*$"))+labs(caption = "(b) the second dimension of the (radius corrected) multivariate functional data.") )+
  plot_layout(guides = "collect") & theme(axis.text=element_text(size=18), axis.title =element_text(size=20), plot.caption =element_text(size=24, hjust=0.5), legend.text=element_text(size=22), legend.title =element_text(size=22), legend.background = element_rect(fill = 'white', colour = 'black'))

ggsave(function_plot_best_angle,file="latent_functions.pdf",width=18,height=22)

# DD plot, 79-117th
set.seed(1)
chamber = c(Sall_members[pmatch(nodes,Sall_members$icpsr),"chamber"],recursive=TRUE)

sim=cosine_similarity(eval.fd(seq(79,117,length=100), result[[1]]), seq(79,117,length=100))
depth_stats = metric_depth(1-sim)
depth_reindex = order(depth_stats,decreasing =TRUE)
tsne_result <- Rtsne(1-sim, dims = 2, perplexity = 30, is_distance=TRUE)
tsne_embeddings <- tsne_result$Y
tsne_embeddings=as.data.frame(tsne_embeddings)
colnames(tsne_embeddings) = c("x1","x2")
tsne_embeddings$label=party
tsne_embeddings$label[tsne_embeddings$chamber=="President"]="President"
tsne_embeddings$label =  factor(tsne_embeddings$label, levels = c("Democratic", "Republican", "Independent", "President"))
tsne_embeddings[depth_reindex,"depth_order"] = as.numeric(1:dim(tsne_embeddings)[1])
tsne_embeddings$name = Sall_members$bioname[pmatch(nodes, Sall_members$icpsr)]
tsne_embeddings$name_label = paste0(tsne_embeddings$name,",",tsne_embeddings$depth_order)
tsne_embeddings$depth = depth_stats
tsne_embeddings$degree = apply(dynamic_network_adjacency,1,sum)

tsne_pos1 = ggplot(tsne_embeddings,aes_string(x="depth",y="degree",color="label"))+geom_point( aes(size=label,alpha=label,shape=label))+annotation_custom(
  grob = grid::linesGrob(y = unit(c(1, 0), "npc")), xmin = -Inf, ymin =Inf, xmax = Inf, ymax = -Inf)+
  ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:15)), aes(label = name_label),size=5,force=80,max.iter=1000000,show.legend=FALSE)+
  scale_color_manual(values = c("#386CB0","#C82423","gray","black"))+
  scale_size_manual(values=c(3,3,3,3,3,3,3,3))+
  scale_alpha_manual(values= c(1,1,1,1,1,1,1,1,1,1))+guides(color=guide_legend(override.aes = list(size=7)))+
  theme_bw()+theme(axis.text=element_text(size=28), axis.title =element_text(size=24),  plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))

## DD plot, Republican & 110-117th
set.seed(1)
active_republicans = intersect(which(sapply(active_points, condition <- function(x) {
  any(x%in%(110:117))
})),which(party=="Republican"))
sim=cosine_similarity(eval.fd(seq(110,117,length=100), result[[1]][active_republicans]), seq(110,117,length=100))
depth_stats = metric_depth(1-sim)
depth_reindex = order(depth_stats,decreasing =TRUE)
tsne_result <- Rtsne(1-sim, dims = 2, perplexity = 10, is_distance=TRUE)
tsne_embeddings <- tsne_result$Y
tsne_embeddings=as.data.frame(tsne_embeddings)
colnames(tsne_embeddings) = c("x1","x2")
rownames(tsne_embeddings) = active_republicans
tsne_embeddings$label = "Republicans"
tsne_embeddings$label = factor(tsne_embeddings$label)
tsne_embeddings$name = c(Sall_members[pmatch(nodes[active_republicans], Sall_members$icpsr),"bioname"],recursive=TRUE)
tsne_embeddings[depth_reindex,"depth_order"] = as.numeric(1:dim(tsne_embeddings)[1])
tsne_embeddings$name_label = paste0(c(tsne_embeddings$name,recursive=TRUE),",",tsne_embeddings$depth_order)
tsne_embeddings$degree = apply(dynamic_network_adjacency[active_republicans,active_republicans,32:39],1,sum)
tsne_embeddings$depth = depth_stats
tsne_pos2 = ggplot(tsne_embeddings,aes_string(x="depth",y="degree",color="label"))+geom_point( aes(size=label,alpha=label),shape=17)+annotation_custom(
  grob = grid::linesGrob(y = unit(c(1, 0), "npc")), xmin = -Inf, ymin =Inf, xmax = Inf, ymax = -Inf)+
  ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(1:10)), aes(label = name_label),size=5,force=50,max.iter=1000000,max.overlaps=15,show.legend = FALSE)+
  ggrepel::geom_label_repel(data = subset(tsne_embeddings, depth_order%in%c(57,74,63,72)), aes(label = name),size=5,force=80,fontface="bold",show.legend = FALSE)+
  scale_color_manual(values = c("#C82423","gray","black"))+
  scale_size_manual(values=c(3,3,3,3,3,3,3,3))+
  scale_alpha_manual(values= c(1,1,1,1,1,1,1,1,1,1))+guides(color=guide_legend(override.aes = list(size=7)))+
  theme_bw()+theme(axis.text=element_text(size=28), axis.title =element_text(size=24),  plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))

layout <- "EAAFBBG"
tsne_pos = tsne_pos1+ labs(caption =  "(a) DD plot for all senators.")+tsne_pos2+guides(color="none",size="none",shape="none",alpha="none")+ labs(caption =  "(b) DD plot for Republican\n senators from 110 to 117th Congress.")+
  plot_spacer()+plot_spacer()+plot_spacer()+
  plot_layout(guides = "collect", design = layout, heights = 8,width = c(4.5,7,7,6,7,7,4.5)) &
  theme(axis.text=element_text(size=24), axis.title =element_text(size=24), plot.caption =element_text(size=30, hjust=0.5), legend.text=element_text(size=28), legend.title =element_text(size=30), legend.background = element_rect(fill = 'white', colour = 'black'))
ggsave(tsne_pos, file = "centralityofsenators.pdf", height = 7.5, width = 24)

