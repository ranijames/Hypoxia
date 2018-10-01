library(ggplot2)
library(plotly)
library(plyr)
args         = commandArgs(TRUE)
bayer       = read.table("/Users/alvajames/Mount_projects/data/sample_bayes.tsv", header=TRUE, sep="\t")
gtex_bg0    = read.table("/Users/alvajames/Mount_projects/data/gtex_metadata.tsv", header=TRUE, sep="\t")
bayes_factor= read.table("/Users/alvajames/Mount_projects/data/Bayes_score.tsv", header=TRUE)
colnames(bayes_factor)[1] = "sample_id"
tcga_bg1    = read.table("/Users/alvajames/Mount_projects/data/tcga_metadata.tsv", header=TRUE, sep="\t")
bg0         = read.table("/Users/alvajames/Mount_projects/data/gtex_background.tsv", header=TRUE)
bg1         = read.table("/Users/alvajames/Mount_projects/data/tcga_background.tsv", header=TRUE)
colnames(bg0)=c("pathway_score")
colnames(bg1)=c("pathway_score")
bg0         = merge(bg0,bayes_factor,on="pathway_score")
bg1         = merge(bg1,bayes_factor,on="pathway_score")
bg1         = bg1[bg1$Y == "1",]
bg1         = bg1[1:3]
bg0         = bg0[1:3]
bg0_gtex_meta = merge(bg0,unique(gtex_bg0[c("sample_id","study")]),on="sample_id") 
bg1_tcga_meta = merge(bg1,unique(tcga_bg1[c("sample_id","study")]),on="sample_id") 
bg1_tcga_meta$type = "TCGA"
bg0_gtex_meta$type = "GTEX"
gtex_tcga          = rbind(bg1_tcga_meta,bg0_gtex_meta)
gtex_tcga_mean     = ddply(gtex_tcga, "type", summarise, rating.mean = mean(bayes_factor))

############################################################
############ The plotter function ##########################
###########################################################

plotter <- function(bg1, bg2, matrix, sample_of_interest){
  if (sample_of_interest = matrix$sample_id) 
  ggplot(data = gtex_tcga, aes(x=pathway_score,fill = type)) + geom_density (alpha = .3)  + geom_label(babel=sprintf('n = %d',sample_of_interest,gtex_tcga$sample_id))
    geom_vline(data = gtex_tcga_mean, aes(xintercept = rating.mean, colour = type),linetype = "longdash", size=1)
  
  
}
  
p1  = ggplot(data = gtex_tcga, aes(x=pathway_score,fill = type)) + geom_density (alpha = .3)  + 
  geom_label(aes(y=gtex_tcga$pathway_score,label= gtex_tcga$sample_id), size = 3.5)
  
  geom_vline(data = gtex_tcga_mean, aes(xintercept = rating.mean, colour = type),linetype = "longdash", size=1)
p1 + annotate('text', x = gtex_tcga$sample_id, y = gtex_tcga$bayes_factor,label = sprintf('n = %s',gtex_tcga$sample_id), vjust = 0)
ggplotly(p1)
