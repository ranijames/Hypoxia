#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
args         = commandArgs(TRUE)
background1  = read.table(args[1], sep="\t",header=TRUE,check.names=FALSE)
background1  = as.numeric(unlist(background1))
background2  = read.table(args[2], sep="\t",header=TRUE,check.names=FALSE)
background2  = as.numeric(unlist(background2))
score        = read.table(args[3],sep='\t',header=TRUE,check.names = FALSE)
dca_matrix   = read.table("/Users/alvajames/Mount_projects/data/sample_input.tsv", header=TRUE, sep="\t",check.names=FALSE)
proj_vector  = read.table("/Users/alvajames/Mount_projects/data/gene_weights.tsv", header=TRUE, sep="\t",check.names=FALSE)
Output       = read.table("/Users/alvajames/Mount_projects/data/sample_scores.tsv", header=TRUE, sep="\t",check.names=FALSE)
dca_matrix            = dca_matrix[!duplicated(dca_matrix$sample_id), ]
rownames(dca_matrix)  = dca_matrix[,1]
dca_matrix$sample_id  = NULL
dca_data              = as.matrix(t(scale(log10(t(dca_matrix[2:35] +1)))))
dca_data[is.nan(dca_data)] = 0
  ## ordering the column names based on the gene_ids in vector_proj
  matches                    = order(match(colnames(dca_data),proj_vector$gene_ids,))
  dca_data                   = dca_data[,matches]
  dim(dca_data)
  vector_score               = as.vector(t((proj_vector$gene_weight)))
  data_proj                  = dca_data[,1:34] %*% vector_score
  data_proj                  = as.data.frame(data_proj)
  setDT(data_proj, keep.rownames= TRUE)[]
  colnames(data_proj)           = c("sample_id","hypoxia_score_new")
  data_proj                     = merge(data_proj,Output,on='sample_id')
  data_proj$Diff_hypoxia        = (data_proj$hypoxia_score_new - data_proj$pathway_score)
  sc                            = data_proj
  sc                            = sc[!duplicated(sc$sample_id), ]
  write.table(sc,file = "/Users/alvajames/Mount_projects/data/Hypoxia_score", sep="\t",row.names=F)
  ####################################
  ####### For bayes score ############
  ####################################
  bg0  = read.table("/Users/alvajames/Mount_projects/data/gtex_background.tsv", header=TRUE)
  bg1  = read.table("/Users/alvajames/Mount_projects/data/tcga_background.tsv", header=TRUE)
  sc   = read.table("/Users/alvajames/Mount_projects/data/Hypoxia_score", header=TRUE, sep="\t",check.names=FALSE,row.names=1)
  bayer= read.table("/Users/alvajames/Mount_projects/data/sample_bayes.tsv", header=TRUE, sep="\t",row.names=1)
  score_var=sc[1]
  # calculate odds ratio
  mean_0 = mean(as.matrix(bg0))
  sd_0   = sd(as.matrix(bg0))
  mean_1 = mean(as.matrix(bg1))
  sd_1   = sd(as.matrix(bg1))
  liklihood_1               = dnorm(score_var, mean_1, sd_1)
  liklihood_0               = dnorm(score_var, mean_0, sd_0)
  bayes_factor              = liklihood_1/liklihood_0
  classifications           = 2*log(bayes_factor)
  colnames(classifications) = c("bayes_score_new")
  head(classifications)
  classification_new        = merge(bayer,classifications,by="row.names")
  colnames(classification_new)=c("sampl_id","bayes_factor","bayes_factor_new")
  setDT(sc, keep.rownames= TRUE)[]
  colnames(sc)              =c("sampl_id","pathway_score_new","Y","pathway_score","Difference")
  classification_new_pathway= merge(classification_new,sc,on="sample_id")
  classification_new_pathway$bayes_factor_diff =   classification_new_pathway$bayes_factor -classification_new_pathway$bayes_factor_new
  write.table(classification_new_pathway,file = "/Users/alvajames/Mount_projects/data/Bayes_score.tsv", sep="\t",row.names=F)
  
