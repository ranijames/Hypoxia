#!/usr/bin/env Rscript
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
# The script can be executed as: Rscript pathway_score.R ../data/sample_input.tsv "cancertype of interest" " sample of interest"
# The gene expression matrix is called #
# The duplicated samples or removed here
args             = commandArgs(TRUE)
Gene_expression  = read.table(args[1], sep="\t",header=TRUE,check.names=FALSE)
#print(string)
cancertype       = as.character(args[2])
string           = as.character(args[3])

#print(cancertype)
# Output directory
results_dir      = "/Users/alvajames/Mount_projects/Results/"
expected_output  = paste(results_dir, "sample_output.tsv", sep="")
outfile          = paste(results_dir, "sample_pathway_bayes_core.pdf", sep="")

# The gene weight is called here #
Gene_weights               = read.table("/Users/alvajames/Mount_projects/data/gene_weights.tsv", header=TRUE, sep="\t",check.names=FALSE)

# The pathway score function #

Get_pathwayscore<- function(Gene_weights, Gene_expression){
  
  if(sum(colnames(Gene_expression) %in% c("Y", "sample_id")) != 2){
    print("Need column names Y and sample_id")
    return(NA)
  }
  
  if(sum(colnames(Gene_weights) %in% c("gene_weight", "gene_ids")) != 2){
    print("Need column names gene_weight and gene_ids")
    return(NA)
  }
  
  Gene_expression            = Gene_expression[!duplicated(Gene_expression$sample_id), ]
  gene_ids = setdiff(colnames(Gene_expression), c("Y", "sample_id"))
  
  # The gene expression matrix normalization
  Gene_expression_data       = Gene_expression[,gene_ids]
  Gene_expression_data       = as.matrix(t(scale(log10(t(Gene_expression_data+1)))))
  Gene_expression_data[is.nan(Gene_expression_data)] <- 0
  row.names(Gene_expression_data) = Gene_expression$sample_id
  
  # format the projection vector
  Gene_weigh_vector              = Gene_weights$gene_weight
  names(Gene_weigh_vector)       = Gene_weights$gene_ids
  Gene_weigh_vector              = as.vector(t((Gene_weights$gene_weight)))
  
  
  if(sum(names(Gene_weigh_vector) != colnames(Gene_expression_data)) > 0){
    print("Error Matching Gene ids")
    return(NA)
  }
  
  # get projection
  Pathway_score = Gene_expression_data %*% Gene_weigh_vector
  Pathway_score = data.frame(pathway_score=Pathway_score, sample_id=row.names(Pathway_score))
  Pathway_score = merge(Gene_expression[,c("sample_id", "Y")], Pathway_score)
  Pathway_score = Pathway_score[order(Pathway_score$pathway_score),]
  
  return(Pathway_score)
  
}
Pathway_score = Get_pathwayscore(Gene_weights, Gene_expression)
write.table(Pathway_score,file = expected_output, sep="\t",row.names=F)

# For SKCM cancer type #
meta_data_gtex = read.table("/Users/alvajames/Mount_projects/data/gtex_metadata.tsv", header=TRUE, sep="\t",check.names=FALSE)
meta_data_tcga = read.table("/Users/alvajames/Mount_projects/data/tcga_metadata.tsv", header=TRUE, sep="\t",check.names=FALSE)
Metadata       = merge(as.data.frame.matrix(meta_data_tcga), as.data.frame.matrix(meta_data_gtex), all = TRUE, sort = FALSE)
sample_id_SCKM = Metadata$sample_id[Metadata$study %in% cancertype]
SCKM           = Gene_expression[Gene_expression$sample_id %in% sample_id_SCKM,]

# Get the normal (0) and cancer  (1) backgrounds for SKCM cancer type from the pathway scores,
Pathway_score_SCKM = merge(SCKM[,1:2],Pathway_score,on="sample_id")

  Cancer = Pathway_score_SCKM[Pathway_score_SCKM$Y == 1, ][3]
  Normal = Pathway_score_SCKM[Pathway_score_SCKM$Y == 0, ][3]

   #Finding the Bayes score using the above information
  classify_bayes_score <- function(Normal, Cancer, Pathway_score_SCKM){
    
    # calculate odds ratio
    mean_1 = mean(as.matrix(Cancer))
    sd_1   = sd(as.matrix(Cancer))
    
    mean_0 = mean(as.matrix(Normal))
    sd_0   = sd(as.matrix(Normal))
    
    liklihood_1 = dnorm(Pathway_score_SCKM$pathway_score, mean_1, sd_1)
    liklihood_0 = dnorm(Pathway_score_SCKM$pathway_score, mean_0, sd_0)
    
    bayes_factor           = liklihood_1/liklihood_0
    bayes_factor           = as.matrix(bayes_factor)
    classifications        = 2*log(bayes_factor)
    colnames(classifications) = c("bayes_factor")
    classification_new        = cbind(classifications, Pathway_score_SCKM)
    classification_new        = classification_new[c("sample_id", "bayes_factor", "pathway_score","Y")]
    return(classification_new)
    
  }
  
  Bayes_score = classify_bayes_score(Normal,Cancer,Pathway_score_SCKM)
  write.table(Bayes_score,file = expected_output,sep="\t",row.names = FALSE)
  
  plotter <- function(Cancer, Normal, Bayes_score, string){
    if (nrow(Bayes_score[which(Bayes_score$sample_id==string),])!=0) {
      mylabel      = paste('Bayes factor =',Bayes_score[which(Bayes_score$sample_id==string),]$bayes_factor)
      pathlabel    = paste('Pathway score =',Bayes_score[which(Bayes_score$sample_id==string),]$pathway_score)
      sample_label = paste('Sample =',string)
      p1 = ggplot(data = Normal, aes(x=pathway_score,fill = "blue")) + 
        geom_density(alpha=0.6) + 
        geom_density(data = Cancer, aes(x=pathway_score,fill = "green"),alpha=0.6) + 
        scale_fill_manual(labels = c("bg0", "bg1"), values = c("blue", "#f45c42")) +
        geom_vline(xintercept = Bayes_score[which(Bayes_score$sample_id==string),]$bayes_factor, 
                   linetype="dotted", color = "black", size=0.5)+
        geom_text(aes(x=Bayes_score[which(Bayes_score$sample_id==string),]$bayes_factor,
                      label=paste(mylabel,sep="\n",
                                  pathlabel, sample_label),y=1,hjust = 1.2), colour="black",hjust = -0.12,vjust=0.8) + 
        ggtitle("Plot of pathway score distribution")
      
      #################### histogram ###############
      p2 = ggplot(data = Normal, aes(x=pathway_score,fill = "blue")) + 
        geom_histogram(alpha=0.6,binwidth =0.2,bins=70) + 
        geom_histogram(data = Cancer, aes(x=pathway_score,fill = "green"),alpha=0.6,binwidth =0.2,bins=70) + 
        scale_fill_manual(labels = c("bg0", "bg1"), values = c("blue", "#f45c42")) +
        geom_vline(xintercept = Bayes_score[which(Bayes_score$sample_id==string),]$bayes_factor, 
                   linetype="dotted", color = "black", size=0.5) +
        geom_text(aes(x=Bayes_score[which(Bayes_score$sample_id==string),]$bayes_factor, 
                      label=paste(mylabel,sep="\n",pathlabel, sample_label), y = max(Normal$pathway_score), hjust = 1, vjust = 20),
                  colour="black",hjust = -0.025,vjust=-10.1) 
      pdf(outfile)
      print(p1)
      print(p2)
      dev.off() 
    }
  }

  plotter(Cancer, Normal, Bayes_score, string)
  

