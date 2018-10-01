library(ggplot2)
library(plotly)
library(plyr)
args        = commandArgs(TRUE)
results_dir = "/Users/alvajames/Mount_projects/Results/"
#matrix      = read.table(args[3], sep="\t",header=TRUE,check.names=FALSE)
sample_of_interest =args[4]
outfile            =paste(results_dir, "test.pdf", sep="")

############################################################
############ The plotter function ##########################
############################################################

plotter <- function(bg0, bg1, dataset1, string){
  if (nrow(dataset1[which(dataset1$sample_id==string),])!=0) {
    mylabel      = paste('Bayes factor =',dataset1[which(dataset1$sample_id==string),]$bayes_factor)
    pathlabel    = paste('Pathway score =',dataset1[which(dataset1$sample_id==string),]$pathway_score)
    sample_label = paste('Sample =',string)
    p1 = ggplot(data = bg0, aes(x=pathway_score,fill = "blue")) + 
      geom_density(alpha=0.6) + 
      geom_density(data = bg1, aes(x=pathway_score,fill = "green"),alpha=0.6) + 
      scale_fill_manual(labels = c("bg0", "bg1"), values = c("blue", "#f45c42")) +
      geom_vline(xintercept = dataset1[which(dataset1$sample_id==string),]$bayes_factor, 
                   linetype="dotted", color = "black", size=0.5)+
      geom_text(aes(x=dataset1[which(dataset1$sample_id==string),]$bayes_factor,
                    label=paste(mylabel,sep="\n",
      pathlabel, sample_label),y=1,hjust = 1.2), colour="black",hjust = -0.12,vjust=0.8) + 
      ggtitle("Plot of pathway score distribution")
    
    #################### histogram ###############
    p2 = ggplot(data = bg0, aes(x=pathway_score,fill = "blue")) + 
      geom_histogram(alpha=0.6,binwidth =0.2,bins=70) + 
      geom_histogram(data = bg1, aes(x=pathway_score,fill = "green"),alpha=0.6,binwidth =0.2,bins=70) + 
      scale_fill_manual(labels = c("bg0", "bg1"), values = c("blue", "#f45c42")) +
      geom_vline(xintercept = dataset1[which(dataset1$sample_id==string),]$bayes_factor, 
                 linetype="dotted", color = "black", size=0.5) +
      geom_text(aes(x=dataset1[which(dataset1$sample_id==string),]$bayes_factor, 
                    label=paste(mylabel,sep="\n",pathlabel, sample_label), y = max(bg0$pathway_score), hjust = 1, vjust = 20),
                colour="black",hjust = -0.025,vjust=-10.1) 
    pdf(outfile)
    print(p1)
    print(p2)
    dev.off() 
  }
}

  