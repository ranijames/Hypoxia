# get projection
#mat_scaled on Z-score formulae; for normalization optional
args <- commandArgs(TRUE)
dca_matrix <- read.table(args[1], sep="\t",header=TRUE, sep="\t",check.names=FALSE)
proj_vector<- read.table(args[2], sep="\t",header=TRUE, sep="\t",check.names=FALSE)
dca_matrix  = read.table("sample_input-2.tsv", header=TRUE, sep="\t",check.names=FALSE)
proj_vector  =read.table("gene_weights_hypoxia-2.tsv", header=TRUE, sep="\t",check.names=FALSE)
Output =read.table("sample_output-2.tsv", header=TRUE, sep="\t",check.names=FALSE)
dca_matrix  = dca_matrix[!duplicated(dca_matrix$sample_id), ]
rownames(dca_matrix)=dca_matrix[,1]
dca_matrix$sample_id<-NULL
#dca_data    =t(apply(log10(dca_data), 1, function(x) (x - mean(x)) / sd(x))) ## I used it for now
dca_data     =as.matrix(t(scale(log10(t(dca_data +1)))))
dca_data[is.nan(dca_data)]<-0
vector_score=as.vector(t((proj_vector$gene_weight)))
matches <- match(proj_vector$gene_id, colnames(dca_data))
dca_data <- dca_data[,matches]
data_proj       =dca_data %*% vector_score
data_proj =as.data.frame(data_proj)
setDT(data_proj, keep.rownames = TRUE)[]
colnames(data_proj)=c("sample_id","hypoxia_score_new")
data_proj=merge(data_proj,Output, on='sample_id')
data_proj=merge(data_proj,Output, on='sample_id')
data_proj$Diff_hypoxia= data_proj$hypoxia_score_new- data_proj$hypoxia_score



