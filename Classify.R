args        = commandArgs(TRUE)
background1 = read.table(args[1], sep="\t",header=TRUE,check.names=FALSE)
background2 = read.table(args[2], sep="\t",header=TRUE,check.names=FALSE)
score       = read.table(args[3],sep='\t',header=TRUE,check.names = FALSE)
bg1  = read.table("sample_input-2.tsv", header=TRUE, sep="\t",check.names=FALSE)
bg2  = read.table("sample_input-2.tsv", header=TRUE, sep="\t",check.names=FALSE)
sc  = read.table("sample_input-2.tsv", header=TRUE, sep="\t",check.names=FALSE)

  # calculate odds ratio
  mean_1 = mean(background_1)
  sd_1   = sd(background_1)

  mean_0 = mean(background_0)
  sd_0   = sd(background_0)

  liklihood_1 = dnorm(scores_to_eval, mean_1, sd_1)
  liklihood_0 = dnorm(scores_to_eval, mean_0, sd_0)

  bayes_factor           = liklihood_1/liklihood_0
  classifications        = 2*log(bayes_factor)
  names(classifications) = nam
