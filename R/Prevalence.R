Prevalence  = function(U.0, c, nobs)
{
  classfied = matrix(0, nrow=c, ncol=2)
  rownames(classfied)<- paste("Latent Class",1:c)
  colnames(classfied) <- c("Percent", "Count")

  nclass <- table(U.0$label)
  classfied[names(nclass),2] = nclass
  classfied[ ,1] = classfied[ ,2]/nobs*100

  classfied
}
