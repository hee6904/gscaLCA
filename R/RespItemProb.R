RespItemProb <- function(dat, varnames, mem.dat, LEVELs)
{
  dat = dat[,varnames] 

  mem.dat <- data.frame(mem.dat)
  dat$label <- mem.dat$label

  mem.mat.total <- list()

  for (i in 1:(ncol(dat)-1)) # i: variable
  {

    mem.mat <- table(dat$label, dat[,i])/rowSums(table(dat$label, dat[,i]))
    dimnames(mem.mat)[[2]] <-  LEVELs[[i]]
    Categories = expand.grid(dimnames(mem.mat),stringsAsFactors=F )

    mat = cbind(Categories[order(Categories$Var1),], as.vector(t(mem.mat)))

    names(mat) = c("Class", "Category", "Estimate")
    rownames(mat) = 1:nrow(mat)

    mem.mat.total[[i]] <- mat
  }

  names(mem.mat.total) <- colnames(dat)[1:(ncol(dat)-1)]
  return(mem.mat.total)

}
