test_binomial = function(dat.cov, COVNAMES, membership.1 , num.cluster, partition, BC.method = "BHC")
{
  U.cov = data.frame(cbind(membership.1, dat.cov[,COVNAMES]))
  # nCOV = length(COVNAMES)
  # names(U.cov)[(ncol(membership.1)+1):ncol(U.cov)]=

  names(U.cov)[(ncol(membership.1)+1):ncol(U.cov)] = COVNAMES

  output = list()
  output.coff = list()
  for(k in 1:num.cluster)
  {
    U.cov$ClusM.k = ifelse(U.cov$label==paste0("Latent Class ",k), 1, 0 ) # select probability column

    formula.k= stats::as.formula(paste("ClusM.k ~", paste(COVNAMES, collapse = "+")))
    if(partition=="soft"){
      result.k = stats::glm(formula.k, data = U.cov, weights =  U.cov[,k],
                     family = stats::binomial())

    }else{
      result.k = stats::glm(formula.k, data = U.cov, family = stats::binomial())
    }


    result.k_coeff = data.frame(summary(result.k)$coefficients)

    # output.summary = capture.output(summary(result.k))
    # Coeff.which = which(output.summary=="Coefficients:" )
    # result.k_coeff = output.summary[(Coeff.which+2):(Coeff.which+2+nCOV)]
    #
    # result.k_coeff = matrix(unlist(lapply(result.k_coeff,  function(x) as.numeric(strsplit(x, "\\s+")[[1]][2:5]))), byrow = T, ncol=4)
    # result.k_coeff = data.frame(result.k_coeff)
    colnames(result.k_coeff) = c("Estimate", "Std.error", "z value", "Pr(>|z|)")
    #rownames(result.k_coeff) = c("(Intercept)", COVNAMES)

    output.coff[[k]] = result.k_coeff
    output[[k]] = result.k
  }

  names(output.coff) = paste("Latent Class", 1:num.cluster)
  names(output) = paste("Latent Class", 1:num.cluster)

  return(list(test_results= output.coff,
              binomial_raw= output))

}
