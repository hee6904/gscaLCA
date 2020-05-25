test_multinomial = function(dat, covnames, membership.1 , num.cluster, ref.clsuter = "MAX",
                            partition )
{
  U.cov = data.frame(cbind(membership.1, dat[,covnames]))

  names(U.cov)[(ncol(membership.1)+1):ncol(U.cov)] = covnames

  if(ref.clsuter=="MAX"){
    Ref.Cluster  = which.max(table(membership.1$label))
  }else if(ref.clsuter=="MIN"){
    Ref.Cluster  = which.min(table(membership.1$label))
  }else if(ref.clsuter=="FIRST"){
    Ref.Cluster  = 1
  }else if(ref.clsuter=="LAST"){
    Ref.Cluster  = num.cluster
  }


  U.cov$label = factor(U.cov$label, c(paste0("Latent Class ", Ref.Cluster),
                                      paste0("Latent Class ", setdiff(1:num.cluster, Ref.Cluster ))),
                       c(paste0("Latent Class ", Ref.Cluster),
                         paste0("Latent Class ", setdiff(1:num.cluster, Ref.Cluster ))))

  formula.lm= stats::as.formula(paste("label ~", paste(covnames, collapse = " + ")))


  if(partition == "hard"){

      test.lm <- nnet::multinom(formula.lm, data = U.cov,
                                trace = FALSE)

  }else {

    wt = apply(U.cov[,grep("Class", names(U.cov))], 1, function(x) x[which.max(x)])

    test.lm <- nnet::multinom(formula.lm, data = U.cov, weights = wt,
                                trace = FALSE)

  }


  test_results = list()


  if(num.cluster==2){

    test_result_tem = cbind(summary(test.lm)$coefficients,
                            summary(test.lm)$standard.errors)
    test_result_tem = data.frame(test_result_tem)

    test_result_tem$V3 = test_result_tem[,1]/test_result_tem[,2]
    test_result_tem$V4 = (1 - stats::pnorm(abs(test_result_tem$V3), 0, 1)) * 2
    names(test_result_tem) = c("Estimate", "Std.error", "z value", "Pr(>|z|)")

    test_results[[1]] = test_result_tem


  }else{

    for(kt in 1:(num.cluster-1))
    {
      test_result_tem = cbind(summary(test.lm)$coefficients[kt,],
                              summary(test.lm)$standard.errors[kt,])
      test_result_tem = data.frame(test_result_tem)

      test_result_tem$V3 = test_result_tem[,1]/test_result_tem[,2]
      test_result_tem$V4 = (1 - stats::pnorm(abs(test_result_tem$V3), 0, 1)) * 2
      names(test_result_tem) = c("Estimate", "Std.error", "z value", "Pr(>|z|)")
      rownames(test_result_tem) = dimnames(summary(test.lm)$coefficients)[[2]]

      test_results[[kt]] = test_result_tem
    }

  }

  label.level = levels(U.cov$label)
  names(test_results) = paste0(label.level[-1]," / ", label.level[1])
  return(list(test_results= test_results,
              multinom_raw= test.lm))

}
