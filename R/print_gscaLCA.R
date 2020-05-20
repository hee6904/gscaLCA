print_gscaLCA = function(c, nobs, nobs.origin, Boot.num, Boot.num.im, model.fit.result,LCprevalence.result,
                        RespProb.1,  
                        cov_results.multi.hard = NULL, cov_results.bin.hard = NULL, #cov_results.lm = NULL,   
                        cov_results.multi.soft = NULL, cov_results.bin.soft = NULL, 
                        print.cov.output = NULL){
  cat("=========================================================\n")
  cat("LCA by using Fuzzing Clusterwise GSCA\n")
  cat("=========================================================\n")
  cat(paste("Fit for", c, "latent classes:"), "\n",
      paste0("number of used observations: ", nobs), "\n",
      paste0("number of deleted observation: ", nobs.origin - nobs),"\n",
      paste0("number of bootstrap for SE: ",Boot.num.im ),"/",Boot.num, "\n")
  #     paste0("number of bootstrap for SE: ",Boot.num - (length(model.fit)-1)),"\n")

  cat("\n")

  cat("MODEL FIT -----------------------------------------------\n",
      "FIT      : ", sprintf("%.4f", model.fit.result[1,"Estimate"]), "\n",
      "AFIT     : ", sprintf("%.4f", model.fit.result[2,"Estimate"]), "\n",
      "FPI      : ", sprintf("%.4f", model.fit.result[3,"Estimate"]), "\n",
      "NCE      : ", sprintf("%.4f", model.fit.result[4,"Estimate"]), "\n", "\n")

  cat("Estimated Latent Class Prevalances (%) ------------------\n",
      paste0(sprintf("%.2f", LCprevalence.result[,"Percent"]), "%"), "\n", "\n")


  cat("Conditional Item Response Probability -------------------\n ")
  print(lapply(RespProb.1, function(x) {
                                  x[, "Estimate"] <- sprintf("%.4f", x[, "Estimate"])
                                  x[,1:3] }))
  
  if(!is.null(print.cov.output) & !is.null(cov_results.multi.hard)){
    cat("Relationship Between Prevalances and Covariates -------\n ")

    if(print.cov.output == "multinomial.hard"){
    cat("Multinomial logistic regression is applied with hard partitioning \n ")  
      cov_results.print = lapply(cov_results.multi.hard, function(y) 
        apply(y, 2, function(x){sprintf("%.4f", x)}))
      

    }else if(print.cov.output == "binomial.hard"){
    cat("Binomial logistic regression is applied with hard partitioning \n ")  
     
      cov_results.print = lapply(cov_results.bin.hard, function(y) 
        apply(y, 2, function(x){sprintf("%.4f", x)}))
     
    }else if (print.cov.output == "multinomial.soft"){
      cat("Multinomial logistic regression is applied with soft partitioning \n ")  
      cov_results.print = lapply(cov_results.multi.soft, function(y) 
        apply(y, 2, function(x){sprintf("%.4f", x)}))
      
      
    }else if(print.cov.output == "binomial.soft"){
      cat("Binomial logistic regression is applied with soft partitioning \n ")  
      
      cov_results.print = lapply(cov_results.bin.soft, function(y) 
        apply(y, 2, function(x){sprintf("%.4f", x)}))
      
    } 
    # }else if(print.cov.output == "lm"){
    # cat("Linear regression is applied with Soft partitioning\n ")  
    #   
    #   cov_results.print = lapply(cov_results.lm, function(y)
    #     apply(y, 2, function(x){sprintf("%.4f", x)}))
    #   
    # }
    
    for(i in 1:length(cov_results.print))
    {
      rownames(cov_results.print[[i]])= rownames(cov_results.multi.hard[[1]])
      cov_results.print[[i]] =  data.frame(cov_results.print[[i]])
      colnames(cov_results.print[[i]])[4]= 'P-value'
    }
    
    print(cov_results.print)

  }
  
}

