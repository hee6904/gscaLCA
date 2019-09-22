print_gscaLCA= function(c, z0, dat.origin, Boot.num, model.fit, model.fit.result,LCprevalence.result,
                        RespProb.1){
  cat("=========================================================\n")
  cat("LCA by using Fuzzing Clusterwise GSCA\n")
  cat("=========================================================\n")
  cat(paste("Fit for", c, "latent classes:"), "\n",
      paste0("number of used observations: ", nrow(z0)),"\n",
      paste0("number of deleted observation: ", nrow(dat.origin)-nrow(z0)),"\n",
      paste0("number of bootstrap for SE: ",Boot.num - sum(unlist(lapply(model.fit, function(x){is.null(x)})))),"/",Boot.num, "\n")
  #     paste0("number of bootstrap for SE: ",Boot.num - (length(model.fit)-1)),"\n")

  cat("\n")

  cat("MODEL FIT -----------------------------------------------\n",
      "FIT      : ", sprintf("%.4f", model.fit.result[1,"Estimate"]), "\n",
      "AFIT     : ", sprintf("%.4f", model.fit.result[2,"Estimate"]), "\n",
      "FPI      : ", sprintf("%.4f", model.fit.result[3,"Estimate"]), "\n",
      "NCE      : ", sprintf("%.4f", model.fit.result[4,"Estimate"]), "\n", "\n")

  cat("Estimated Latent Class Prevalnces (%) -------------------\n",
      paste0(sprintf("%.2f", LCprevalence.result[,"Percent"]), "%"), "\n", "\n")


  cat("Conditional item response probability -------------------\n ")
  print(lapply(RespProb.1, function(x) {
    x[, "Estimate"] <- sprintf("%.4f", x[, "Estimate"])
    x }))
}

