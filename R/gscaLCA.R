#' @title Main function of LCA by using fuzzy clustering GSCA
#'
#' @description This function enables to run LCA based on GSCA algorithm.
#'
#' @param dat Data that you want to implement GSCA.
#' @param varnames A character vector. The names of columns to be used for gscaLCA.
#' @param ID.var  A character element. The name of ID variable. If ID variable is not specified, gscaLCA will find a ID variable in a given data. The ID of observation are automatically created when data set does not have any ID variable. The default is NULL.
#' @param num.cluster A numeric element. Number of cluster to be analyzed. The default is 2.
#' @param num.factor Either "EACH" or "ALLin1"."EACH" indicates that each variable assumes to have latent variable. "ALLin1" indicates that all variables assumes to share one latent variable. The default is "EACH".
#' @param Boot.num   Number of bootstrap. The standard errors of parameters are obtained by bootstrap in GSCA algorithm. The default is 20.
#' @param multiple.Core A logical element. TRUE enables to use multiple cores for the bootstrap. The default is FASLE.
#'
#' @return A list of the used sample size (N), the number of cluster (C), the number of Bootstrap actually used (Boot.num.im), the model fit indices(model.fit), the latent class prevalence (LCprevalence), the item response probability (RespRrob),  the prosterior membership & the predicted class membership (membership), and the graphs of item response probability (plot).
#'
#' @import parallel
#' @import devtools
#' @import ggplot2
#' @import gridExtra
#' @import progress
#' @import doSNOW
#'
#' @importFrom psych tr
#' @importFrom MASS ginv
#' @importFrom stringr str_extract
#' @importFrom fastDummies dummy_cols
#' @importFrom fclust FKM
#' @importFrom stats complete.cases quantile runif sd
#' @importFrom foreach "%do%" "%dopar%" foreach
#' @export
#'
#'
#' @examples
#' \dontrun{
#' # AddHealth data with 2 clusters
#' R2 = gscaLCA(AddHealth, varnames = names(AddHealth)[2:6], num.cluster = 2, Boot.num=0)
#' R2$model.fit      # Model fit
#' R2$LCprevalence   # Latent Class Prevalence
#' R2$RespProb       # Item Reponse Probability
#' R2$membership     # Membership for all observations
#'
#' # TALIS data with 3 clusters
#' T3 = gscaLCA(TALIS, names(TALIS)[2:6], num.factor = "ALLin1", num.cluster = 3, Boot.num=0)
#'
#' }
#'
#' @references Ryoo, J. H., Park, S., & Kim, S. (2019). Categorical latent variable modeling utilizing fuzzy clustering generalized structured component analysis as an alternative to latent class analysis. Behaviormetrika, 1-16.
#'
gscaLCA <- function(dat, varnames=NULL, ID.var=NULL, num.cluster=2,
                    num.factor="EACH", Boot.num=20, multiple.Core = FALSE)
{

  dat.origin = dat

  #options(warn=1)

  ### Input data

  if(is.null(varnames)) stop ("Variable names for analysis are not specified.")
  if(is.null(num.cluster)) stop ("the number of cluster is not specified.")


  # Check whether data is completed or not. If not, use listwise delection was conducted.
  if(sum(complete.cases(dat))!=nrow(dat)){
    warning('Listwise delection was used. Uncompleted data is not available in the current version')
    dat = dat[complete.cases(dat[, varnames]),  ]
  }


  if(!is.null(ID.var)){
    ID = dat[, ID.var]
  }else if(length(grep("id",names(dat), ignore.case=TRUE))==1){
    ID = dat[, grep("id",names(dat), ignore.case=TRUE)]
  }else{

    if(length(grep("id",names(dat), ignore.case=TRUE)) > 1){

      warning('Observation ID was created automatically.')
    }

    ID = 1:nrow(dat)

  }

  Z00 <- dat[, varnames]

  LEVELs = list()
  for(j in 1:ncol(Z00))
  {
    LEVELs[[j]] = levels(as.factor(Z00[,j]))

    if(is.factor(Z00[,j])){
      Z00[,j] = as.numeric(Z00[,j])
    }else if(is.character(Z00[,j])){
      Z00[,j] = as.numeric(as.factor(Z00[,j]))
    }else if(is.numeric(Z00[,j])){
      Z00[,j] =  as.numeric(factor(Z00[,j],
                                   levels=as.numeric(levels(as.factor(Z00[,j]))),
                                   labels= letters[1:length(levels(as.factor(Z00[,j])))])) }
  }




  z0 = data.frame(Z00)

  ######################################################
  # input
  ######################################################


  MS <- 0             # 0 = No mean structure 1 = mean structure



  c <- num.cluster
  const <- 2         # usually 2
  alpha <- 1/(const-1)


  if(num.factor=="EACH"){
    ####### ncol(z0) factor ###################
    loadtype <- rep(1, ncol(z0))             #0 = formative, 1 = reflective

    W0 <- diag(ncol(z0))
    diag(W0)<- 99

    B0 <- matrix(99,  ncol=ncol(z0), nrow=ncol(z0))
    diag(B0) = 0

    # BO.vect <- c()
    # for (r in 1:ncol(z0))
    # {
    #   BO.vect <- c( BO.vect, rep(99, r-1), rep(0, ncol(z0)+1-r ))
    # }
    #
    # B0 <- matrix(BO.vect, nrow=ncol(z0))

  }else if(num.factor=="ALLin1"){

    ###### 1 factor ##############
    loadtype <- c(1)             #0 = formative, 1 = reflective

    W0 <- matrix(rep(99, ncol(z0)), nrow=ncol(z0),byrow=TRUE)

    B0 <- matrix(c(0), nrow=1, byrow=TRUE)
  }



  vb <- 1:ncol(z0) # num of variable (columns of data)

  A0 <- t(W0)

  nobs<- nrow(z0);  nvar <- ncol(z0)     # NOBS = NUM OF OBSERVATIONS, NVAR = NUM OF MANIFESTS
  nlv <- length(loadtype)                 # NUM OF LATENTS
  ntv <- nvar + nlv                       # NUM OF MANIFESTS AND LATENTS

  nnv <- c()          # num of nonmissing values in each variable
  for (j in 1:nvar)
  {
    nnv <- c(nnv, length(!is.na(z0[,j])))
  }


  ######################################################
  # WEIGHT, LOADING AND PATH COEFFICIENT MATRICES
  ######################################################
  for (p in 1:nlv)
  {
    if(loadtype[p] == 0){
      A0[p, ] <- rep(0, nvar)
    }
  }

  T0 <- cbind(A0,B0)
  vect0 <- as.vector(T0)
  nzct <- which(vect0 == 99)           # FREE PARAMETERS IN T
  nzt <- length(nzct)                  # nzt = NUM OF FREE PARAMETERS IN T

  #nzct1 <- which(vect0 != 99)           # FIXED ELEMENTS IN T
  #nzct1.v <- vect0[which(vect0 != 99) ]
  #nzt1 <- length(nzct1)                # nzt1 = NUM OF FIXED ELEMENTS IN T


  EST = EST_ft(T0, nzt, vect0, ID, LEVELs, loadtype,
               MS, z0, c, nobs, nvar, ntv, nlv, nzct, const, W0,vb, alpha)
  model.fit.1    = EST$model.fit.1
  LCprevalence.1 = EST$LCprevalence.1
  RespProb.1     = EST$RespProb.1
  membership.1   = EST$membership.1


  if (Boot.num > 0){
    Boot.Gen = function(z0){
      z0[sample(1:nrow(z0),replace = TRUE),]
    }



   if(isTRUE(multiple.Core)){

     BZ0 = foreach(i=1:Boot.num) %do% Boot.Gen(z0)

     chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

     if (nzchar(chk) && chk == "TRUE") {
       # use 2 cores in CRAN/Travis/AppVeyor
       numCores <- 2L
     } else {
       # use all cores in devtools::test()
       numCores <- parallel::detectCores()
     }



     cl <- makeCluster(numCores)
      registerDoSNOW(cl)


     pb <- progress_bar$new(
       #format = "letter = :letter [:bar] :elapsed | eta: :eta",
       total = Boot.num,    # 100
       width = 60)

     progress_letter <- rep(LETTERS[1:10], 10)  # token reported in progress bar

     # allowing progress bar to be used in foreach -----------------------------
     progress <- function(n){
       pb$tick(tokens = list(letter = progress_letter[n]))
     }

     opts <- list(progress = progress)

     # foreach loop ------------------------------------------------------------


     BOOT.result <-
       foreach(i=1:Boot.num, .options.snow = opts) %dopar% Boot_ft(T0, nzt, vect0, ID, LEVELs, loadtype, LCprevalence.1, RespProb.1, varnames,
                                                                   MS, z0, BZ0[[i]], c, nobs, nvar, ntv, nlv, nzct, const, W0,vb, alpha)

     stopCluster(cl)

     model.fit = lapply(BOOT.result, function(x){x[["model.fit.b"]]})
     RespProb = lapply(BOOT.result, function(x){x[["RespProb.b"]]})
     LCprevalence = lapply(BOOT.result, function(x){x[["LCprevalence.b"]]})

   }else{
     model.fit =list()
     RespProb=list()
     LCprevalence= list()

    pb <- progress_bar$new(total = (Boot.num))
     for (b in 1: Boot.num)
     {
       bz0 = Boot.Gen(z0)

       BOOT.result = Boot_ft(T0, nzt, vect0, ID, LEVELs, loadtype, LCprevalence.1, RespProb.1, varnames,
                   MS, z0, bz0, c, nobs, nvar, ntv, nlv, nzct, const, W0,vb, alpha)


       model.fit[[b]] =  BOOT.result$model.fit.b
       RespProb[[b]] = BOOT.result$RespProb.b
       LCprevalence[[b]] = BOOT.result$LCprevalence.b

       pb$tick()
       Sys.sleep(1 /(Boot.num) )

     }

   }

  }else{

    model.fit= list()
    LCprevalence=list()
    RespProb= list()
    model.fit[[1]]= model.fit.1
    LCprevalence[[1]]= LCprevalence.1
    RespProb[[1]]= RespProb.1
  }


  if(all(unlist(lapply(model.fit, is.null)))|
     all(unlist(lapply(LCprevalence, is.null)))|
     all(unlist(lapply(RespProb, is.null)))){
    model.fit[[1]]= model.fit.1
    LCprevalence[[1]]= LCprevalence.1
    RespProb[[1]]= RespProb.1
  }


  ### Model Fit ####
  model.fit.mat  <- rbind(model.fit.1,
                          matrix(unlist(model.fit), ncol = 4, byrow = TRUE))

  model.fit.result <- cbind(model.fit.1,
                            apply(model.fit.mat, 2, sd),
                            apply(model.fit.mat, 2, function(x){stats::quantile(x, probs=0.025)}),
                            apply(model.fit.mat, 2, function(x){stats::quantile(x, probs=0.975)}))

  colnames(model.fit.result) <- c("Estimate", "SE", "95CI.lower", "95CI.upper")


  ### Latent class Prevalence ####
  LCprevalence.mat <-  rbind ( LCprevalence.1[,1],
                              matrix(unlist(lapply(LCprevalence,
                                           function(x){
                                             x[1:c]})),
                             ncol=c, byrow = TRUE))


  LCprevalence.result <- cbind(LCprevalence.1 ,
                               apply(LCprevalence.mat, 2, sd),
                               apply(LCprevalence.mat, 2, function(x){stats::quantile(x, probs=0.025)}),
                               apply(LCprevalence.mat, 2, function(x){stats::quantile(x, probs=0.975)}))

  colnames(LCprevalence.result) <- c("Percent", "Count","SE",  "95.CI.lower", "95.CI.upper")


  RespProb.results <- RespProb.1

  if(all(unlist(lapply(RespProb, function(x){nrow(x[[1]])}))!=nrow(RespProb.1[[1]]))) stop("An optimalized solution cannot be found; it maybe due to different numbers of classes over bootstrap.")

  RespProb[[Boot.num+1]]=RespProb.1
  RespProb.mat <- matrix( unlist(lapply(RespProb, function(x){
    vect.varnames <- c()
    if(!is.null(x)){
      for(i in varnames)
      {
        vect.varnames <- c(vect.varnames,  unlist(x[[i]][,"Estimate"]))
      }
      vect.varnames
    }

  })), ncol= (sum(LCprevalence.1[,1]!=0)*length(unlist(LEVELs))),byrow=TRUE)




  tem.resprob <-  cbind(apply(RespProb.mat, 2, sd),
                        apply(RespProb.mat, 2, function(x){stats::quantile(x, probs=0.025)}),
                        apply(RespProb.mat, 2, function(x){stats::quantile(x, probs=0.975)}))
  colnames(tem.resprob) <-c("SE",  "95.CI.lower", "95.CI.upper")

  nrow.mat = unlist(lapply(RespProb.results, nrow))
  for(i in 1:length(varnames))
  {
    if(i==1){
      RespProb.results[[i]] <- cbind(RespProb.results[[i]], tem.resprob[1:nrow.mat[1], ])
    }else{
      RespProb.results[[i]] <- cbind(RespProb.results[[i]], tem.resprob[(sum(nrow.mat[1:(i-1)])+1):sum(nrow.mat[1:i]), ])
    }

  }


  Iden.vect= c()

  for (j in 1:(length(LEVELs)-1))
  {
    Iden.vect= c(Iden.vect, identical(LEVELs[[j]], LEVELs[[j+1]]))
  }

  if(all(Iden.vect!=TRUE)){
    P ="Sorry, We don't provide the plot because the variable does not have the same number of categories."
  }else{

    P = list()
    for (j in 1:length(LEVELs[[1]]))
    {
      matYes = c()
      for(l in 1:length(LEVELs))
      {

        matYes = rbind(matYes,
                       subset(RespProb.results[[l]],  RespProb.results[[l]]$Category==LEVELs[[l]][j])[,1:3])
      }

      matYes = cbind(rep(varnames, each=sum(LCprevalence.1[,1]!=0)), matYes)
      matYes$Category=NULL
      names(matYes)[1] ="Type"
      matYes$Type = factor(matYes$Type, levels= varnames, labels=varnames)


      class.numeric = as.numeric(str_extract(unique(matYes$Class), "\\-*\\d+\\.*\\d*"))


      matYes$Class = factor(matYes$Class , unique(matYes$Class),
        paste0(class.numeric,
               " (",sprintf("%.2f", LCprevalence.result[class.numeric,"Percent"]), "%)") )




      ## Line plot by using Aggregated Data (with or without error bar)
      P[[j]]= ggplot(matYes, aes(x= matYes$Type , y= matYes$Estimate, colour= matYes$Class, group=matYes$Class)) +
        geom_line(size=1, aes(linetype=matYes$Class)) +
        geom_point(size=3, aes(shape =matYes$Class))+
        theme_light()+
        ylim(0, 1)+
        #ylab("Probability of Yes")+
        ggtitle(paste("Response:", LEVELs[[1]][j])) +
        theme(
          plot.title = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #axis.title.y = element_text(size =13),
          axis.text.x = element_text(size = 15),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15),
          axis.text.y = element_text(size = 15))



      # print(j)
      # print(p)
    }
  }



  RESULT <- list( N = nrow(z0), C = c,
                  Boot.num.im = Boot.num - sum(unlist(lapply(model.fit, function(x){is.null(x)}))),
                  model.fit = model.fit.result,
                  LCprevalence = LCprevalence.result,
                  RespProb=RespProb.results,
                  membership = membership.1, plot = P)


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

  if(all(Iden.vect==TRUE)){
    if(length(LEVELs[[1]])==2){
      print(P[[1]])
    }else{

      get_legend<-function(myggplot){
        tmp <- ggplot_gtable(ggplot_build(myggplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }

      legend <- get_legend(P[[1]]+theme(legend.position="bottom"))

      P.1 =lapply(P, function(x) x + theme(legend.position="none"))
      P.1[[length(LEVELs[[1]])+1]] =legend


      grid.arrange(grobs = P.1, layout_matrix = matrix(
        c(rep(1:length(LEVELs[[1]]),each=5),
          length(LEVELs)+1),ncol=1))


      #do.call(grid.arrange, c(P.1, nrow= length(LEVELs[[1]])))
    }
  }



  return(RESULT)


}
