#' @title Main function of gscaLCA by using fuzzy clustering GSCA
#'
#' @description Fitting a component-based LCA by utilizing fuzzy clustering GSCA algorithm.
#'
#' @param dat         Data that you want to fit the gscaLCA function into.
#' @param varnames    A character vector. The names of columns to be used in the gscaLCA function.
#' @param ID.var      A character element. The name of ID variable. If ID variable is not specified, the gscaLCA function will search an ID variable in the given data. The ID of observations will be automatically generated as a numeric variable if the data set does not include any ID variable. The default is NULL.
#' @param num.class A numeric element. The number of classes to be identified The default is 2.
#' @param num.factor  Either "EACH" or "ALLin1"."EACH" specifies the sitatuion that each indicator is assumed to be its phantom latent variable. "ALLin1" indicates that all variables are assumed to be explained by a common latent variable. The default is "EACH".
#' @param Boot.num    The number of bootstraps. The standard errors of parameters are computed from the bootstrap within the gscaLCA algorithm. The default is 20.
#' @param multiple.Core A logical element. TRUE enables to use multiple cores for the bootstrap wehn they are available. The default is \code{FASLE}.
#' @param covnames    A character vector of covariates. The covariates are used when latent class regression (LCR) is fitted.
#' @param cov.model   A numeric vector. The indicator function of latent class regression (LCR) that covariates are involved in fitting the fuzzy clustering GSCA. 1 if gscaLCA is for LCR and otherwise 0.
#' @param multinomial.ref A character element. Options of \code{MAX}, \code{MIX}, \code{FIRST}, and \code{LAST} are available for setting a reference group. The default is \code{MAX}.
#'
#' @return A list of the sample size (N), the number of cluster (C), the number of bootstraps (Boot.num/Boot.num.im), the model fit indices (model.fit), the latent class prevalence (LCprevalence), the item response probability (RespProb), the posterior membership & the predicted class membership (membership), and the graphs of item response probability (plot). When it include covariates, the regression results are also provided.
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
#' @importFrom ggplot2 ggplot
#' @export
#'
#'
#' @examples
#'
#' #AddHealth data with 3 clusters with 500 samples
#' AH.sample= AddHealth[1:500,]
#' R3 = gscaLCA (dat = AH.sample,
#'                varnames = names(AddHealth)[2:6],
#'                ID.var = "AID",
#'                num.class = 3,
#'                num.factor = "EACH",
#'                Boot.num = 0)
#' summary(R3)
#' R3$model.fit      # Model fit
#' R3$LCprevalence   # Latent Class Prevalence
#' R3$RespProb       # Item Response Probability
#' head(R3$membership)     # Membership for all observations
#'
#' # AddHealth data with 3 clusters with 500 samples with two covariates
#' R3_2C = gscaLCA (dat = AH.sample,
#'                  varnames = names(AddHealth)[2:6],
#'                  ID.var = "AID",
#'                  num.class = 3,
#'                  num.factor = "EACH",
#'                  Boot.num = 0,
#'                  multiple.Core = FALSE,
#'                  covnames = names(AddHealth)[7:8], # Gender and Edu
#'                  cov.model = c(1, 0),   # Only Gender varaible is added to the gscaLCR.
#'                  multinomial.ref = "MAX")
#'
#' # To print with the results of multinomial regression with hard partitioning of the gscaLCR,
#' # use the option of "multinomial.hard".
#' summary(R3_2C, "multinomial.hard")
#'
#' \donttest{
#' # AddHealth data with 2 clusters with 20 bootstraps
#' R2 = gscaLCA(AddHealth,
#'              varnames = names(AddHealth)[2:6],
#'              num.class = 2,
#'              Boot.num = 20,
#'              multiple.Core = FALSE) # "multiple.Core = TRUE" is recommended.
#' # TALIS data with 3 clusters with 20 bootstraps and the "ALLin1" option
#' T3 = gscaLCA(TALIS,
#'              varnames = names(TALIS)[2:6],
#'              num.class = 3,
#'              num.factor = "ALLin1",
#'              Boot.num = 20,
#'              multiple.Core = FALSE) # "multiple.Core = TRUE" is recommended.
#'
#'}
#'
#' @references Ryoo, J. H., Park, S., & Kim, S. (2019). Categorical latent variable modeling utilizing fuzzy clustering generalized structured component analysis as an alternative to latent class analysis. Behaviormetrika, 47, 291-306. https://doi.org/10.1007/s41237-019-00084-6
#'
gscaLCA <- function(dat, varnames = NULL,  ID.var = NULL, num.class = 2,
                    num.factor = "EACH", Boot.num = 20, multiple.Core = FALSE,
                    covnames = NULL, cov.model = NULL, multinomial.ref = "MAX")
{
  num.cluster = num.class

  if(is.null(varnames)) stop ("Variable names for analysis are not specified.")
  if(length(intersect(varnames, names(dat)))!=length(varnames)) stop ("Variable names for analysis are not in the data set")
  if(!is.null(covnames)){
    if(length(intersect(varnames, names(dat)))!=length(varnames)) stop ("Covariates names for analysis are not in the data set")
    if(length(grep("Class", covnames)) != 0) stop ("Please change covariates names which do not including \"Class\".")
  }

  if(is.null(num.cluster)) stop ("The number of class is not specified.")
  if(!(num.factor =="EACH" | num.factor =="ALLin1")) stop ("Please check the option `num.factor`")
  if(!is.null(covnames)){
    if(length(covnames) != length(cov.model)) stop ("Please check the objects covnames and cov.model")
  }

  dat.origin = dat
  nobs.origin = nrow(dat.origin)

  # Check whether data is completed or not. If not, use listwise delection was conducted.
  if(sum(complete.cases(dat))!=nrow(dat)){
    print('Listwise deletion was applied for cases whose variables used in gscaLCA have any missing.')
    dat = dat[complete.cases(dat[, c(varnames,covnames)]),  ]
  }


  # Set up ID names
  if(!is.null(ID.var)){
    ID = dat[, ID.var]
  }else if(length(grep("id",names(dat), ignore.case=TRUE))==1){
    ID = dat[, grep("id",names(dat), ignore.case=TRUE)]
  }else{

    if(length(grep("id",names(dat), ignore.case=TRUE)) > 1){

      warning('Observation ID was created automatically.')
    }

    ID = 1:nrow(dat)
    dat$ID = ID

  }


  dat.cov = dat


  # Variables' location for optimization
  # vb = c(which(names(dat) %in% varnames),
  #        which(names(dat) %in% covnames)[sapply(dat[,covnames],is.factor)])
  #vb <- 1:ncol(z0) # num of variable (columns of data)



  # Set up the number of variables
  nInd <- length(varnames)

  if(is.null(covnames)){
    nCov= 0
    nCOV = 0}else{
      COVNAMES <- covnames
      nCOV <- length(COVNAMES)
      covnames <- covnames[cov.model]
      nCov <- length(covnames)}

  dat = dat[, c(varnames,covnames)]
  vb = c(which(names(dat) %in% varnames))

  nvar = nInd + nCov
  nobs = nrow(dat)



  Z00 <- dat

  LEVELs = list()
  for(j in 1:ncol(Z00))
  {
    if(j %in% 1:nInd){
      LEVELs[[j]] = levels(as.factor(Z00[,j]))
    }

#    LEVELs[[j]] = levels(as.factor(Z00[,j]))

    if(is.factor(Z00[,j])){
      Z00[,j] = as.numeric(Z00[,j])
    }else if(is.character(Z00[,j])){
      Z00[,j] = as.numeric(as.factor(Z00[,j]))
    }else if(is.numeric(Z00[,j])){
      Z00[,j] =  as.numeric(factor(Z00[,j],
                                   levels=as.numeric(levels(as.factor(Z00[,j]))),
                                   labels= letters[1:length(levels(as.factor(Z00[,j])))])) }
  }

  names(LEVELs) = varnames


  z0 = data.frame(Z00)

  ######################################################
  # input
  ######################################################


  MS <- 0             # 0 = No mean structure 1 = mean structure



  c <- num.cluster
  const <- 2         # usually 2
  alpha <- 1/(const-1)


  if(num.factor=="EACH"){
    nlv <- nInd + nCov               # number of latent variable without covarite

    ####### num.ind factor ###################
    #    loadtype <- c(rep(1, nInd), rep(0, nCov))             #0 = formative, 1 = reflective
    loadtype <- rep(1, nlv)
    W0 <- diag(nInd)
    diag(W0)<- 99
    B0 <- matrix(99,  ncol=nInd, nrow=nInd)
    diag(B0) <- 0

    if(nCov!=0)
    {
      W0.C = diag(nCov)
      diag(W0.C) = 99

      W0 <- cbind(rbind(W0, matrix(0, nrow = nCov, ncol = nInd)),
                  rbind(matrix (0, nrow = nInd, ncol = nCov), W0.C))
      B0 <- cbind(rbind(B0, matrix(99, nrow = nCov, ncol = nInd)),
                  matrix(0, nrow = nlv, ncol = nCov))

      # B0 <- matrix(99,  ncol=nCov+nInd, nrow=nCov+nInd)
      # diag(B0) = 0
      # B0[1:nInd, (nInd+1):(nCov+nInd)] = 0
    }

  }else if(num.factor=="ALLin1"){
    nlv <- 1 + nCov

    ###### 1 factor ##############
    #    loadtype <- c(1, rep(0, nCov))             #0 = formative, 1 = reflective
    loadtype <- rep(1, nlv)

    W0 <- matrix(rep(99, nInd), nrow=nInd, byrow=TRUE)

    B0 <- matrix(c(0), nrow=1, byrow=TRUE)
    #diag(B0) <- 0

    if(nCov!=0)
    {
      W0.C = diag(nCov)
      diag(W0.C) = 99

      W0 <- cbind(rbind(W0, matrix(0, nrow = nCov, ncol = 1)),
                  rbind(matrix (0, nrow = nInd, ncol = nCov), W0.C))
      B0 <- cbind(rbind(B0, matrix(99, nrow = nCov, ncol = 1)),
                  matrix(0, nrow = nlv, ncol = nCov))

    #   W0 <- cbind(c(as.numeric(W0), rep(0,  nCov)),
    #               rbind(matrix (0, nrow = nInd, ncol = nCov), diag(nCov)))
    # B0 <-  matrix(99, nrow = nlv, ncol = nlv)
    # diag(B0) <- 0
    # B0[1,] <- 0
    }

  }



  A0 <- t(W0)

  nobs <- nrow(z0);  nvar <- ncol(z0)     # NOBS = NUM OF OBSERVATIONS, NVAR = NUM OF MANIFESTS
  nlv <- length(loadtype)                 # NUM OF LATENTS
  ntv <- nvar + nlv                       # NUM OF MANIFESTS AND LATENTS

  # nnv <- c()          # num of nonmissing values in each variable
  # for (j in 1:nvar)
  # {
  #   nnv <- c(nnv, length(!is.na(z0[,j])))
  # }


  ######################################################
  # WEIGHT, LOADING AND PATH COEFFICIENT MATRICES
  ######################################################
  # for (p in 1:nlv)
  # {
  #   if(loadtype[p] == 0){
  #     A0[p, ] <- rep(0, nvar)
  #   }
  # }

  T0 <- cbind(A0,B0)
  vect0 <- as.vector(T0)
  nzct <- which(vect0 == 99)           # FREE PARAMETERS IN T
  nzt <- length(nzct)                  # nzt = NUM OF FREE PARAMETERS IN T

  #nzct1 <- which(vect0 != 99)           # FIXED ELEMENTS IN T
  #nzct1.v <- vect0[which(vect0 != 99) ]
  #nzt1 <- length(nzct1)                # nzt1 = NUM OF FIXED ELEMENTS IN T


  EST = EST_ft(T0, nzt, vect0, ID, LEVELs, loadtype,
               MS, z0, c, nobs, nvar, ntv, nlv, nzct, const, W0,vb, alpha, varnames, A0,
               num.factor, nInd)
  model.fit.1    = EST$model.fit.1
  LCprevalence.1 = EST$LCprevalence.1
  RespProb.1     = EST$RespProb.1
  membership.1   = EST$membership.1
  it.out = EST$it.out
  it.in = EST$it.in
  A.mat = EST$A.mat
  B.mat = EST$B.mat
  W.mat = EST$W.mat

  names(A.mat) = paste0("Class", 1:num.cluster)
  names(B.mat) = paste0("Class", 1:num.cluster)
  names(W.mat) = paste0("Class", 1:num.cluster)

  if(length(which(LCprevalence.1[, "Count"]==0))>0){
    if(Boot.num != 0 & nCOV != 0){
      print("The estimated number of classes is not consistent with the assigned number of classes; thus, the bootstrap for SE and the fitting a regression model are not implemented.")
      Boot.num = 0
    }else if(Boot.num == 0 & nCOV != 0){
      print("The estimated number of classes is not consistent with the assigned number of classes; thus, the fitting a regression model is not implemented.")
    }else if(Boot.num != 0 & nCOV == 0){
      print("The estimated number of classes is not consistent with the assigned number of classes; thus, the bootstrap for SE is not implemented.")
      Boot.num = 0
    }

  }

  if (Boot.num > 0){

    Boot.Gen = function(z0){
      z0[sample(1:nrow(z0),replace = TRUE),]
    }

   if(multiple.Core){

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


     pb <- progress::progress_bar$new(
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
                                                                   MS, z0, BZ0[[i]], c, nobs, nvar, ntv, nlv, nzct, const, W0,vb, alpha,A0,
                                                                   num.factor, nInd)

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
                   MS, z0, bz0, c, nobs, nvar, ntv, nlv, nzct, const, W0,vb, alpha, A0, num.factor, nInd)


       model.fit[[b]] =  BOOT.result$model.fit.b
       RespProb[[b]] = BOOT.result$RespProb.b
       LCprevalence[[b]] = BOOT.result$LCprevalence.b

       pb$tick()
       Sys.sleep(1 /(Boot.num) )

     }# b in 1:Boot.num

   }# multiple core or not

  }else{

    model.fit= list()
    LCprevalence=list()
    RespProb= list()
    model.fit[[1]]= model.fit.1
    LCprevalence[[1]]= LCprevalence.1
    RespProb[[1]]= RespProb.1
  } # Boot.num > 0

  names(membership.1)[1:num.cluster] = paste0("Class", 1:num.cluster)


  if(all(unlist(lapply(model.fit, is.null)))|
     all(unlist(lapply(LCprevalence, is.null)))|
     all(unlist(lapply(RespProb, is.null)))){
    model.fit[[1]]= model.fit.1
    LCprevalence[[1]]= LCprevalence.1
    RespProb[[1]]= RespProb.1
  }

  # if(nCOV > 0 & length(unique(membership.1$label)) != num.cluster){
  #   print("The estimated number of classes is not consistent with the assigned number of classes, thus the fitting a regression model is not implemented.")
  # }

  if(nCOV > 0 & length(unique(membership.1$label)) == num.cluster){

    ## hard ##
      multinom_result.hard = test_multinomial(dat.cov, COVNAMES, membership.1, num.cluster, multinomial.ref,
                                        partition = "hard")
      cov_results.multi.hard = multinom_result.hard$test_results
      cov_results_raw.multi.hard =multinom_result.hard$multinom_raw

      binom_result.hard = test_binomial(dat.cov, COVNAMES, membership.1 , num.cluster,
                                        partition = "hard")
      cov_results.bin.hard = binom_result.hard$test_results
      cov_results_raw.bin.hard = binom_result.hard$binomial_raw

    ## soft ##
      multinom_result.soft = test_multinomial(dat.cov, COVNAMES, membership.1, num.cluster, multinomial.ref,
                                              partition = "soft")
      cov_results.multi.soft = multinom_result.soft$test_results
      cov_results_raw.multi.soft = multinom_result.soft$multinom_raw


      binom_result.soft = suppressWarnings( test_binomial(dat.cov, COVNAMES, membership.1 , num.cluster,
                                                         partition = "soft"))
      cov_results.bin.soft = binom_result.soft$test_results
      cov_results_raw.bin.soft = binom_result.soft$binomial_raw

      # soft_result = test_soft(dat.cov, COVNAMES, membership.1 , num.cluster)
      # cov_results.lm = soft_result$test_results
      # cov_results_raw.lm = soft_result$glm_raw
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

#  if(all(unlist(lapply(RespProb, function(x){nrow(x[[1]])}))!=nrow(RespProb.1[[1]]))) stop("An optimalized solution cannot be found; it maybe due to different numbers of classes over bootstrap.")

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
  colnames(tem.resprob) <- c("SE",  "95.CI.lower", "95.CI.upper")

  nrow.mat = unlist(lapply(RespProb.results, nrow))
  for(i in 1:length(varnames))
  {
    if(i==1){
      RespProb.results[[i]] <- cbind(RespProb.results[[i]], tem.resprob[1:nrow.mat[1], ])
    }else{
      RespProb.results[[i]] <- cbind(RespProb.results[[i]], tem.resprob[(sum(nrow.mat[1:(i-1)])+1):sum(nrow.mat[1:i]), ])
    }

  }


  Iden.vect = c()

  for (j in 1:(length(LEVELs)-1))
  {
    Iden.vect = c(Iden.vect, identical(LEVELs[[j]], LEVELs[[j+1]]))
  }
  all.Levels.equal = all(Iden.vect==TRUE)

  if(all(Iden.vect!=TRUE)){
    P ="Sorry, We don't provide the plot because the variable does not have the same number of categories."
  }else{
    Mat.YES = list()

    for (j in 1:length(LEVELs[[1]]))
    {
      matYes = c()
      for(l in 1:length(LEVELs))
      {

        matYes = rbind(matYes,
                       subset(RespProb.results[[l]],  RespProb.results[[l]]$Category==LEVELs[[l]][j])[,1:3])
      }

      matYes = cbind(rep(varnames, each=sum(LCprevalence.1[,1]!=0)), matYes)
      names(matYes)[1] ="Type"
      matYes$Type = factor(matYes$Type, levels= varnames, labels=varnames)


      class.numeric = as.numeric(stringr::str_extract(unique(matYes$Class), "\\-*\\d+\\.*\\d*"))


      matYes$Class = factor(matYes$Class , unique(matYes$Class),
                            paste0(class.numeric,
                                   " (",sprintf("%.2f", LCprevalence.result[class.numeric,"Percent"]), "%)") )

      # Class= matYes$Class

      Mat.YES[[j]]= matYes
    }

    P = lapply(Mat.YES, function(mat.yes){
        Type = mat.yes$Type
        Estimate = mat.yes$Estimate
        Class = mat.yes$Class
        Category = mat.yes$Category

        ggplot(data = mat.yes , aes(x= Type , y= Estimate,
                      colour= Class, group=Class)) +
        geom_line(size=1, aes(linetype=Class)) +
        geom_point(size=3, aes(shape =Class))+
        theme_light()+
        ylim(0, 1)+
        #ylab("Probability of Yes")+
        ggtitle(paste("Response:", unique(Category))) +
        theme(
          plot.title = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #axis.title.y = element_text(size =13),
          axis.text.x = element_text(size = 15),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15),
          axis.text.y = element_text(size = 15))+
        labs(color='Latent Class', shape="Latent Class", linetype="Latent Class")
       })



  }

  Boot.num.im = Boot.num - sum(unlist(lapply(model.fit, function(x){is.null(x)})))

  if(nCOV > 0 & length(unique(membership.1$label)) == num.cluster){
    RESULT <- list( N = nobs, N.origin = nobs.origin, LEVELs = LEVELs, all.Levels.equal = all.Levels.equal,
                    num.class = num.cluster,
                    Boot.num = Boot.num,
                    Boot.num.im = Boot.num.im,
                    model.fit = model.fit.result,
                    LCprevalence = LCprevalence.result,
                    RespProb = RespProb.results,
                    it.in= it.in,
                    it.out = it.out,
                    membership = membership.1, plot = P,
                    cov_results.multi.hard = cov_results.multi.hard,
                    cov_results_raw.multi.hard  = cov_results_raw.multi.hard,

                    cov_results.bin.hard = cov_results.bin.hard,
                    cov_results_raw.bin.hard = cov_results_raw.bin.hard,

                    cov_results.multi.soft = cov_results.multi.soft,
                    cov_results_raw.multi.soft  = cov_results_raw.multi.soft,

                    cov_results.bin.soft = cov_results.bin.soft,
                    cov_results_raw.bin.soft = cov_results_raw.bin.soft,
                    #cov_results.lm = cov_results.lm,
                    #cov_results_raw = cov_results_raw,
                    A.mat = A.mat,
                    B.mat = B.mat,
                    W.mat = W.mat,
                    used.dat = dat.cov)

    # if(verbose) print_gscaLCA(c, nobs, nobs.origin, Boot.num, Boot.num.im, model.fit.result,LCprevalence.result,
    #                           RespProb.1, cov_results.multi, cov_results.bin, cov_results.lm)

  }else{
    RESULT <- list( N = nobs, N.origin = nobs.origin, LEVELs = LEVELs, all.Levels.equal = all.Levels.equal,
                    num.class = num.cluster,
                    Boot.num = Boot.num,
                    Boot.num.im = Boot.num.im,
                    model.fit = model.fit.result,
                    LCprevalence = LCprevalence.result,
                    RespProb = RespProb.results,
                    it.in= it.in,
                    it.out = it.out,
                    membership = membership.1, plot = P,
                    A.mat = A.mat,
                    B.mat = B.mat,
                    W.mat = W.mat,
                    used.dat = dat.cov)
    # if(verbose) print_gscaLCA(c, nobs, nobs.origin, Boot.num, Boot.num.im,  model.fit.result,LCprevalence.result,
    #                           RespProb.1)
  }


  # if(graphs_print) print_graph_gscaLCA (all.Levels.equal, LEVELs, P)

  class(RESULT) <- 'gscaLCA'
  #RESULT$call <- match.call()

  return(RESULT)

}
