
Boot_ft = function(T0, nzt, vect0, ID, LEVELs, loadtype, LCprevalence.1, RespProb.1, varnames,
                   MS,z0, bz0, c, nobs, nvar, ntv, nlv, nzct, const,W0,vb,alpha,A0, num.factor, nInd)
{

  model.fit.b= NULL
  RespProb.b = NULL
  LCprevalence.b = NULL
  ######################################################
  # INITIALIZATIONS - RANDOM STARTS
  ######################################################

  ## initial matrix W T V
  W <-  list()
  T.mat <- list()
  for ( k in 1:c)
  {
    W[[k]]<- matrix(ifelse(as.vector(W0)==99, runif(length(as.vector(W0)==99)),as.vector(W0)), nrow=nrow(W0))
    T.mat[[k]]<- matrix(ifelse(vect0==99, runif(nzt), vect0), ncol=ncol(T0))
  }


  VV <- c()
  for (p in 1:nlv)
  {
    w0 = W[[1]][,p]
    nnw <- sum(w0!=0)

    if(loadtype[p] == 0){
      VV <- c(VV, rep(0,nnw))
    }else{
      VV <- c(VV, rep(1,nnw))
    }
  }
  VV <- diag(VV)

  V <- list()
  for(k in 1:c)
  {
    V[[k]]<- cbind(VV, W[[k]])
  }
  ## input data Z (unstandized or STADANDIZED RAW DATA )
  ## bz0 <- z0[sample(1:nrow(z0),replace = TRUE),]


  ##############
  AL_gscaLCA = al_gscaLCA(MS,z0, bz0, c, nobs, nvar, ntv,nlv, nzct, const,V, W,W0, T.mat,vb,alpha,A0,
                          num.factor, nInd)
  U = AL_gscaLCA$U
  bi = AL_gscaLCA$bi
  f1 = AL_gscaLCA$f1
  f2 = AL_gscaLCA$f2

  # results -----------------------------------------------------------------
  ############################
  # convert to crispy clustering
  ############################

  U.0 <- data.frame(U)
  U.0$label <- paste("Latent Class", apply(U.0, 1, which.max))

  # find cluster order
  U.1.0 <- U.0
  red.c = c

  classfied.tem = Prevalence(U.0, c, nobs)

  if(length(which(LCprevalence.1[,"Percent"]==0))== length(which(classfied.tem[,"Percent"]==0))) {

    if(any(LCprevalence.1[,"Count"]==0)){

      classfied = Prevalence(U.0, c, nobs)

      # if(length(which(LCprevalence.1[,"Percent"]==0))!= length(which(classfied[,"Percent"]==0))) stop("An optimalized solution cannot be found; it maybe due to different numbers of classes over bootstrap.")
      index0 = rep(0, c)

      index0.raw =  which(classfied[,"Count"]==0)
      index0.b1 = which(LCprevalence.1[,"Count"]==0)

      index0[index0.b1]= index0.raw

      all.index = 1:c
      index0[all.index[-intersect(1:c,index0.b1)]]= all.index[-intersect(1:c,index0.raw)]

      U.1.0  = data.frame(U[, index0 ])
      U.1.0$label <- paste("Latent Class", apply(U.1.0, 1, which.max))

      LCprevalence.b = Prevalence(U.1.0, c, nobs)
      red.c = c - length(index0.b1)
    }

    RespProb.raw  <- RespItemProb(bz0, varnames, U.1.0, LEVELs)

#    if(nrow(RespProb.raw[[1]])!= nrow(RespProb.1[[1]])) stop("An optimalized solution cannot be found; it maybe due to different numbers of classes over bootstrap.")


    if(red.c > 1){

      index.mat  <- c()
      for(i in 1:length(varnames))
      {
        for (j in 1:length(LEVELs[[i]]))
        {
          STD = RespProb.1[[varnames[i]]]
          STD.1 = STD[STD$Category==LEVELs[[i]][j], ]

          RAW = RespProb.raw[[i]]
          RAW.1 = RAW[RAW$Category==LEVELs[[i]][j], ]

          STD_RAW = abs(matrix(rep(STD.1$Estimate, each= nrow(STD.1)), nrow=nrow(STD.1)) - RAW.1$Estimate)
          index.mat = rbind(index.mat, STD.1[apply(STD_RAW, 2, which.min),"Class"])

        } # for j
      } # for i

      index.mat = matrix(as.numeric(apply(index.mat, 2,
                                          function(x){stringr::str_extract ( x,"\\-*\\d+\\.*\\d*")})),
                         ncol=ncol(index.mat))

      new.index = as.numeric(apply(index.mat,2,  function(x){ tem = table(x)
      names(tem)[which.max(tem)]}))

      if(length(unique(new.index))!=red.c){

        max.count = apply(index.mat, 2,
                          function(x){tem = table(x)
                          max.tem = max(tem)
                          length(as.numeric(names(tem)[tem==max.tem])) })
        suppressWarnings({
          if(!all(max.count==1)){
            new.index.1 = ifelse(max.count > 1 , NA ,  new.index)

            new.index.prob  = index.mat[, which(is.na(new.index.1))]

            for( p in 1:length(which(is.na(new.index.1))))
            {
              if(length(which(is.na(new.index.1)))==1){
                new.index.prob.vec =new.index.prob
              }else{
                new.index.prob.vec = new.index.prob[,p]
              }

              tem = table(new.index.prob.vec)
              max.tem = max(new.index.prob.vec)
              new.index.prob.1 = new.index.prob.vec[-which(new.index.prob.vec==intersect(new.index.1 , as.numeric(names(tem)[tem==max.tem])))]

              tem.1 = table(new.index.prob.1)
              max.tem.1 = max(tem.1)

              new.latent.class = as.numeric(names(tem.1)[tem.1==max.tem.1])
              if(length(new.latent.class)==1)
              {
                new.index.1[which(is.na(new.index.1))[p]] <- new.latent.class
              }

            }# for loop of p
            new.index = new.index.1
          }

        })

      }else{
        max.count <- rep(2, red.c)
      }  # length(unique(new.index))!=c

      if(sum(is.na(new.index)) < 1 & !all(max.count==1) & length(unique(new.index))==red.c){
        new.index.full <-c()
        for(k in 1:c)
        {
          new.index.full <- c(new.index.full, grep(new.index[k], RespProb.raw[[1]]$Class))
        }

        RespProb.raw.new = lapply(RespProb.raw, function(x){x[new.index.full,]})

        for(i in 1:length(varnames))
        {
          STD = RespProb.1[[varnames[i]]]
          RespProb.raw.new[[i]]$Class= STD$Class
        }

        RespProb.b <- RespProb.raw.new

        if(red.c !=c){
          index0[all.index[-intersect(1:c,index0.b1)]] = new.index
          U.1 = data.frame(U[, index0 ])
        }else{
          U.1 = data.frame(U[, new.index ])
        }

        U.1$label <- paste("Latent Class", apply(U.1, 1, which.max))

        LCprevalence.b <-   Prevalence(U.1, c, nobs)
      }

    }else{ # if red.c > 1

      new.index= rep(1, red.c)
      max.count <- rep(2, red.c)

      for(i in 1:length(varnames))
      {
        STD = RespProb.1[[varnames[i]]]
        RespProb.raw[[i]]$Class= STD$Class
      }
      RespProb.b = RespProb.raw

    }




    ############################
    #    % MODEL FIT MEASURES
    ##############################

    ## cluster validity criteria
    bibj <- t(bi)%*%bi - diag(c)
    bibj <- as.vector(bibj)[which(bibj!=0)]
    CI <- f1/(nobs*max(bibj))
    PC <-  sum(U^2)/nobs

    #print(paste0("Num of Cluster=", c, "; CI=", CI,"; PC=",PC ))

    FPI <- 1-(c*PC-1)/(c-1)
    PE <- -sum(U*log(U))/nobs           # Partition Entropy
    MPE <- PE/log(c)                    # Modified Partition Entropy


    ## MODEL FIT
    #NITER <- it

    DF <- nobs*nvar
    npar <- c*(nzt + sum(W[[1]]!=0)) # not sure

    FIT <- 1-f1/f2
    AFIT = 1 - ((1-FIT)*DF/(DF - npar))

    if((sum(is.na(new.index)) < 1 &  !all(max.count==1) & length(unique(new.index))==red.c)|red.c==1){
      model.fit.b <- c(FIT, AFIT, FPI, MPE)
      names(model.fit.b) <- c("FIT", "AFIT", "FPI", "NCE")
    }

    return(list(model.fit.b= model.fit.b,
                RespProb.b = RespProb.b,
                LCprevalence.b = LCprevalence.b))
  }else{

    list(model.fit.b= NULL,
         RespProb.b = NULL,
         LCprevalence.b = NULL)
  }



}#Boot_ft

