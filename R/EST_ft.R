EST_ft = function(T0, nzt, vect0, ID, LEVELs, loadtype,
                  MS, Z, z0, c, nobs, nvar, ntv, nlv, nzct, const,V, W, W0, T,vb,alpha)
{
  ######################################################
  # INITIALIZATIONS - RANDOM STARTS
  ######################################################

  ## initial matrix W T V
  W <-  list()
  T <- list()
  for ( k in 1:c)
  {
    W[[k]]<- matrix(ifelse(as.vector(W0)==99, runif(length(as.vector(W0)==99)),as.vector(W0)), nrow=nrow(W0))
    T[[k]]<- matrix(ifelse(vect0==99, runif(nzt), vect0), ncol=ncol(T0))
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
  bz0 <- z0     # ORIGINAL DATA

  ##############

  AL_gscaLCA = al_gscaLCA(MS,Z,z0, bz0, c, nobs, nvar, ntv,nlv, nzct, const,V, W,W0, T,vb,alpha)
  U = AL_gscaLCA$U
  bi = AL_gscaLCA$bi
  f1 = AL_gscaLCA$f1
  f2 = AL_gscaLCA$f2

  ################
  ## Results ----------------------------------
  U.0 <- data.frame(U)#[,class.order])
  U.0$label <- paste("Latent Class", apply(U.0, 1, which.max))
  rownames(U.0) = ID

  LCprevalence.1 = Prevalence(U.0, c, nobs)

  membership.1 <- U.0
  RespProb.1 <-  RespItemProb(bz0, membership.1, LEVELs)


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

  model.fit.1 <- c(FIT, AFIT, FPI, MPE)
  names(model.fit.1) <- c("FIT", "AFIT", "FPI", "NCE")

  return (list(model.fit.1 = model.fit.1,
               LCprevalence.1 = LCprevalence.1,
               RespProb.1 = RespProb.1,
               membership.1 = membership.1))

}# EST_ft



