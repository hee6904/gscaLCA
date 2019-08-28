al_gscaLCA = function(MS,z0, bz0, c, nobs, nvar, ntv,nlv, nzct, const,V, W,W0, T,vb,alpha){
  ## MS=1 (mean structure); MS=0 (no mean structure)
  if (MS == 1){
    Z <- bz0
  }else{

    Z <- matrix(ncol=ncol(bz0), nrow=nrow(bz0))

    for (j in 1:nvar){
      Z[,j] <- scale(bz0[,j])*sqrt(length(bz0[,j]))/sqrt(length(bz0[,j])-1)
    }

  }


  ## initial membership
  if( c < 2 ) stop ("The number of cluster should be larger than 1")

  index0 <- ceiling(runif(nobs)*c)

  U <- c()
  for (i in 1:c)
  {
    U <- cbind(U, ifelse(index0 ==i,1,0))
  }

  results <- fclust::FKM(Z, k=c, m=2, startU=U) # k : num of cluster; m:parameter of fuzziness
  U <- results$U

  bi <- matrix(0, nvar*ntv, c) # to compute CI

  #print(paste("Number of Cluster =", c))

  #########################################
  # ALS ALGORITHM FOR GSCA
  #########################################
  itmax = 1000              # maximum number of iterations
  it = 0                   # iteration counter
  ceps = 0.00001           # convergence tolerance
  imp = 10^50             # initial improvement
  f0 = 10^50


  while((it <= itmax) & (imp > ceps))
  {
    it<- it + 1
    U_old <- U

    # Update GSCA parameters
    itmax1 =1000
    it1=0
    ceps1=0.01
    imp1=1000000000000000
    f10 =10^46

    while((it1 <= itmax1) & (imp1 > ceps1))
    {
      it1 <- it1+1
      M <- c(); f1=0; f2=0; #snzw=0; snzt=0;

      for (k in 1:c)
      {
        D_c <- sqrt(U[,k]^const)%*% t(rep(1,ncol(Z)))
        Z_c <- Z*D_c

        Psi_c <- Z_c %*% V[[k]]
        Gamma_c <- Z_c %*% W[[k]]

        ##### STEP 1: UPDATE T (LOADINGS AND PATH COEFFICIETNS)

        vect0 <- as.numeric(T[[k]])
        vecPsi <- as.numeric(Psi_c)

        m0 <- kronecker(diag(ntv), Gamma_c)
        m  <- m0[,nzct]
        #mm <- m0[,nzct1]%*% nzct1.v # mm would be necessary
        #vect <- MASS::ginv(t(m)%*%m)%*%t(m)%*%(vecPsi-mm)

        vect <- MASS::ginv(t(m)%*%m)%*%t(m)%*%(vecPsi)
        vect0[nzct] <- vect

        T[[k]]<- matrix(vect0, nrow=nlv, ncol=ntv )

        ##### STEP 2: UPDATE W (COMPONENT WEIGHTS)
        for(p in 1:nlv)
        {
          w0 <- W[[k]][,p]
          w00 <- W0[,p]

          P <- rep(0, ntv)
          P[nvar+p] <- 1
          Q <- T[[k]][p, ]
          beta <- P-Q

          H <- diag(nlv)
          H[p,p]<- 0

          HH <- diag(ntv)
          HH[nvar+p, nvar+p] <-  0

          Delta <- W[[k]] %*% H %*% T[[k]] - V[[k]] %*% HH
          ZDelta <- Z_c%*%Delta
          vecZDel <- as.vector(ZDelta)

          nzcw <- which(w00 == 99)
          nzw <- length(nzcw)

          if(nzw > 0){
            N <- kronecker(beta, Z_c)[,nzcw]
            w <-  MASS::ginv(t(N)%*% N) %*%t(N) %*% vecZDel
            w0[nzcw]<- w

            zw <- Z_c %*% w0

            if(MS==0){
              w0 <- w0/sqrt(sum(zw^2))
            }

            W[[k]][,p] <- w0
            V[[k]][,nvar+p] <- w0

            #snzw <- snzw +nzw
          }
        }# p in 1:nlv for STEP 2

        ##### STEP 3: OPTIMAL SCALING

        VWT <-  V[[k]] - W[[k]]%*%T[[k]]
        mz0 <- Z_c

        for (i in 1:length(vb))
        {
          ZH <-  diag(nvar)
          J <- vb[i]
          ZH[J, J] <- 0

          K <- -Z_c %*% ZH %*% VWT
          q <- VWT[J, ]
          mz0[, J] <-(K %*% q)/as.numeric(t(q)%*%q) # model prediction for each variable

          mz <- mz0[!is.na(mz0[,J]), J]
          # z0_mn <- z0[!is.na(z0[,J]), J]  #not necessary

          u <- fastDummies::dummy_cols(z0[!is.na(z0[,J]),J])
          u <- as.matrix(u[,paste0(".data_", 1:max(z0[,J]))])

          ## Optimizing based on type
          # NORMINAL
          invind <- MASS::ginv(t(u) %*% u) %*% t(u) %*% mz
          sz <- u %*% invind

          # common part after obtaining sz based on type

          sz0 <- ifelse(is.na(z0[,J]), mz0[,J],sz ) # would be different for missing data; should be checked
          Z_c[,J] <- sz0

          s <- as.numeric(1/sqrt(t(Z_c[,J])%*%Z_c[,J]/nobs))
          Z_c[,J] <- s*Z_c[,J]

        } # for i in 1:leng(vb)


        Psi_c <- Z_c %*% V[[k]]
        Gamma_c <- Z_c %*% W[[k]]

        dif <- Z %*% V[[k]]-  Z%*% W[[k]]%*%T[[k]]  # Z was used in matlab code; Z_c was not work well
        obj_func <- Psi_c - Gamma_c %*% T[[k]]

        f1 <- f1 + psych::tr(t(obj_func)%*%obj_func)
        f2 <- f2 + psych::tr(t(Psi_c)%*%Psi_c)

        M <- cbind(M, rowSums(dif^2))

        bi0 <- V[[k]] - W[[k]]%*%T[[k]]
        bi[,k] <- as.vector(bi0)
        bi[,k] <- bi[,k]/norm( bi[,k], type="2") ## why type 2?

      }# k in 1:c

      imp1 <- abs(f10 - f1)
      #print(paste0("I_LOOP: it1=",it1, ", f10=", round(f10,2), ", f1=", round(f1,2), ", imp1=", round(imp1,7)))
      f10 <- f1

      #snzt<- snzt + nzt # should be before loop for k
    } # second while

    ## STEP 4 Update U

    for (k in 1:c)
    {
      T00 <- (M[,k] %*% t(rep(1, c)))/M
      TT0<- rowSums(T00^alpha)
      U[,k]<- 1/TT0

    }

    imp= f0-f1
    #print(paste0("O.LOOP: it=",it, ", f1=", round(f1,2), ", f0=", round(f0,2), ", imp=", round(imp,7)))
    f0 <- f1

  }#whole while (end ALS)
  #print(b)

  return(list(U=U, bi=bi, f1=f1, f2=f2))
}
