#' The 2nd and 3rd step of gscaLCA, which are the partitioning and fitting regression
#' @description The 2nd and 3rd step of gscaLCA, which are the partitioning and fitting regression in the latent class regression.
#'
#' @param results.obj the results of gscaLCA.
#' @param covnames A character vector of covariates. The covariates are used when latent class regression (LCR) is fitted.
#' @param multinomial.test.ref A character element. Options of \code{MAX}, \code{MIX}, \code{FIRST}, and \code{LAST} are available for setting a reference group. The default is \code{MAX}.
#'
#' @return Results of the gscaLCR, fitting regression after partioning in addtion to gscaLCA results.
#' @export
#'
#' @examples
#' R2 = gscaLCA (dat = AddHealth[1:500, ], # Data has to include the possible covarite to run gscaLCR
#'                varnames = names(AddHealth)[2:6],
#'                ID.var = "AID",
#'                num.cluster = 3,
#'                num.factor = "EACH",
#'                Boot.num = 0,
#'                multiple.Core = F)
#'
#' R2.gender = gscaLCR (R2, covnames = "Gender")
#' summary(R2.gender,  "multinomial.hard") # hard partitioning with multinomial regression
#' summary(R2.gender,  "multinomial.soft") # soft partitioning with multinomial regression
#' summary(R2.gender,  "binomial.hard")    # hard partitioning with binomial regression
#' summary(R2.gender,  "binomial.soft")    # soft partitioning with binomial regression
#'
gscaLCR = function(results.obj, covnames, multinomial.test.ref = "MAX")
{


  dat.cov = results.obj$used.dat
  membership.1 = results.obj$membership
  num.cluster = results.obj$num.cluster
  COVNAMES = covnames
  if(!all(COVNAMES %in% names(dat.cov)))stop("Please check the data has the covariates that you assign")
  # if(all(COVNAMES %in% names(used.dat))){
  #   dat.cov = used.dat
  # }else{
  #
  #   if(!is.null(ID.var)){
  #    dat.cov = dat[which(dat[, ID.var], rownames(membership.1)),]
  #   }else if(length(grep("id",names(dat), ignore.case=TRUE))==1){
  #     ID = dat[, grep("id",names(dat), ignore.case=TRUE)]
  #     dat.cov = dat[which(ID, rownames(membership.1)),]
  #   }else{
  #
  #    stop("something... more ")
  #   }
  #
  # }

  # # Check whether data is completed or not. If not, use listwise delection was conducted.
  # if(sum(complete.cases(dat.cov))!=nrow(dat.cov)){
  #   print('Listwise deletion was used. A option for incompleted data is not available in the current version')
  #   dat.cov = dat.cov[complete.cases(dat.cov[, c(varnames,covnames)]),  ]
  # }
  #dat.cov = dat.cov


  ## hard ##
  multinom_result.hard = test_multinomial(dat.cov, COVNAMES, membership.1, num.cluster, multinomial.test.ref,
                                          partition = "hard")
  cov_results.multi.hard = multinom_result.hard$test_results
  cov_results_raw.multi.hard =multinom_result.hard$multinom_raw

  binom_result.hard = test_binomial(dat.cov, COVNAMES, membership.1 , num.cluster,
                                    partition = "hard")
  cov_results.bin.hard = binom_result.hard$test_results
  cov_results_raw.bin.hard = binom_result.hard$binomial_raw

  ## soft ##
  multinom_result.soft = test_multinomial(dat.cov, COVNAMES, membership.1, num.cluster, multinomial.test.ref,
                                          partition = "soft")
  cov_results.multi.soft = multinom_result.soft$test_results
  cov_results_raw.multi.soft = multinom_result.soft$multinom_raw


  binom_result.soft = suppressWarnings( test_binomial(dat.cov, COVNAMES, membership.1 , num.cluster,
                                                      partition = "soft"))
  cov_results.bin.soft = binom_result.soft$test_results
  cov_results_raw.bin.soft = binom_result.soft$binomial_raw

  RESULT = list(    N = results.obj$N, N.origin = results.obj$N.origin,
                    LEVELs = results.obj$LEVELs,
                    all.Levels.equal = results.obj$all.Levels.equal,
                    num.cluster = results.obj$num.cluster,
                    Boot.num.im = results.obj$Boot.num.im,
                    model.fit = results.obj$model.fit,
                    LCprevalence = results.obj$LCprevalence,
                    RespProb = results.obj$RespProb,
                    it.in= results.obj$it.in,
                    it.out = results.obj$it.out,
                    membership = results.obj$membership,
                    plot = results.obj$plot,
                    A.mat = results.obj$A.mat,
                    B.mat = results.obj$B.mat,
                    W.mat = results.obj$W.mat,
                    used.dat = results.obj$used.dat,

                    cov_results.multi.hard = cov_results.multi.hard,
                    cov_results_raw.multi.hard  = cov_results_raw.multi.hard,

                    cov_results.bin.hard = cov_results.bin.hard,
                    cov_results_raw.bin.hard = cov_results_raw.bin.hard,

                    cov_results.multi.soft = cov_results.multi.soft,
                    cov_results_raw.multi.soft  = cov_results_raw.multi.soft,

                    cov_results.bin.soft = cov_results.bin.soft,
                    cov_results_raw.bin.soft = cov_results_raw.bin.soft)

  class(RESULT) <- 'gscaLCA'
  return(RESULT)
}
