#' Summary of gscaLCA output or gscaLCR output
#'
#' @param object the object of gscaLCA or gscaLCR
#' @param print.cov.output a character of what type partitioning and regression. Four possible option are possible "multinomial.hard", "multinomial.soft", "binomial.hard", and "binomial.soft".
#' @param ... Additional arguments affecting the summary produced.
#'
#' @return print model fit, prevalence, item probabilities, and regression results
#' @export
#'
#' @examples
#' # summary(R2)
summary.gscaLCA = function(object, print.cov.output = NULL, ...)
{
  if(class(object) !="gscaLCA") stop ("The object is not from gscaLCA")

  if(is.null(print.cov.output)){
    print_gscaLCA(object$num.class,object$N, object$N.origin, object$Boot.num,
                  object$Boot.num.im,
                  object$model.fit,
                  object$LCprevalence,
                  object$RespProb,
                  cov_results.multi.hard = NULL,
                  cov_results.bin.hard = NULL,
                  cov_results.multi.soft = NULL,
                  cov_results.bin.soft = NULL,
                  print.cov.output = NULL)
  }else{
    print_gscaLCA(object$num.class,object$N, object$N.origin, object$Boot.num,
                  object$Boot.num.im,
                  object$model.fit,
                  object$LCprevalence,
                  object$RespProb,
                  cov_results.multi.hard = object$cov_results.multi.hard,
                  cov_results.bin.hard = object$cov_results.bin.hard,
                  cov_results.multi.soft = object$cov_results.multi.soft,
                  cov_results.bin.soft = object$cov_results.bin.soft,
                  print.cov.output = print.cov.output)
  }

  print_graph_gscaLCA (object$all.Levels.equal, object$LEVELs, object$plot)
}

