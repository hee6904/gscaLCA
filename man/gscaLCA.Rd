\name{gscaLCA}
\alias{gscaLCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function of LCA by using fuzzy clustering GSCA
}
\description{
This function enables to run LCA based on GSCA algorithm.
}
\usage{
gscaLCA(dat, varnames, ID.var=NULL, num.cluster=2, num.factor="EACH", Boot.num=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  Data that you want to implement GSCA.
  }
  \item{varnames}{
  A character vector. The names of columns to be used for gscaLCA.
  }
 \item{ID.var}{
  A character element. The name of ID variable. If ID variable is not specified, gscaLCA will find a ID variable in a given data. The ID of observation are automatically created when data set does not have any ID variable. The default is NULL.
  }
  \item{num.cluster}{
  A numeric element. Number of cluster to be analyzed.
  }
  \item{num.factor}{
  Either "EACH" or "ALLin1".
  "EACH" indicates that each variable assumes to have latent variable.
  "ALLin1" indicates that all variables assumes to share one latent variable. Default is "EACH".
  }
  \item{Boot.num}{
  Number of bootstrap. The standard errors of parameters are obtained by bootstrap in GSCA algorithm. Default is 20.
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{N }{The number of observations used after listwise delection when missing values exist.}
 \item{C }{The number of cluster.}
 \item{Boot.num.im}{The number of bootstrap implemented. Not all iterations of bootstrap is appropriate to estimate the standard error.}
  \item{model.fit}{The model fit indices. FIT, AFIT, FPI and NCE are provided with the standard error and 95\% credible interval lower and upper bounds.}
  \item{LCprevalence}{The latent class prevalence. The percent of class, the number of observation for each class, standard error, and 95\% credible interval lower and upper bounds are provided.}
  \item{RespRrob}{The item response probabilities for each variable are reported as elements of a list. Each element consists of a table containing the probabilities with respect to the possible categories of each variable. The standard error and 95\% credible interval lower and upper bounds are also reported.}
  \item{membership}{A data frame of the posterior probabilities for each examinee with the predicted class membership.}
  \item{plot}{Graphs for each categories. For example, with two categories, each graph is stored as "p1" and "p2" in the list of "plot". When the number of options across variables are different, the graphs are not provided.}
%% ...
}
\references{
Ryoo, J. H., Park, S., & Kim, S. (2019). Categorical latent variable modeling utilizing fuzzy clustering generalized structured component analysis as an alternative to latent class analysis. Behaviormetrika, 1-16.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# AddHealth data with 2 clusters
R2 = gscaLCA(AddHealth, varnames = names(AddHealth)[2:6], num.cluster = 2)
R2$model.fit      # Model fit
R2$LCprevalence   # Latent Class Prevalence
R2$RespProb       # Item Reponse Probability
R2$membership     # Membership for all observations

# TALIS data with 3 clusters
T3 = gscaLCA(TALIS, varnames = names(TALIS)[2:6], num.factor = "ALLin1", num.cluster = 3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
