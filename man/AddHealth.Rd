\name{AddHealth}
\alias{AddHealth}
\docType{data}
\title{
Add Health data about substance use and abuse
}
\description{
%5,114 observation with 6 variables including ID. The five variables are binary.
}
\usage{data("Addhealth")}
\format{
  A data frame with 5114 observations on the following 6 variables.
  \describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{Smoking}}{a factor with levels \code{Yes} or \code{No}; H4TO1: Have you ever smoked an entire cigarette? }
    \item{\code{Alcohol}}{a factor with levels \code{Yes} or \code{No}; H4TO33: Have you had a drink of beer, wine, or liquor more than two or three times? Do not include sips or tastes from someone else’s drink.}
    \item{\code{Drug}}{a factor with levels \code{Yes} or \code{No}; H4TO65E: Have you every used any of the following drugs? (65E) Other types of illegal drugs, such as LSD, PCP, ecstasy, heroin, or mushrooms; or inhalants. }
    \item{\code{Marijuana}}{a factor with levels \code{Yes} or \code{No}; H4TO65B: Have you every used any of the following drugs? (65B) Marijuana (hash, bhang, ganja) }
    \item{\code{Cocaine}}{a factor with levels \code{Yes} or \code{No}; H4TO65C: Have you every used any of the following drugs? (65C) Cocaine (crack, coca leaves) }
  }
}
\details{
This AddHealth data consist of 5,144 participants' responses with a randomly generated ID variable and five item variables, such as Smoking, Alcohol, Other Types of Illegal Drug, Marijuana, and Cocaine. The responses of the five items are dichotomous as either ``Yes'' or ``No'' and are treated the other missing codes as systematic missing. This data can be obtained from the National Longitudinal Study of Adolescent to Adult Health (Add Health; Harris et al., 2009) where the study have mainly focused on the investigation of how health factors in childhood affect adult outcomes. In terms of data collection, there have been four additional waves since 1994. In this package, the data of a specific section of substance use at the wave IV is pre-installed.
}

\source{
https://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/21600/datadocumentation#
}
\references{
%%Harris, KM (2009) The national longitudinal study of adolescent to adult health (Add Health), Waves I &II, 1994-1996; Wave III, 2001-2002; Wave IV, 2007-2009 (Machine-readable data file and documentation). Chapel Hill: Carolina Population Center, University of North Carolina at Chapel Hill. Retrieved from http://doi.org/10.3886/ICPSR21600.v21

Harris, Kathleen Mullan, and Udry, J. Richard. National Longitudinal Study of Adolescent to Adult Health (Add Health), 1994-2008 [Public Use]. Ann Arbor, MI: Carolina Population Center, University of North Carolina-Chapel Hill [distributor], Inter-university Consortium for Political and Social Research [distributor], 2018-08-06. https://doi.org/10.3886/ICPSR21600.v21

}
\examples{
str(AddHealth)
head(AddHealth)
}
\keyword{datasets}
