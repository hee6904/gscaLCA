#' @title  Teaching and Learning International Survey
#'
#' @docType data
#'
#' @usage data(TALIS)
#'
#' @format A data frame with 2560 observations on the following 6 variables.
#' \describe{
#'   \item{IDTEACH}{a numeric vector of teachers' ID.}
#'   \item{Mtv_1}{Integers with levels from 1 to 3 (1: not/low important, 2: moderate important, 3: high important); TT3G07A: To become a teacher, teaching offered a steady career path.}
#'   \item{Mtv_2}{Integers with levels from 1 to 3 (1: not/low important, 2: moderate important, 3: high important); TT3G07D: To become a teacher, teaching schedule fit with responsibilities in my personal life.}
#'   \item{Pdgg_1}{Integers with levels from 1 to 3 (1: not at all/to some extent, 2: quite a bit 3: a lot); TT3G34B: What extend you can do help my students value learning.}
#'   \item{Pdgg_2}{Integers with levels from 1 to 3 (1: not at all/to some extent, 2: quite a bit 3: a lot); TT3G34D: What extend you can do control disruptive behavior in the classroom.}
#'   \item{Stsf}{Integers with levels from 1 to 3 (1: strongly disagree/disagree, 2: agree, 3: strongly agree); TT3G53E: Feeling I enjoy working at this school.}
#' }
#'
#' @keywords datasets
#'
#' @references OECD (2019), TALIS 2018 Results (Volume I): Teachers and School Leaders as Lifelong Learners, TALIS, OECD Publishing, Paris, https://doi.org/10.1787/1d0bc92a-en.
#'
#' @source \href{http://www.oecd.org/education/talis/talis-2018-data.htm}{TALIS 2018 data}
#'
#' @examples
#' str(TALIS)
#' head(TALIS)
#'
#' @details The Teaching and Learning International Survey (TALIS) 2018 focusing on teachers, school leaders, and the learning environment in schools was conducted by the Organization for Economic Cooperation and Development (OECD). There have been three cycles, TALIS 2008, TALIS 2013, and TALIS 2018. In this study, we utilize publicly available TALIS 2018 U.S. Data, 2,560 teachers’ responses. We focused on five items: two items are on motivation, two items are on pedagogy, and the last item is on satisfaction. Items’ responses are originally four ordered categorical data. Due to too small frequencies, we modified them into three ordered categories.
#'
#'
"TALIS"
