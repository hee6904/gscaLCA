TALIS <- read.csv("data/TALIS_JSS1.csv")[,c(1,3:4,6:8)]

colnames(TALIS)[1] = c("IDTEACH")
colnames(TALIS)[2:6] = c("Mtv_1","Mtv_2",
                         "Pdgg_1","Pdgg_2" ,"Stsf")


if(!("devtools"%in%attr(installed.packages(), "dimnames")[[1]]))
  install.packages("devtools")
suppressMessages(library(devtools))

use_data(TALIS,  overwrite = T)

