AddHealth <- read.csv("data/TAD_DATA_062419.csv")

for (i in 2:6)
{
  AddHealth[,i]= ifelse(AddHealth[,i]==1, "Yes",
                      ifelse(AddHealth[,i]==0, "No",NA))

  AddHealth[,i] = factor(AddHealth[,i], levels = c("Yes", "No"))
}

rownames(AddHealth) = 1:nrow(AddHealth)

if(!("devtools"%in%attr(installed.packages(), "dimnames")[[1]]))
  install.packages("devtools")
suppressMessages(library(devtools))

use_data(AddHealth,  overwrite = T)

