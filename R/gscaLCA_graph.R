gscaLCA_graph= function( RespProb.results, LEVELs, LCprevalence.result)
{

  varnames = names(RespProb.results)


  P = list()
  for (j in 1:length(LEVELs[[1]]))
  {
    matYes = c()
    for(l in 1:length(RespProb.results))
    {

      matYes = rbind(matYes,
                     subset(RespProb.results[[l]],  RespProb.results[[l]]$Category==LEVELs[[l]][j])[,1:3])
    }

    matYes = cbind(rep(varnames, each = sum(LCprevalence.result[,1]!=0)), matYes)
    matYes$Category=NULL
    names(matYes)[1] ="Type"
    matYes$Type = factor(matYes$Type, levels= varnames, labels=varnames)


    class.numeric = as.numeric(str_extract(unique(matYes$Class), "\\-*\\d+\\.*\\d*"))


    matYes$Class = factor(matYes$Class , unique(matYes$Class),
                          paste0(class.numeric,
                                 " (",sprintf("%.2f", LCprevalence.result[class.numeric,"Percent"]), "%)") )

    # Class= matYes$Class


    ## Line plot by using Aggregated Data (with or without error bar)

    P[[j]]= ggplot(matYes, aes(x= matYes$Type , y= matYes$Estimate, colour= matYes$Class, group=matYes$Class)) +
      geom_line(size=1, aes(linetype=matYes$Class)) +
      geom_point(size=3, aes(shape =matYes$Class))+
      theme_light()+
      ylim(0, 1)+
      #ylab("Probability of Yes")+
      ggtitle(paste("Response:", LEVELs[[1]][j])) +
      theme(
        plot.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.title.y = element_text(size =13),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.text.y = element_text(size = 15))+
      labs(color='Class', shape="Class", linetype="Class")

  }

}


