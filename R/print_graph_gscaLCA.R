print_graph_gscaLCA= function(Iden.vect, LEVELs, P)
{
  if(all(Iden.vect==TRUE)){
    if(length(LEVELs[[1]])==2){
      print(P[[1]])
    }else{

      get_legend <- function(myggplot){
        tmp <- ggplot_gtable(ggplot_build(myggplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }

      legend <- get_legend(P[[1]]+theme(legend.position="bottom"))

      P.1 =lapply(P, function(x) x + theme(legend.position="none"))
      P.1[[length(LEVELs[[1]])+1]] =legend


      grid.arrange(grobs = P.1, layout_matrix = matrix(
        c(rep(1:length(LEVELs[[1]]),each=5),
          length(LEVELs)+1),ncol=1))


      #do.call(grid.arrange, c(P.1, nrow= length(LEVELs[[1]])))
    }
  }
}
