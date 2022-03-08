PlotCombinedTimelines = function(mdatas, legend=TRUE, plot=TRUE, legend_ypos=NULL, legend_xpos=NULL){

  tl_comb = mdatas[[1]]$time_line_cohorts
  tl_comb$AutotrophBiomass = mdatas[[1]]$time_line_stocks$TotalStockBiomass
  vlines = c(max(mdatas[[1]]$time_line_cohorts$Month))

  for(i in 2:length(mdatas)){

    max_year = max(tl_comb$Year)
    max_month = max(tl_comb$Month)
    vlines = c(vlines,max_month)

    mdatas[[i]]$time_line_cohorts$Year = mdatas[[i]]$time_line_cohorts$Year + max_year
    mdatas[[i]]$time_line_cohorts$Month = mdatas[[i]]$time_line_cohorts$Month + max_month

    tdf = data.frame(matrix(0,ncol=ncol(tl_comb),nrow=nrow(mdatas[[i]]$time_line_cohorts)))
    colnames(tdf) = colnames(tl_comb)

    tl_comb = rbind(tl_comb,tdf)
    new_row_idx = tl_comb$Month == 0
    tl_comb$Month[new_row_idx] = mdatas[[i]]$time_line_cohorts$Month
    tl_comb$Year[new_row_idx] = mdatas[[i]]$time_line_cohorts$Year
    tl_comb$AutotrophBiomass[new_row_idx] = mdatas[[i]]$time_line_stocks$TotalStockBiomass

    FGs = colnames(tl_comb)[grep("FG",colnames(tl_comb))]
    for(fg in FGs){
      col_idx = which(colnames(mdatas[[i]]$time_line_cohorts)==fg)
      if(length(col_idx)>0){
        tl_comb[,fg][new_row_idx] = mdatas[[i]]$time_line_cohorts[,col_idx]
      }
    }

  }

  FGs_idx = colnames(tl_comb)[grep("FG",colnames(tl_comb))]
  tl_temp = tl_comb[tl_comb==0] = NA
  min_val = floor(log10(min(tl_comb[,FGs_idx], na.rm=TRUE)))
  max_val = ceiling(log10(max(tl_comb$AutotrophBiomass, na.rm=TRUE)))
  colors=c('#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

  if(plot){
    plot(tl_comb$Month,log10(tl_comb$AutotrophBiomass), type="l", ylim=c(min_val,max_val),col="darkgreen",
         ylab="Log10 biomass",xlab="Time in months")
    c=0
    for(fg in FGs_idx){
      c=c+1
      lines(tl_comb$Month,log10(tl_comb[,fg]),col=colors[c], lty = c)
    }

    for(vv in vlines) abline(v = vv, col="grey")

    names_leg = c("Autotroph biomass",FGs_idx)
    if(legend){

      if(is.null(legend_ypos) & is.null(legend_xpos)) {
        legend(1, min_val+2, legend=names_leg,col=c("darkgreen",colors[1:c]), lty=c(1,1:c), cex=0.8)
      }else if(!is.null(legend_ypos) & is.null(legend_xpos)){
        legend(1, legend_ypos, legend=names_leg,col=c("darkgreen",colors[1:c]), lty=c(1,1:c), cex=0.8)
      }else if(!is.null(legend_ypos) & is.null(legend_xpos)){
        legend(legend_xpos, min_val+2, legend=names_leg,col=c("darkgreen",colors[1:c]), lty=c(1,1:c), cex=0.8)
      }else if(!is.null(legend_ypos) & !is.null(legend_xpos)){
        legend(legend_xpos, legend_ypos, legend=names_leg,col=c("darkgreen",colors[1:c]), lty=c(1,1:c), cex=0.8)
      }

    }

    }else{
      return(tl_comb)
    }

  }
