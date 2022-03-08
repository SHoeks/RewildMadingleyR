
plotBinnedTimelines = 
  function(madingley_data,plot=TRUE,only_above_zero=FALSE){
    
    cdef = madingley_data$cohort_def
    
    if(class(madingley_data$out_path)=="NULL"){
      dir_path = paste0(tempdir(),madingley_data$out_dir_name,"/cohort_properties/") # make dir name
    }else{
      dir_path = paste0(madingley_data$out_path,madingley_data$out_dir_name,"/cohort_properties/") # make dir name
    }
    
    out_paths = list.files(dir_path,pattern='.csv',full.names = T) # list cohort output files
    out_paths = grep("UserDefinedBinnedCohortStatistics",out_paths,value = T) # grep only BinnedCohortStatistics
    data = lapply(out_paths,read.csv) # read all data
    
    print(dir_path)
    
    for(i in 1:length(data)){ 
      
      data[[i]]$Month = i # add months
      data[[i]]$Group = 1:nrow(data[[i]]) # add group id
      
    }
    
    data = do.call("rbind", data) # rbind all data in list
    n_groups = max(data$Group) # get number of groups
    n_months = max(data$Month) # total number of monhts
    max_fg = max(data$FunctionalGroupIndex) # fg end
    min_fg = min(data$FunctionalGroupIndex) # fg start
    
    output_data = list()
    
    # plot data
    for(fg in min_fg:max_fg){
      
      data_fg = data[data$FunctionalGroupIndex == fg,]
      par(mfrow=c(2,4))
      
      for(bin in unique(data_fg$Group)){
        
        data_fg_bin = data_fg[data_fg$Group == bin,]
        TotalBiomass_kg = log10(data_fg_bin$TotalBiomass_kg)
        TotalBiomass_kg = ifelse(is.infinite(TotalBiomass_kg),0,TotalBiomass_kg)
        
        if(plot){
          if(min_fg==0) {fg_data = fg + 1} else {fg_data = fg}
          diet = cdef$DEFINITION_Nutrition.source[fg_data]
          therm = cdef$DEFINITION_Endo.Ectotherm[fg_data]
          main1 = paste(diet," ",therm," (",fg,")",'\n bin = ',
                        data_fg_bin$LowerBodyMassBin[1],' to ',
                        data_fg_bin$UpperBodyMassBin[1],'kg',sep='')
        }
        if(plot) {
          if(only_above_zero){
            if( sum(TotalBiomass_kg)>0){
              plot(data_fg_bin$Month,TotalBiomass_kg,type='l',main=main1,xlab="month",ylab='log10 biomass kg',col = 'darkgrey')
            }
          }else{
            plot(data_fg_bin$Month,TotalBiomass_kg,type='l',main=main1,xlab="month",ylab='log10 biomass kg',col = 'darkgrey')
          }
        }
      }
      
      if(!plot){
        if(min_fg==0) {fg_data = fg + 1} else {fg_data = fg}
        data_fg$Group = NULL
        data_fg$LowerLog10Bin = log10(data_fg$LowerBodyMassBin)
        data_fg$LowerBodyMassBin = NULL
        data_fg$UpperLog10Bin = log10(data_fg$UpperBodyMassBin)
        data_fg$UpperBodyMassBin = NULL
        output_data[[fg_data]] = data_fg
      }
    }
    
    par(mfrow=c(1,1))
    
    if(!plot){
      return(output_data)
    }
    
    
  }