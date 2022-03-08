get_binned_fw_data_abs = function(madingley_data, years_end = 5){
  
  # check if out_dir was specified manually within madingley_run()
  if(!is.null(madingley_data$out_path)){ # out_dir specified manually
    tdo = madingley_data$out_path
    # remove slashes from out_dir
    if(substr(tdo,(nchar(tdo)+1)-1,nchar(tdo))=='/')  tdo=substr(tdo,1,nchar(tdo)-1)
    if(substr(tdo,(nchar(tdo)+1)-1,nchar(tdo))=='\\') tdo=substr(tdo,1,nchar(tdo)-1)
    if(dir.exists(paste0(tdo,madingley_data$out_dir_name))) {
      out_dir = tdo
      cat(paste0("loading inputs from: ",out_dir,madingley_data$out_dir_name))
    }
  }else{ # use default output dir
    out_dir = tempdir()
    cat(paste0("loading inputs from: ",out_dir,madingley_data$out_dir_name))
  }
  
  
  # check if dir exists
  if(!dir.exists(paste0(out_dir,madingley_data$out_dir_name))){ # not exist
    stop("Unable to find output folder")
  }
  
  # paths
  # needed madingley_data!!
  if(length(list.files(paste0('/private',gsub("//", "/", out_dir, fixed=TRUE),madingley_data$out_dir_name)))==0){
    base_path = paste0(gsub("//", "/", out_dir, fixed=TRUE),madingley_data$out_dir_name)
  }else{
    base_path = paste0('/private',gsub("//", "/", out_dir, fixed=TRUE),madingley_data$out_dir_name)
  }

  # make file paths
  cohort_path = paste0(base_path,'cohort_properties')
  bins_files = list.files(cohort_path,pattern="UserDefinedBinnedCohortStatistics",full.names = T)
    
  # check if required files were exported
  if(!length(bins_files)>1){ # not exist
    stop("Required files were not exported during model run")
  }
  
  # init loop
  start_files = length(bins_files) - years_end*12
  end_files = length(bins_files)
  iterator = start_files:end_files
  bins_l = list()
  
  # loop over time steps
  for(timestep in iterator){
    
    # list files and load the last 120 csvs
    bins_files = list.files(cohort_path,pattern="UserDefinedBinnedCohortStatistics",full.names = T)
    bins = read.csv(bins_files[timestep])

    # process bins
    bins = aggregate(bins,by=list(bins$FunctionalGroupIndex,bins$LowerBodyMassBin,bins$UpperBodyMassBin),FUN=sum)
    bins = bins[bins$TotalBiomass_kg!=0,]
    bins[,c("FunctionalGroupIndex","LowerBodyMassBin","UpperBodyMassBin")] = NULL
    names(bins)[1:3] = c("FunctionalGroupIndex","LowerBodyMassBin","UpperBodyMassBin")
    
    guilds = as.vector(madingley_data$cohort_def$DEFINITION_Nutrition.source)
    guilds = guilds[bins$FunctionalGroupIndex+1]
    bins$guilds = guilds
    thermo = as.vector(madingley_data$cohort_def$DEFINITION_Endo.Ectotherm)
    thermo = thermo[bins$FunctionalGroupIndex+1]
    bins$thermo = thermo
    guilds[guilds == "Herbivore"] = 1
    guilds[guilds == "Omnivore"] = 2
    guilds[guilds == "Carnivore"] = 3
    thermo[thermo == "Endotherm"] = 0.5
    thermo[thermo == "Ectotherm"] = 0
    bins$y = as.numeric(guilds) + as.numeric(thermo)
    bins$x = bins$LowerBodyMassBin
    
    # put in list
    bins_l[[timestep]] = bins
    
  }
  
  bins = do.call(rbind,bins_l)
  bins = bins[,c("LowerBodyMassBin","UpperBodyMassBin","TotalAbundance","TotalBiomass_kg","guilds","thermo")]
  
  # average over years
  bins = aggregate(.~LowerBodyMassBin+UpperBodyMassBin+guilds+thermo, data=bins,FUN=mean)
  
  bins$ID = paste0(bins$thermo,"_",bins$guilds,"_",log10(bins$LowerBodyMassBin))
  
  return(bins)
  
}
