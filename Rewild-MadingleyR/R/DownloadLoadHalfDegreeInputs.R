DownloadLoadHalfDegreeInputs = function(wd)
{

  if(substr(wd,(nchar(wd)+1)-1,nchar(wd))=='/')  wd=substr(wd,1,nchar(wd)-1)
  if(substr(wd,(nchar(wd)+1)-1,nchar(wd))=='\\') wd=substr(wd,1,nchar(wd)-1)

  zip_path = paste0(wd,"/0.5degree.zip")
  rasters_path = paste0(wd,"/MadingleyR_0.5degree_inputs-master")

  if(!file.exists(paste0(wd,"/0.5degree.zip"))){
      cat('Downloading zip from github repository: ')
      download.file(url = "https://github.com/SHoeks/MadingleyR_0.5degree_inputs/archive/master.zip", destfile = zip_path);
  }else{
    cat('Zip already downloaded \n')
  }
  if(!file.exists(paste0(wd,"/MadingleyR_0.5degree_inputs-master"))){
      cat('Extracting zip \n')
      setwd(wd)
      unzip('0.5degree.zip')
  }else{
    cat('Zip already extracted \n')
  }

  if (!"rgdal" %in% installed.packages()[, "Package"]) {
    stop("Package 'rgdal' not installed")
  }
  else {
    require(rgdal)
  }
  if (!"raster" %in% installed.packages()[, "Package"]) {
    stop("Package 'raster' not installed")
  }
  else {
    require(raster)
  }
  input = list()

  spatial_path = paste0(rasters_path)
  file_names = list.files(spatial_path, full.names = T, recursive = T)
  list_names = gsub("\\..*", "", list.files(spatial_path, full.names = F, recursive = T))
  FILES_sp = c("realm_classification", "land_mask", "hanpp",
               "available_water_capacity") #, "Ecto_max", "Endo_C_max", "Endo_H_max", "Endo_O_max")
  FILES_sp_temp = c("terrestrial_net_primary_productivity",
                    "near-surface_temperature", "precipitation", "ground_frost_frequency",
                    "diurnal_temperature_range")
  for (i in FILES_sp) {
    if (length(grep(i, file_names, value = T)) != 1) {
      stop("Could not find raster: ", i, ".tif \n")
    }
  }
  for (i in FILES_sp_temp) {
    if (length(grep(i, file_names, value = T)) != 12) {
      stop("Could not find raster all 12 monthly rasters containing data on: ",
           i, "\n")
    }
  }
  cat("Reading default input rasters from: ", spatial_path)
  for (i in FILES_sp) {
    file_name = grep(i, file_names, value = T)
    cat(".")
    input[[i]] = raster(file_name)
  }
  for (i in FILES_sp_temp) {
    file_name = grep(i, file_names, value = T)
    if (length(grep("_1.tif", file_name, value = T)) ==
        0) {
      file_name_sort = file_name
    }
    else {
      if (substr(spatial_path, nchar(spatial_path),
                 nchar(spatial_path)) == "/") {
        file_name_sort = paste0(spatial_path, i, "_",
                                1:12)
      }
      else {
        file_name_sort = paste0(spatial_path, "/",
                                i, "_", 1:12, ".tif")
      }
    }
    if (length(file_name_sort) == 12) {
      input[[i]] = brick(lapply(file_name_sort, raster))
    }
    cat(".")
  }
  cat("\n")

  # no max body mass raster available yet, setting single value
  input$Ecto_max = input$Endo_C_max = input$Endo_H_max = input$Endo_O_max = input$land_mask
  input$Ecto_max[] = 50000

  # Set max body masses
  input$Endo_H_max[] = 700000
  input$Endo_O_max[] = 200000
  input$Endo_C_max[] = 50000

  return(input)

}
