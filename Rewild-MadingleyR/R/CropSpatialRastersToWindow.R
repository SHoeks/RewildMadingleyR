CropSpatialRastersToWindow = function(inputs, spatial_window){

  if(class(inputs)=="list"){
    for(i in 1:length(inputs)) { inputs[[i]] = raster::crop(inputs[[i]], spatial_window) }
  }else{
    inputs = raster::crop(inputs, spatial_window)
  }

  return(inputs)
}
