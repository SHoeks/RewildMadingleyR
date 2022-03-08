# set value of select_location to entire raster
SetSingleValueForAllRasters = function(sp_inputs_crop,select_location){

  sp_input_select_values = list()
  for(i in 1:length(sp_inputs)) {
    
    # get value at select_location
    sp_input_select_values[[i]] = extract(sp_inputs_crop[[i]],select_location)
    
    # if brick
    if(class(sp_inputs_crop[[i]])=="RasterBrick"){
      for(j in 1:12){
        sp_inputs_crop[[i]][[j]][] = sp_input_select_values[[i]][j]
      }
    }else{
      sp_inputs_crop[[i]][] = sp_input_select_values[[i]][1]
    }
  }
  return(sp_inputs_crop)
}