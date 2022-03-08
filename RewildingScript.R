library(MadingleyRewilding) # load to run the altered version of the package
library(MadingleyR) # load to make plots at the end of simulation
library(raster) # for spatial inputs
library(dplyr) # for data manipulation

# Please reset the baseDIR path to a folder of your choice
baseDir <- "/Users/osx/Desktop/"
envDIR <- paste0(baseDir,'/env/')
outDIR <- paste0(baseDir,'/out_','_',format(Sys.time(), "%Y_%m_%d_%H_%M"),'/')
outDIRtemp <- paste0(baseDir,'/tempdir_',sample(1000000:8000000,1),'/')
dir.create(envDIR)
dir.create(outDIR)
dir.create(outDIRtemp)

# Simulation settings
coords_LL <- cbind(-5, 40) # Low NPP, Low seasonilty
coords_HL <- cbind(35, 42) # High NPP, Low seasonilty
coords_LH <- cbind(16, 69) # Low NPP, High seasonilty
coords_HH <- cbind(32, 60) # High NPP, High seasonilty
select_location <- coords_HL # select the desired location here
spatial_window <- c(select_location[1],select_location[1]+3,select_location[2],select_location[2]+3) # Set spatial window
Years_SpinUp <- 50 #1000 #SpinUp
Years_PostRemoval <- 50 #500 #Stabilisation post removal 1
Years_PostReintroduction1 <- 10 #50 #Stabilisation post reintroduction 1 (herbivores and omnivores)
Years_PostReintroduction2 <- 50 #500 #Stabilisation post reintroduction 2 (carnivores)
BMafterRemoval_H <- 200*1000 # max mass herbivores after removal (g) (Bison)
BMafterRemoval_O <- 100*1000 # max mass omnivores after removal (g) (Bear)
BMafterRemoval_C <- 10*1000  # max mass carnivores after removal (g) (Wolf)
mxch <- 500 #700 #max cohort number
SampRemoved <- 2 # cohorts to remove per functional group per grid cell per year 
nAnimalsToReintroduce <- 10 #n animals to introduce per functional group
outbins <- c(0.001, 0.01, 0.1, 1, 10, 20, 100, 200, 1000) # bins in KG, outputs for plotting
MnPreyDens <- 0.01#min density of small species (0.01 corresponds to 1/km2, which is extremely low density for small endotherms)
MinPreyDensThresh <- 0.15 # 0.15 #body mass threshold below which MinPreyDensityPerHect applies (150 g seems sensible)

# Downloads and imports the 0.5 degree inputs from a zip into the selected dir
# Skips download if files already exist in selected dir
# After downloading, rasters get cropped to spatial window 
# This helps to increase loading times of the spatial rasters during the model runs
sp_inputs <- DownloadLoadHalfDegreeInputs(envDIR) %>% 
  CropSpatialRastersToWindow(spatial_window) %>%
  SetSingleValueForAllRasters(select_location)

# Get Rewilding Europe specific cohort definitions
chrt_def <- GetRewildingEuropeCohortDefs()

# Get other default (non-spatial) model params
stck_def <- MadingleyInputs('stock definition')
mdl_prms <- MadingleyInputs('model parameters')

###################################
###### INITIALISATION PHASE #######
###################################

# Create list that holds all the Madingley output data
mdata <- list()

# Initialize the Madingley model
init <- MadingleyInit( 
  cohort_def=chrt_def, 
  stock_def=stck_def, 
  spatial_inputs=sp_inputs, 
  spatial_window=spatial_window, 
  max_cohort = mxch)

# Run spin-up model
years <- Years_SpinUp
export <- years - 5
mdata$spinup <- MadingleyRun( 
  madingley_data = init,
  cohort_def=chrt_def, 
  stock_def=stck_def, 
  spatial_inputs=sp_inputs, 
  threshold_prey_density_per_hectare=MnPreyDens, 
  threshold_prey_body_mass_kg=MinPreyDensThresh, 
  out_dir = outDIR, 
  years = years, 
  max_cohort = mxch, 
  output_timestep = c(0,export,export,export),
  model_parameters = mdl_prms,
  silenced=FALSE, parallel = TRUE,
  cohort_output_bins = outbins
) 

###################################
######## REMOVAL PHASE ###########
###################################

# Stop cohorts from "evolving" (changing their body mass between generations) 
# Settings this to 1 fixes their body mass, this helps in the removal process
# The juvenile body mass is always identical to the parent cohort
mdl_prms[51,2] <- 1

# Create data needed for removal phase
mdata_tmp <- mdata$spinup

# Set max body mass
sp_inputs$Endo_H_max[] <- BMafterRemoval_H
sp_inputs$Endo_C_max[] <- BMafterRemoval_C
sp_inputs$Endo_O_max[] <- BMafterRemoval_O

# Get cohorts to remove, used for the reintroduction later
Cohorts_to_remove <- mdata_tmp$cohorts %>% 
  filter(FunctionalGroupIndex == 0 & AdultMass>BMafterRemoval_H | 
         FunctionalGroupIndex == 1 & AdultMass>BMafterRemoval_C | 
         FunctionalGroupIndex == 2 & AdultMass>BMafterRemoval_O )

# This repeat goes on until all cohorts above a max body mass are gone
# It removes one random cohort above the threshold per endothermic functional group per year
tracker <- 0
mdata$spinup <- mdata_tmp # reset spinup
repeat {
  
  TotalRemaining <- 0
  tracker <- tracker+1
  print(paste('Year', tracker))
  
  # loop over model grid cells
  for(j in 1:length(unique(mdata_tmp$cohorts$GridcellIndex))){
    
    # Get cell index
    cell <- unique(mdata_tmp$cohorts$GridcellIndex)[j]
    
    #for each endothermic functional group
    for (d in 0:2) { 
      
      # functional group index to feeding guild category
      guild <- switch(as.character(d),"0"="H","1"="C","2"="O")
      
      # Get unique body masses of cohorts to remove
      BM_toRemove <- mdata[[tracker]]$cohorts %>% 
        filter(GridcellIndex==cell & FunctionalGroupIndex==d) %>% 
        filter(AdultMass>get(paste0('BMafterRemoval_',guild))) %>% 
        arrange(desc(AdultMass)) %>% 
        pull(AdultMass) %>%
        unique() 
      
      # If nothing to remove, skip
      if(length(BM_toRemove)==0){
        assign(paste0('BM_toRemove_', d), c()); next
      } 
      
      # Remove the first X starting from the largest
      if(length(BM_toRemove)>SampRemoved) { 
        BM_toRemove2 <- BM_toRemove[1:SampRemoved]
      } else {
        BM_toRemove2 <- BM_toRemove
      } 

      ind_rem <- mdata[[tracker]]$cohorts %>% 
        with(which(AdultMass %in% BM_toRemove2 & GridcellIndex==cell & FunctionalGroupIndex==d))
      
      if(length(ind_rem)>0) {
        mdata[[tracker]]$cohorts <- mdata[[tracker]]$cohorts[-ind_rem,]
      }
  
      BM_toRemove <- BM_toRemove[!BM_toRemove %in% BM_toRemove2]
      
      assign(paste0('BM_toRemove_', d), BM_toRemove)
      
    }# Closes loop through endothermic functional groups
    
    TotalRemaining <- TotalRemaining + length(BM_toRemove_0) + length(BM_toRemove_1) + length(BM_toRemove_2)
    
    
  } # End loop through grid cell cells
  
  print(paste0('Large-bodied endothermic cohorts remaining: ', TotalRemaining))
  
  # If there's nothing else to remove, break the repeat
  if(TotalRemaining==0) {
    break
  } 
  
  # Empty outDIRtemp except for output needed for next simulation which is stored in mdata[[tracker]]
  system(paste0("cd ",outDIRtemp, " && ls | grep -v ",sub("/", "", sub("/", "", mdata[[tracker]]$out_dir_name))," | xargs rm -rf"))

  # Run the model with removed cohorts for 1 year
  mdata[[tracker+1]] <- MadingleyRun(
    cohort_def=chrt_def,
    stock_def=stck_def,
    spatial_inputs=sp_inputs, 
    threshold_prey_density_per_hectare=MnPreyDens, 
    threshold_prey_body_mass_kg=MinPreyDensThresh, 
    madingley_data = mdata[[tracker]],
    dispersal_off = TRUE,
    output_timestep = c(0, 999, 999, 999), 
    out_dir = outDIRtemp,
    years = 1,
    max_cohort = mxch,
    model_parameters = mdl_prms,
    silenced=TRUE,
    cohort_output_bins = outbins,
    parallel = TRUE
  )
  
}

# Run another year and save the output
mdata[[length(mdata)+1]] <- MadingleyRun( 
  madingley_data = mdata[[length(mdata)]],
  cohort_def=chrt_def, 
  stock_def=stck_def,
  spatial_inputs=sp_inputs,
  threshold_prey_density_per_hectare=MnPreyDens, 
  threshold_prey_body_mass_kg=MinPreyDensThresh,  
  out_dir = outDIR, 
  years = 1,
  max_cohort = mxch, 
  model_parameters = mdl_prms, 
  cohort_output_bins = outbins,
  parallel = TRUE
)

###################################
###### STABILIZATION PHASE ########
###################################

years <- Years_PostRemoval
export <- years - 5
mdata$post_removal <- MadingleyRun( 
  madingley_data = mdata[[length(mdata)]],
  cohort_def=chrt_def, stock_def=stck_def,
  spatial_inputs=sp_inputs, 
  threshold_prey_density_per_hectare=MnPreyDens, 
  threshold_prey_body_mass_kg=MinPreyDensThresh,   
  out_dir = outDIR, 
  years = years,
  max_cohort = mxch, 
  model_parameters = mdl_prms, 
  output_timestep = c(0,export,export,export),
  parallel = TRUE,
  cohort_output_bins = outbins
)

###################################
###### REINTRODUCTION PHASE #######
###################################

# Get cohort properties to reintroduce
# This function reduces the number of cohorts to reintroduce
# It does so by aggregating by adult body mass, gridcell and functional group
Cohorts_to_reintroduce <- GetCohortsToReintroduce(Cohorts_to_remove)

# set the maximum allowed body masses 
#sp_inputs$Endo_H_max[] = chrt_def$PROPERTY_Maximum.mass[chrt_def$NOTES_group.description=="Bison"]
#sp_inputs$Endo_O_max[] = chrt_def$PROPERTY_Maximum.mass[chrt_def$NOTES_group.description=="Bear"]

# Split cohorts to reintroduce and sort by adult body mass
H_Cohorts_to_reintro <- Cohorts_to_reintroduce %>% 
  filter(FunctionalGroupIndex==0) %>% 
  arrange(GridcellIndex, AdultMass)
C_Cohorts_to_reintro <- Cohorts_to_reintroduce %>% 
  filter(FunctionalGroupIndex==1) %>% 
  arrange(GridcellIndex, AdultMass)
O_Cohorts_to_reintro <- Cohorts_to_reintroduce %>% 
  filter(FunctionalGroupIndex==2) %>% 
  arrange(GridcellIndex, AdultMass)


## Reintroduce herbivores+omnivores
tracker <- 0
repeat{
  
  # if all cohorts are inserted break
  if(nrow(H_Cohorts_to_reintro)  + nrow(O_Cohorts_to_reintro)  == 0) {break}
  
  # if there are still herbivores or omnivores to reintroduce
  if(nrow(H_Cohorts_to_reintro)  + nrow(O_Cohorts_to_reintro) > 0) {
    
    # update tracker
    tracker <- tracker + 1
    print(paste("year:",tracker))
    
    # get herbivore cohorts to insert
    if(nrow(H_Cohorts_to_reintro)>0){
      H_select <- !duplicated(H_Cohorts_to_reintro$GridcellIndex) # select fist (smallest) bm cohort from each cell
      H_reintro_now <- H_Cohorts_to_reintro[H_select,] # put selected in own data.frame
      H_Cohorts_to_reintro <- H_Cohorts_to_reintro[!H_select,] # remove selected from cohorts to reintroduce
      H_reintro_now$BirthTimeStep <- round(mean(mdata[[length(mdata)]]$cohorts$BirthTimeStep)) # update cohort properties
      print(paste("reintroducing:",nrow(H_reintro_now),"herbivore cohorts, remaining:",nrow(H_Cohorts_to_reintro)))
    }else{
      H_reintro_now <- H_Cohorts_to_reintro
      print("no more herbivores to reintroduce")
    }
    
    
    # get omnivore cohorts to insert
    if(nrow(O_Cohorts_to_reintro)>0){
      O_select <- !duplicated(O_Cohorts_to_reintro$GridcellIndex) # select fist (smallest) bm cohort from each cell
      O_reintro_now <- O_Cohorts_to_reintro[O_select,] # put selected in own data.frame
      O_Cohorts_to_reintro <- O_Cohorts_to_reintro[!O_select,] # remove selected from cohorts to reintroduce
      O_reintro_now$BirthTimeStep <- round(mean(mdata[[length(mdata)]]$cohorts$BirthTimeStep)) # update cohort properties
      print(paste("reintroducing:",nrow(O_reintro_now),"omnivore cohorts, remaining:",nrow(O_Cohorts_to_reintro)))
    }else{
      O_reintro_now <- O_Cohorts_to_reintro
      print("no more omnivores to reintroduce")
    }
    
    # insert cohorts
    NewCohorts <- rbind(mdata[[length(mdata)]]$cohorts,H_reintro_now,O_reintro_now) 
    NewCohorts <- NewCohorts[order(NewCohorts$GridcellIndex, NewCohorts$FunctionalGroupIndex, decreasing=FALSE),]
    mdata[[length(mdata)]]$cohorts = NewCohorts
  
    # run the model with reintroduced cohorts for 1 year
    print("running model for 1 year")
    mdata[[length(mdata) + 1]] <- MadingleyRun(
      madingley_data = mdata[[length(mdata)]],
      cohort_def=chrt_def, 
      stock_def=stck_def, 
      spatial_inputs=sp_inputs, 
      threshold_prey_density_per_hectare=MnPreyDens, 
      threshold_prey_body_mass_kg=MinPreyDensThresh, 
      output_timestep = c(0, 999, 999, 999), 
      out_dir = outDIRtemp, years = 1, parallel = TRUE,
      max_cohort = mxch, model_parameters = mdl_prms, silenced=TRUE, cohort_output_bins = outbins)

  }
}


# Stabilize before reintroducing carnivores
years <- Years_PostReintroduction1
export <- years-1
mdata$betweenReintroductionPhase <- MadingleyRun(
    madingley_data = mdata[[length(mdata)]],
    cohort_def=chrt_def, stock_def=stck_def,
    spatial_inputs=sp_inputs,
    threshold_prey_density_per_hectare=MnPreyDens, 
    threshold_prey_body_mass_kg=MinPreyDensThresh,
    output_timestep = c(0, export,export,export),
    out_dir = outDIR,
    years = years,
    max_cohort = mxch,
    model_parameters = mdl_prms,
    silenced=FALSE,
    cohort_output_bins = outbins,
    parallel = TRUE
)

## Reintroduce carnivores
#reintroduction_m_data2 = list() # list to store the m_data in created over the removal period
#reintroduction_m_data2[[1]] = m_data4 # put m_data2 in list as starting point for the removal, will be overwritten at first iteration
tracker <- 0
repeat{

    # if all cohorts are inserted break
    if(nrow(C_Cohorts_to_reintro) == 0) {break}

    # If there are still carnivores to reintroduce
    if(nrow(C_Cohorts_to_reintro) > 0) {
        
        # Update tracker
        tracker <- tracker + 1
        print(paste("year:",tracker))
        
        # Get carnivore cohorts to insert
        C_select <- !duplicated(C_Cohorts_to_reintro$GridcellIndex) # select fist (smallest) bm cohort from each cell
        C_reintro_now <- C_Cohorts_to_reintro[C_select,] # put selected in own data.frame
        C_Cohorts_to_reintro <- C_Cohorts_to_reintro[!C_select,] # remove selected from cohorts to reintroduce
        C_reintro_now$BirthTimeStep <- round(mean(mdata[[length(mdata)]]$cohorts$BirthTimeStep)) # update cohort properties
        print(paste("reintroducing:",nrow(C_reintro_now),"carnivore cohorts, remaining:",nrow(C_Cohorts_to_reintro)))
        
        # Insert cohorts
        NewCohorts <- rbind(mdata[[length(mdata)]]$cohorts,C_reintro_now)
        NewCohorts <- NewCohorts[order(NewCohorts$GridcellIndex, NewCohorts$FunctionalGroupIndex, decreasing=FALSE),]
        mdata[[length(mdata)]]$cohorts <- NewCohorts
        
        # Run the model with reintroduced cohorts for 1 year
        print("running model for 1 year")
        mdata[[length(mdata)+1]] <- MadingleyRun(
            madingley_data = mdata[[length(mdata)]],
            cohort_def=chrt_def, stock_def=stck_def, spatial_inputs=sp_inputs,
            threshold_prey_density_per_hectare=MnPreyDens, 
            threshold_prey_body_mass_kg=MinPreyDensThresh,
            parallel = TRUE,
            output_timestep = c(0, 999, 999, 999),
            out_dir = outDIRtemp, years = 1,
            max_cohort = mxch, model_parameters = mdl_prms, silenced=TRUE, cohort_output_bins = outbins)
        
    }

}

# Stabilize after carnivore reintroduction
years <- Years_PostReintroduction2
export <- years-5
mdata$post_reintroduction <- MadingleyRun(
  madingley_data = mdata[[length(mdata)]],
  cohort_def=chrt_def, stock_def=stck_def, 
  spatial_inputs=sp_inputs, 
  threshold_prey_density_per_hectare=MnPreyDens, 
  threshold_prey_body_mass_kg=MinPreyDensThresh,   
  output_timestep = c(0, export,export,export),
  out_dir = outDIR,
  years = years,
  max_cohort = mxch,
  model_parameters = mdl_prms,
  parallel = TRUE,
  cohort_output_bins = outbins)

# Empty the temporary folder: outDIRtemp
system(paste0('rm -rf ',outDIRtemp,'*')) 

# Plot entire timeline
par(mfrow=c(1,1))
PlotCombinedTimelines(mdata, legend_ypos = 6, legend_xpos = 0)

# Compare food-web between initial state and post reintroduction
plot_foodweb(mdata$spinup)
plot_foodweb(mdata$post_reintroduction)


