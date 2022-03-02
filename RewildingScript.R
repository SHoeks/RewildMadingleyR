library(raster)
library(ggplot2)
library(ggpubr)
library(plyr)
# No library(MadingleyRewilding) is needed, all function calls are with MadingleyRewilding::madingley_run()
# This is done to avoid the mix up with function of the default MadingleyR package

# source functions
source('https://raw.githubusercontent.com/SHoeks/MadingleyR_0.5degree_inputs/master/DownloadLoadHalfDegreeInputs.R')
source("https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/SetSingleValueForAllRasters.R")
source('https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/createSingleCohort.R')
source('https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/plotBinnedTimelines_v4.R')
source('https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/BinnedTimeLinesFunctions.R')
source('https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/FW_abs.r')
source("https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/get_binned_fw_data_abs4.R")
source('https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/crop_spatial_rasters_to_window.r')
source("https://raw.githubusercontent.com/SHoeks/RandomMadingleyRFunctions/master/plot_combined_timelines7.R")

# Please reset the baseDIR path to yours
baseDir="/Users/osx/Desktop/"
envDIR=paste0(baseDir,'/env/')
outDIR=paste0(baseDir,'/out_','_',format(Sys.time(), "%Y_%m_%d_%H_%M"),'/')
outDIRtemp=paste0(baseDir,'/tempdir_',sample(1000000:8000000,1),'/')

# Make dirs 
dir.create(envDIR)
dir.create(outDIR)
dir.create(outDIRtemp)

# Simulation settings
coords_LL=cbind(-5, 40) # Low NPP, Low seasonilty
coords_HL=cbind(35, 42) # High NPP, Low seasonilty
coords_LH=cbind(16, 69) # Low NPP, High seasonilty
coords_HH=cbind(32, 60) # High NPP, High seasonilty
select_location = coords_HL # select the desired location here
Years_SpinUp=50 #1000 #SpinUp
Years_PostRemoval=50 #500 #Stabilisation post removal 1
Years_PostReintroduction1=10 #50 #Stabilisation post reintroduction 1 (herbivores and omnivores)
Years_PostReintroduction2=50 #500 #Stabilisation post reintroduction 2 (carnivores)
BMinit_H=700*1000 # max mass herbivores (g) (Bison)
BMinit_O=200*1000 # max mass omnivores (g) (Bear)
BMinit_C=50*1000  # max mass carnivores (g) (Wolf)
BMafterRemoval_H=200*1000
BMafterRemoval_O=100*1000
BMafterRemoval_C=10*1000
BMinit_H_ect=1.5*1000 # max mass herbivores (g) (tortoise)
BMinit_O_ect=0.2*1000 # max mass omnivores (g) (lizard)
BMinit_C_ect=10*1000  # max mass carnivores (g) (big snake)
BMinit_H_ect_sem=0.02*1000 # max mass herbivores (g) (insects)
BMinit_O_ect_sem=0.02*1000 # max mass omnivores (g) (insects)
BMinit_C_ect_sem=0.02*1000  # max mass carnivores (g) (insects)
mxch=500 #700 #max cohort number
SampRemoved=2 # cohorts to remove per functional group per grid cell per year 
nAnimalsToReintroduce=10 #n animals to introduce per functional group
outbins = c(0.001, 0.01, 0.1, 1, 10, 20, 100, 200, 1000) # bins in KG, outputs for plotting
MnPreyDens=0.01#min density of small species (0.01 corresponds to 1/km2, which is extremely low density for small endotherms)
MinPreyDensThresh= 0.15 # 0.15 #body mass threshold below which MinPreyDensityPerHect applies (150 g seems sensible)

# Downloads and imports the 0.5 degree inputs from a zip into the selected dir
# Skips download if files already exist in selected dir
sp_inputs<-DownloadLoadHalfDegreeInputs(envDIR)

# Set spatial window
spatial_window = c(select_location[1],select_location[1]+3,select_location[2],select_location[2]+3) 

# Crop rasters, this helps to increase loading times of the spatial rasters during the model runs
sp_inputs = crop_spatial_rasters_to_window(sp_inputs, spatial_window)

# Set single value for select_location to entire raster
sp_inputs = SetSingleValueForAllRasters(sp_inputs,select_location)

# Get default (non-spatial) model params
chrt_def = MadingleyRewilding::madingley_inputs('cohort definition')
stck_def = MadingleyRewilding::madingley_inputs('stock definition')
mdl_prms = MadingleyRewilding::madingley_inputs('model parameters')

# Set max body masses
chrt_def$PROPERTY_Maximum.mass<-c(BMinit_H, BMinit_C, BMinit_O, 
                                  BMinit_H_ect_sem, BMinit_C_ect_sem, BMinit_O_ect_sem, 
                                  BMinit_H_ect, BMinit_C_ect, BMinit_O_ect)
sp_inputs$Endo_H_max[] = BMinit_H
sp_inputs$Endo_O_max[] = BMinit_O
sp_inputs$Endo_C_max[] = BMinit_C

# Create list that holds all the Madingley output data
mdata = list()

###################################
###### INITIALISATION PHASE #######
###################################

# Initialize the Madingley model
init = MadingleyRewilding::madingley_init( 
  cohort_def=chrt_def, 
  stock_def=stck_def, 
  spatial_inputs=sp_inputs, 
  spatial_window=spatial_window, 
  max_cohort = mxch)

# Run spin-up model
years = Years_SpinUp
export = years - 5
mdata$spinup = MadingleyRewilding::madingley_run( 
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
mdl_prms[51,2]<-1

# Create data needed for removal phase
mdata_tmp = mdata$spinup
CelInd<-unique(mdata_tmp$cohorts$GridcellIndex) #id cells to loop on
minBM<-min(BMafterRemoval_H, BMafterRemoval_C, BMafterRemoval_O) #min body mass among H, C and O to remove

# Set max body mass
sp_inputs$Endo_H_max[] = BMafterRemoval_H
sp_inputs$Endo_O_max[] = BMafterRemoval_O
sp_inputs$Endo_C_max[] = BMafterRemoval_C
H_Cohorts_to_remove<-mdata_tmp$cohorts[mdata_tmp$cohorts$FunctionalGroupIndex==0 & mdata_tmp$cohorts$AdultMass>BMafterRemoval_H,]
C_Cohorts_to_remove<-mdata_tmp$cohorts[mdata_tmp$cohorts$FunctionalGroupIndex==1 & mdata_tmp$cohorts$AdultMass>BMafterRemoval_C,]
O_Cohorts_to_remove<-mdata_tmp$cohorts[mdata_tmp$cohorts$FunctionalGroupIndex==2 & mdata_tmp$cohorts$AdultMass>BMafterRemoval_O,]
Cohorts_to_remove<-rbind(H_Cohorts_to_remove, C_Cohorts_to_remove, O_Cohorts_to_remove)

# This repeat goes on until all cohorts above a max body mass are gone
# it removes one random cohort above the threshold per endothermic functional group per year
tracker<-0
mdata$spinup = mdata_tmp # reset spinup
repeat {
  
  TotalRemaining=0
  
  tracker<-tracker+1
  print(paste('Year', tracker))
  
  # loop over model grid cells
  for(j in 1:length(CelInd)){
    
    cell<-CelInd[j]
    
    #for each endothermic functional group
    for (d in 0:2) { 

      Cohorts_toRemove_tmp<-mdata[[tracker]]$cohorts[mdata[[tracker]]$cohorts$GridcellIndex==cell & mdata[[tracker]]$cohorts$FunctionalGroupIndex==d & mdata[[tracker]]$cohorts$AdultMass>get(paste0('BMafterRemoval_', ifelse(d==0, 'H', ifelse(d==1, 'C', 'O')))),]
      
      Cohorts_toRemove_tmp<-Cohorts_toRemove_tmp[order(Cohorts_toRemove_tmp$AdultMass, decreasing=TRUE),]#sort by mass
      BM_toRemove<-unique(Cohorts_toRemove_tmp$AdultMass)
      
      if(length(BM_toRemove)==0){assign(paste0('BM_toRemove_', d), c()); next} #if nothing to remove, skip
      
      #remove the first X starting from the largest
      if(length(BM_toRemove)>SampRemoved) { 
        BM_toRemove2<-BM_toRemove[1:SampRemoved]
      } else {
        BM_toRemove2<-BM_toRemove #if less remaining than those removed, all are removed
      } 
      
      ind_rem<-which(mdata[[tracker]]$cohorts$AdultMass %in% BM_toRemove2 & mdata[[tracker]]$cohorts$GridcellIndex==cell & mdata[[tracker]]$cohorts$FunctionalGroupIndex==d)
      
      if(length(ind_rem)>0){
        mdata[[tracker]]$cohorts = mdata[[tracker]]$cohorts[-ind_rem,]
      }
      BM_toRemove<-BM_toRemove[!BM_toRemove %in% BM_toRemove2]
      
      assign(paste0('BM_toRemove_', d), BM_toRemove)
      
    }#closes loop through endothermic functional groups
    
    TotalRemaining = TotalRemaining + length(BM_toRemove_0) + length(BM_toRemove_1) + length(BM_toRemove_2)
    
    
  } #end loop through grid cell cells
  
  print(paste0('Large-bodied endothermic cohorts remaining: ', TotalRemaining))
  if(TotalRemaining==0) {break} #if there's nothing else to remove, break the repeat
  
  #//!!! empty outDIRtemp exept for output needed for next simulation which is stored in mdata[[tracker]]
  system(paste0("cd ",outDIRtemp, " && ls | grep -v ",sub("/", "", sub("/", "", mdata[[tracker]]$out_dir_name))," | xargs rm -rf"))
  years=1
  
  # run the model with removed cohorts for 1 year
  mdata[[tracker+1]] = MadingleyRewilding::madingley_run(
    cohort_def=chrt_def,
    stock_def=stck_def,
    spatial_inputs=sp_inputs, 
    threshold_prey_density_per_hectare=MnPreyDens, 
    threshold_prey_body_mass_kg=MinPreyDensThresh, 
    madingley_data = mdata[[tracker]],
    dispersal_off = TRUE,
    output_timestep = c(0, 999, 999, 999), 
    out_dir = outDIRtemp,  #//!
    years = years,
    max_cohort = mxch,
    model_parameters = mdl_prms,
    silenced=TRUE,
    cohort_output_bins = outbins,
    parallel = TRUE
  )
  
}

# Run another year and save the output
mdata[[length(mdata)+1]] = MadingleyRewilding::madingley_run( 
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

years = Years_PostRemoval
export = years - 5
mdata$post_removal = MadingleyRewilding::madingley_run( 
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

#REINTRODUCE HERBIVORES AND OMNIVORES and also carnivores? bad idea to reintroduce them later??!
Cohorts_to_reintroduce<-rbind(H_Cohorts_to_remove, O_Cohorts_to_remove, C_Cohorts_to_remove)

#reduce n of cohorts to reintroduce by aggregating by adult body mass, gridcell and functional group
HeaderOrder<-names(Cohorts_to_reintroduce)
Cohorts_to_reintroduce$AdultMass2<-round(Cohorts_to_reintroduce$AdultMass/1000)*1000
Cohorts_to_reintroduce<-aggregate(. ~ GridcellIndex + FunctionalGroupIndex + AdultMass2, FUN=mean, data=Cohorts_to_reintroduce)
Cohorts_to_reintroduce<-Cohorts_to_reintroduce[,HeaderOrder]
Cohorts_to_reintroduce$AdultMass2 = NULL

# change Cohorts_to_reintroduce properties
head(Cohorts_to_reintroduce)
Cohorts_to_reintroduce$IndividualReproductivePotentialMass = 0
Cohorts_to_reintroduce$MaturityTimeStep = 0
Cohorts_to_reintroduce$IsAdult = 0
Cohorts_to_reintroduce$AgeMonths = 0 
Cohorts_to_reintroduce$TimeStepsJuviline = 0
Cohorts_to_reintroduce$TimeStepsAdult = 0
Cohorts_to_reintroduce$IndividualBodyMass = (Cohorts_to_reintroduce$JuvenileMass + Cohorts_to_reintroduce$AdultMass) / 2
Cohorts_to_reintroduce$CohortAbundance = nAnimalsToReintroduce

# set the maximum allowed body masses 
sp_inputs$Endo_H_max[] = BMinit_H
sp_inputs$Endo_O_max[] = BMinit_O

# split and sort by adult bm
H_Cohorts_to_reintro = subset(Cohorts_to_reintroduce, FunctionalGroupIndex==0)
C_Cohorts_to_reintro = subset(Cohorts_to_reintroduce, FunctionalGroupIndex==1)
O_Cohorts_to_reintro = subset(Cohorts_to_reintroduce, FunctionalGroupIndex==2)
H_Cohorts_to_reintro = H_Cohorts_to_reintro[order(H_Cohorts_to_reintro$GridcellIndex, H_Cohorts_to_reintro$AdultMass, decreasing=FALSE),]
C_Cohorts_to_reintro = C_Cohorts_to_reintro[order(C_Cohorts_to_reintro$GridcellIndex, C_Cohorts_to_reintro$AdultMass, decreasing=FALSE),]
O_Cohorts_to_reintro = O_Cohorts_to_reintro[order(O_Cohorts_to_reintro$GridcellIndex, O_Cohorts_to_reintro$AdultMass, decreasing=FALSE),]

## Reintroduce herbivores+omnivores
tracker = 0
repeat{
  
  # if all cohorts are inserted break
  if(nrow(H_Cohorts_to_reintro)  + nrow(O_Cohorts_to_reintro)  == 0) {break}
  
  # if there are still herbivores or omnivores to reintroduce
  if(nrow(H_Cohorts_to_reintro)  + nrow(O_Cohorts_to_reintro) > 0) {
    
    # update tracker
    tracker = tracker + 1
    print(paste("year:",tracker))
    
    # get herbivore cohorts to insert
    if(nrow(H_Cohorts_to_reintro)>0){
      H_select = !duplicated(H_Cohorts_to_reintro$GridcellIndex) # select fist (smallest) bm cohort from each cell
      H_reintro_now = H_Cohorts_to_reintro[H_select,] # put selected in own data.frame
      H_Cohorts_to_reintro = H_Cohorts_to_reintro[!H_select,] # remove selected from cohorts to reintroduce
      H_reintro_now$BirthTimeStep = round(mean(mdata[[length(mdata)]]$cohorts$BirthTimeStep)) # update cohort properties
      print(paste("reintroducing:",nrow(H_reintro_now),"herbivore cohorts, remaining:",nrow(H_Cohorts_to_reintro)))
    }else{
      H_reintro_now = H_Cohorts_to_reintro
      print("no more herbivores to reintroduce")
    }
    
    
    # get omnivore cohorts to insert
    if(nrow(O_Cohorts_to_reintro)>0){
      O_select = !duplicated(O_Cohorts_to_reintro$GridcellIndex) # select fist (smallest) bm cohort from each cell
      O_reintro_now = O_Cohorts_to_reintro[O_select,] # put selected in own data.frame
      O_Cohorts_to_reintro = O_Cohorts_to_reintro[!O_select,] # remove selected from cohorts to reintroduce
      O_reintro_now$BirthTimeStep = round(mean(mdata[[length(mdata)]]$cohorts$BirthTimeStep)) # update cohort properties
      print(paste("reintroducing:",nrow(O_reintro_now),"omnivore cohorts, remaining:",nrow(O_Cohorts_to_reintro)))
    }else{
      O_reintro_now = O_Cohorts_to_reintro
      print("no more omnivores to reintroduce")
    }
    
    # insert cohorts
    NewCohorts = rbind(mdata[[length(mdata)]]$cohorts,H_reintro_now,O_reintro_now) 
    NewCohorts = NewCohorts[order(NewCohorts$GridcellIndex, NewCohorts$FunctionalGroupIndex, decreasing=FALSE),]
    mdata[[length(mdata)]]$cohorts = NewCohorts
  
    # run the model with reintroduced cohorts for 1 year
    print("running model for 1 year")
    mdata[[length(mdata) + 1]] = MadingleyRewilding::madingley_run(
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
years=Years_PostReintroduction1
export=years-1
mdata$betweenReintroductionPhase = MadingleyRewilding::madingley_run(
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
tracker = 0
repeat{

    # if all cohorts are inserted break
    if(nrow(C_Cohorts_to_reintro) == 0) {break}

    # If there are still carnivores to reintroduce
    if(nrow(C_Cohorts_to_reintro) > 0) {
        
        # Update tracker
        tracker = tracker + 1
        print(paste("year:",tracker))
        
        # Get carnivore cohorts to insert
        C_select = !duplicated(C_Cohorts_to_reintro$GridcellIndex) # select fist (smallest) bm cohort from each cell
        C_reintro_now = C_Cohorts_to_reintro[C_select,] # put selected in own data.frame
        C_Cohorts_to_reintro = C_Cohorts_to_reintro[!C_select,] # remove selected from cohorts to reintroduce
        C_reintro_now$BirthTimeStep = round(mean(mdata[[length(mdata)]]$cohorts$BirthTimeStep)) # update cohort properties
        print(paste("reintroducing:",nrow(C_reintro_now),"carnivore cohorts, remaining:",nrow(C_Cohorts_to_reintro)))
        
        # Insert cohorts
        NewCohorts = rbind(mdata[[length(mdata)]]$cohorts,C_reintro_now)
        NewCohorts = NewCohorts[order(NewCohorts$GridcellIndex, NewCohorts$FunctionalGroupIndex, decreasing=FALSE),]
        mdata[[length(mdata)]]$cohorts = NewCohorts
        
        # Run the model with reintroduced cohorts for 1 year
        print("running model for 1 year")
        mdata[[length(mdata)+1]] = MadingleyRewilding::madingley_run(
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
years=Years_PostReintroduction2
export=years-5
mdata$post_reintroduction = MadingleyRewilding::madingley_run(
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
plot_combined_timelines(mdata, legend_ypos = 6, legend_xpos = 0)

# Compare food-web between initial state and post reintroduction
MadingleyR::plot_foodweb(mdata$spinup)
MadingleyR::plot_foodweb(mdata$post_reintroduction)


