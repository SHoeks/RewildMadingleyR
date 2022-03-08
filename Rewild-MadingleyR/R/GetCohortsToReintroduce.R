GetCohortsToReintroduce = function(Cohorts_to_reintroduce){

  #reduce n of cohorts to reintroduce by aggregating by adult body mass, gridcell and functional group
  HeaderOrder = names(Cohorts_to_reintroduce)
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

  return(Cohorts_to_reintroduce)

}
