# Required library
library(tcltk) # Progress bar
library(RHRV)  # HRV indicators

# AnalyseTemporel 
# Compute HRV indicators for a 60 seconds sample with RHRV package.
# Size is on of the parameter used in the computation of the indicators
AnalyseTemporel = function(ligne, size, pb)
{

  hrv.data = CreateHRVData(FALSE)
  hrv.data$Beat = list(
    Time = cumsum(ligne),
    RR = ligne
  ) 
  
  # Computation of the indicators
  hrv.data = CreateTimeAnalysis(hrv.data,
                                size = size)
  
  indicateur = as.numeric(unlist(hrv.data$TimeAnalysis)[-1])
  return(indicateur)
}


#######


# IndicateursTemp returns the indicators for every sample of 60 seconds in the list "data".
# size is the length of the window used in order to compute some of the indicators with RHRV
IndicateursTemp = function(data, size)
{
  cat("Calcul des indicateurs de la dimension temporelle pour une fenêtre de",size,"seconde(s)\n")
  # Parameters for the progress bar
  iteration <<- 0
  total <<- length(data)
  pb <- tkProgressBar("Analyse Temporelle", min = 0, max = length(data), initial = 0)
  
  # Computation of the indicators with AnalyseTemporel function for every sample in data
  T1 = Sys.time()
  res = data.frame(matrix(unlist(lapply(lapply(data, unlist), FUN = function(x){
    iteration <<- iteration+1
    val = as.character(round(iteration/total,4)*100)
    info <- sprintf("%s%% done", val)
    setTkProgressBar(pb, iteration, sprintf("Calcul pour la fenêtre de taille %d",size,"secondes"),info)
    ind.temps = AnalyseTemporel(x, size = size)
    return(ind.temps)
  })),
  ncol = 10,
  nrow = length(data),
  byrow = TRUE))
  T2 = Sys.time()
  # The indicators are returned in a data.frame
  
  # Names of the indicators
  colnames(res) = c("SDNN","SDANN","SDNNIDX",
                    "pNN50","SDSD","RMSSD",
                    "IRRR","MADRR","TINN","HRV.index")
  
  # Execution time
  Temps = as.numeric(difftime(T2,T1, units = "secs"))
  cat("Temps total d'exécution :",Temps%/%60,"minute(s) et ",Temps%%60,"seconde(s)\n")
  cat("\n")
  
  # Renvoie la taille de la fenêtre considérée, les valeurs des indicateurs de chaque période (60sec)
  # et le temps d'exécution pour l'ensemble des lignes
  
  # Give back the size parameter, values of the indicators (in a data.frame) and the execution time
  close(pb)
  return(list(Size = size, Resultats = res, Temps = Temps))
}

# Return the indicators but for a vector of size
IndicateursTempWind = function(data, size)
{
    return(lapply(size, FUN = function(x){IndicateursTemp(data, x)}))
}
