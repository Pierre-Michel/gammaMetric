##### Preambule ######

# Required libraries
install.packages("R.matlab")        # Allow us to read .mat files
install.packages("ggplot2")         # Graphics
install.packages("randomForest")    # Compute a random forest model
install.packages("pROC")            # Performance indicators
install.packages("e1071")           # SVM
install.packages("mltools")         # MCC
install.packages("rlist")           # load/save list with .RData

# Required libraries
library(R.matlab)        # Allow us to read .mat files
library(ggplot2)         # Graphics
library(randomForest)    # Compute a random forest model
library(pROC)            # Performance indicators
library(e1071)           # SVM
library(mltools)         # MCC
library(rlist)           # load/save list with .RData


# Loading the functions
source("Code/derivees.R")                    # Functions for the rate of increase
source("Code/FonctionsIndTemp.R")            # Functions needed to compute HRV indicators
source("Code/FonctionGammaMetric.R")         # Gamma-metric

# Loading the datas
FA <- readMat("Data/RR_1min_FA_25_Sujets_Filtered_Shuffled.mat")$A
RS <- readMat("Data/RR_1min_RS_72_Sujets_Filtered_Shuffled.mat")$A

##### Compute the indicators #####
# Mean and standard deviation of the derivate up to order 10
FA.df <- create_ts_df(FA, ordre = 10)
RS.df <- create_ts_df(RS, ordre = 10)
donnees.df <- rbind(FA.df, RS.df)

# HRV indicators
FA.ind <- IndicateursTempWind(FA, 5)[[1]]$Resultats
RS.ind <- IndicateursTempWind(RS, 5)[[1]]$Resultats
donnees.ind <- rbind(FA.ind, RS.ind)

##### Construction of the data base #####
# Merge
donnees <- cbind(donnees.df, donnees.ind)

# Removing some columns
donnees <- donnees[,!colnames(donnees)%in%c("pNN50","HRV.index","SDNN")]

# Creation of the outcome column (defined as factor (1 : AF, 0 : RS))
donnees$Classe <- as.factor(c(rep("1",length(FA)), rep("0", length(RS))))

# Checks
head(donnees)
str(donnees)
summary(donnees)

# Parameters of the data.frame
n <- dim(donnees)[1]        # 50 028 observations
p <- dim(donnees)[2]-1      # 29 columns (without Classe, the outcome)
table(donnees$Classe)       # 2 900 AF and 49 128 RS
prev <- sum(donnees$Classe=="1")/n  # 5.796% AF


##### Construction of a random forest ####
mod <- randomForest(Classe~., data = donnees, importance = TRUE)

# Importance of each feature
Importance <- data.frame(Variables = row.names(mod$importance),
                         mod$importance,
                         row.names = NULL)

# sorted features in decreasing order of the importance score
MDA.sorted <- Importance[order(Importance$MeanDecreaseAccuracy, decreasing = TRUE),c("Variables","MeanDecreaseAccuracy")]
MDG.sorted <- Importance[order(Importance$MeanDecreaseGini, decreasing = TRUE), c("Variables", "MeanDecreaseGini")]

##### Ranking with the gamma-metric (by bootstrap on 500 replicates) #####
B <- 500
gamma.b <- matrix(rep(NA,p*B), ncol = p)
colnames(gamma.b) <- colnames(donnees)[1:p]
pb <- tkProgressBar(title = "", label = "Gamma-metric", min = 0, max = B, initial = 0)

for(b in 1:B)
{
  # Progress bar
  val <- as.character(round((b/B)*100,2))
  info <- sprintf("%s%% done", val)
  setTkProgressBar(pb = pb, value = b, label = info)
  
  # Computation of the gamma-metric for each features
  id <- sample(1:n, n, replace = TRUE)
  data.b <- donnees[id,]
  gamma.b[b,] <-  apply(data.b[,1:p], 2, gammaMetric, class = data.b$Classe, plot = FALSE)

}
close(pb)

# Computatin of the mean gamma-metric value over the 500 replicates
gamma <- data.frame(Variables = colnames(donnees)[1:p],
                    Gamma = apply(gamma.b, 2, median),
                    row.names = NULL)

# sorted features in decreasing order of the gamma-metric value
gamma.sorted <- gamma[order(gamma$Gamma, decreasing = TRUE),]


###### data.frame the features order with regard to the method (MDA, MDG, Gamma-metric) #####
sorted.var <- data.frame(ordre.MDA = MDA.sorted$Variables,
                         ordre.MDG = MDG.sorted$Variables,
                         ordre.gamma = gamma.sorted$Variables,
                         MDA = MDA.sorted$MeanDecreaseAccuracy, 
                         MDG = MDG.sorted$MeanDecreaseGini,
                         Gamma = gamma.sorted$Gamma)



# Classification performance with a SVM and features selected by MDA, MDG or Gamma-metric
RES.GLOB.GAMMA = list()
RES.GLOB.MDG = list()
RES.GLOB.MDA = list()

##### 5-fold cross validation repeated 10 times #####
K = 5 # number of folds
N = 10 # number of repetition
tot = K*N*p+1
iteration = 0

pb = tkProgressBar(title = "K-fold Cross Validation", label = "",min = 0, max = tot, initial = 0)

repartition = data.frame(Groupe = paste("Groupe",1:K, sep = ".")) # Data frame which will stock the proportion of AF in each fold and for each repetition

for(j in 1:N)
{
  cat("-----------------------------------------------------------\n")
  cat("                     ITERATIOn",j,"\n")
  cat("-----------------------------------------------------------\n")
  n = dim(donnees)[1]
  setTkProgressBar(pb, iteration, title = paste("Repetition",j), label = info)
  # Initialization of the data.frame in which we will stock the performance indicators of each iteration
  res.gamma = data.frame(Nbre.Variables = 1:p,
                   Variables.ajoutees = rep(0, p),
                   Specificite = rep(0, p),
                   Sensibilite = rep(0, p),
                   Accuracy = rep(0,p),
                   Mal.Classees = rep(0,p),
                   MCC = rep(0,p),
                   Precision = rep(0,p),
                   Rappel = rep(0,p),
                   AUC = rep(0,p))
  
  res.MDA = data.frame(Nbre.Variables = 1:p,
                         Variables.ajoutees = rep(0, p),
                         Specificite = rep(0, p),
                         Sensibilite = rep(0, p),
                         Accuracy = rep(0,p),
                         Mal.Classees = rep(0,p),
                         MCC = rep(0,p),
                         Precision = rep(0,p),
                         Rappel = rep(0,p),
                         AUC = rep(0,p))
  
  res.MDG = data.frame(Nbre.Variables = 1:p,
                         Variables.ajoutees = rep(0, p),
                         Specificite = rep(0, p),
                         Sensibilite = rep(0, p),
                         Accuracy = rep(0,p),
                         Mal.Classees = rep(0,p),
                         MCC = rep(0,p),
                         Precision = rep(0,p),
                         Rappel = rep(0,p),
                         AUC = rep(0,p))
  
  # For each observation, we randomly give a value between 1 and the number of folds (K = 5)
  donnees$Groupe = sample(1:K, n, replace = TRUE)
  
  # Compute the proportion of AF in each fold
  repartition[,j+1] <- sapply(1:K, FUN = function(k)
  {
    nb.FA <- sum(donnees[which(donnees$Groupe==k),]$Classe=="1")
    nb.obs <- nrow(donnees[which(donnees$Groupe==k),])
    return(nb.FA/nb.obs)
  })
  
  # Give the names of the features in decreasing order with regard to the method of ranking
  noms.gamma = paste(sorted.var$ordre.gamma)
  noms.mda = paste(sorted.var$ordre.MDA)
  noms.mdg = paste(sorted.var$ordre.MDG)
  
  # Each one of the K fold will be used as a testing set. 
  for(k in 1:K)
  {
    cat("Estimation des pr?dictions sur le groupe",k,"\n")
    train = donnees[donnees$Groupe != k, 1:(p+1)]   
    test = donnees[donnees$Groupe == k, 1:(p+1)]  # The observations from group k will be used as the testing set for this iteration
    
    # Vectors of the features used in the SVM model
    vars.gamma = NULL
    vars.mdg = NULL
    vars.mda = NULL
    
    # At each iteration of this loop we add a feature in the SVM model
    for(i in 1:p)
    {
      cat("Estimation pour",i,"variable(s) \n")

      # We add the next feature in the sorted vectors of features
      cat("La variable",noms.gamma[i],"est ajout?e (gamma) \n")
      cat("La variable", noms.mdg[i], "est ajout?e (mdg)\n")
      cat("La variable", noms.mda[i],"est ajout?e (mda)\n")
      res.gamma$Variables.ajoutees[i] = noms.gamma[i]
      res.MDG$Variables.ajoutees[i] = noms.mdg[i]
      res.MDA$Variables.ajoutees[i] = noms.mda[i]
      
      # Creation of the formula for the features ranked with the gamma-metric
      vars.gamma = c(vars.gamma, noms.gamma[i])
      formule.gamma = as.formula(paste("Classe~",paste(vars.gamma, collapse="+")))
      
      # Creation of the formula for the features ranked with MDG
      vars.mdg = c(vars.mdg, noms.mdg[i])
      formule.mdg = as.formula(paste("Classe~",paste(vars.mdg, collapse="+")))
      
      # Creation of the formula for the features ranked with MDA
      vars.mda = c(vars.mda, noms.mda[i])
      formule.mda = as.formula(paste("Classe~",paste(vars.mda, collapse="+")))
      
      # Creation of the SVM model for each ranking method
      cat("Mod?lisation SVM\n")
      mod.gamma = svm(formule.gamma, data = train, probability = TRUE)
      mod.mdg = svm(formule.mdg, data = train, probability = TRUE)
      mod.mda = svm(formule.mda, data = train, probability = TRUE)
      
      cat("Calcul des pr?dictions\n")
      # Computation of the predictions on the testing set
      pred.gamma = predict(mod.gamma, newdata = test, probability = TRUE)
      pred.mdg = predict(mod.mdg, newdata = test, probability = TRUE)
      pred.mda = predict(mod.mda, newdata = test, probability = TRUE)
      
      # Computation of the probabilities P(Y = 1| X = x)
      proba.gamma = attr(pred.gamma, "probabilities")[,"1"]
      proba.mdg = attr(pred.mdg, "probabilities")[,"1"]
      proba.mda = attr(pred.mda, "probabilities")[,"1"]
      
      
      # roc object will be used in order to derive the performance indicators (pROC package)
      roc1.gamma = roc(test$Classe, proba.gamma)
      roc1.mdg = roc(test$Classe, proba.mdg)
      roc1.mda = roc(test$Classe, proba.mda)
      
  
      cat("Calcul des indicateurs\n")
      # Computation of the performance indicators : Specificities, Sensibility, Accuracy, Error Rate, MCC, Precision, Recall and AUC
      # For the Gamma-metric ranking method
      indicateurs.gamma = as.numeric(coords(roc1.gamma, x="best",best.method="closest.topleft", ret = c("sp","se","accuracy", "tp", "tn", "fp", "fn")))
      res.gamma$Specificite[i] = res.gamma$Specificite[i] + indicateurs.gamma[1]
      res.gamma$Sensibilite[i] = res.gamma$Sensibilite[i] + indicateurs.gamma[2]
      res.gamma$Accuracy[i] = res.gamma$Accuracy[i] + indicateurs.gamma[3]
      res.gamma$Mal.Classees[i] = res.gamma$Mal.Classees[i] + (1-indicateurs.gamma[3])
      res.gamma$MCC[i] = res.gamma$MCC[i] + mcc(TP = indicateurs.gamma[4], FP = indicateurs.gamma[6], TN = indicateurs.gamma[5], FN = indicateurs.gamma[7])
      res.gamma$Precision[i] = res.gamma$Precision[i] + (indicateurs.gamma[4]/(indicateurs.gamma[4]+indicateurs.gamma[6]))
      res.gamma$Rappel[i] = res.gamma$Rappel[i] + (indicateurs.gamma[4]/(indicateurs.gamma[4]+indicateurs.gamma[7]))
      res.gamma$AUC[i] = res.gamma$AUC[i] + auc(roc1.gamma)
      
      # For the MDG ranking method
      indicateurs.mdg = as.numeric(coords(roc1.mdg, x="best",best.method="closest.topleft", ret = c("sp","se","accuracy", "tp", "tn", "fp", "fn")))
      res.MDG$Specificite[i] = res.MDG$Specificite[i] + indicateurs.mdg[1]
      res.MDG$Sensibilite[i] = res.MDG$Sensibilite[i] + indicateurs.mdg[2]
      res.MDG$Accuracy[i] = res.MDG$Accuracy[i] + indicateurs.mdg[3]
      res.MDG$Mal.Classees[i] = res.MDG$Mal.Classees[i] + (1-indicateurs.mdg[3])
      res.MDG$MCC[i] = res.MDG$MCC[i] + mcc(TP = indicateurs.mdg[4], FP = indicateurs.mdg[6], TN = indicateurs.mdg[5], FN = indicateurs.mdg[7])
      res.MDG$Precision[i] = res.MDG$Precision[i] + (indicateurs.mdg[4]/(indicateurs.mdg[4]+indicateurs.mdg[6]))
      res.MDG$Rappel[i] = res.MDG$Rappel[i] + (indicateurs.mdg[4]/(indicateurs.mdg[4]+indicateurs.mdg[7]))
      res.MDG$AUC[i] = res.MDG$AUC[i] + auc(roc1.mdg)
      
      # For the MDA ranking method
      indicateurs.mda = as.numeric(coords(roc1.mda, x="best",best.method="closest.topleft", ret = c("sp","se","accuracy", "tp", "tn", "fp", "fn")))
      res.MDA$Specificite[i] = res.MDA$Specificite[i] + indicateurs.mda[1]
      res.MDA$Sensibilite[i] = res.MDA$Sensibilite[i] + indicateurs.mda[2]
      res.MDA$Accuracy[i] = res.MDA$Accuracy[i] + indicateurs.mda[3]
      res.MDA$Mal.Classees[i] = res.MDA$Mal.Classees[i] + (1-indicateurs.mda[3])
      res.MDA$MCC[i] = res.MDA$MCC[i] + mcc(TP = indicateurs.mda[4], FP = indicateurs.mda[6], TN = indicateurs.mda[5], FN = indicateurs.mda[7])
      res.MDA$Precision[i] = res.MDA$Precision[i] + (indicateurs.mda[4]/(indicateurs.mda[4]+indicateurs.mda[6]))
      res.MDA$Rappel[i] = res.MDA$Rappel[i] + (indicateurs.mda[4]/(indicateurs.mda[4]+indicateurs.mda[7]))
      res.MDA$AUC[i] = res.MDA$AUC[i] + auc(roc1.mda)
      cat("\n")
      
      # Progress bar updates
      iteration = iteration + 1
      val <- as.character(round((iteration/tot)*100,2))
      info <- sprintf("%s%% done", val)
      setTkProgressBar(pb, iteration, title = paste("Repetition",j), label = info)
      
    }
    cat("\n")
  }
  
  # Each fold has been used as a testing set, we compute the mean of each indicators
  # Gamma
  Resultat.SVM.gamma = res.gamma
  Resultat.SVM.gamma[,3:10] = apply(res.gamma[,3:10], 2, FUN = function(x){x/K})
  RES.GLOB.GAMMA = append(RES.GLOB.GAMMA,list(Resultat.SVM.gamma))
  
  # MDG
  Resultat.SVM.MDG = res.MDG
  Resultat.SVM.MDG[,3:10] = apply(res.MDG[,3:10], 2, FUN = function(x){x/K})
  RES.GLOB.MDG = append(RES.GLOB.MDG,list(Resultat.SVM.MDG))
  
  # MDA
  Resultat.SVM.MDA = res.MDA
  Resultat.SVM.MDA[,3:10] = apply(res.MDA[,3:10], 2, FUN = function(x){x/K})
  RES.GLOB.MDA = append(RES.GLOB.MDA,list(Resultat.SVM.MDA))
}
close(pb)

###### Results ######

# Function used in order to create a single data.frame of the results over all the 5-fold cross validation 
aggregation <- function(data, operation, indices)
{
  res <- data.frame(do.call("cbind",lapply(indices, FUN = function(ind)
  {
    # For each repeat we extract the indicator we are interested in and call the operator 'operation' over each repeatition
    res.cv <- data.frame(do.call("cbind", lapply(data, "[[", ind)))
    res<- apply(res.cv, 1, FUN = function(x){do.call(operation, list(x))})
    return(res)
  })))
  colnames(res) <- indices
  return(res)
}

list.save(RES.GLOB.GAMMA, file = paste("RES_GLOB_GAMMA",Sys.Date(),".RData",sep=""))
list.save(RES.GLOB.MDG, file = paste("RES_GLOB_MDG",Sys.Date(),".RData", sep = ""))
list.save(RES.GLOB.MDA, file = paste("RES_GLOB_MDA",Sys.Date(),".RData", sep = ""))
write.table(x = repartition, file = paste("repartition",Sys.Date(),".txt",sep = ""))

# Indicators of performance we want to aggregate
indicateurs <- colnames(RES.GLOB.GAMMA[[1]])[!colnames(RES.GLOB.GAMMA[[1]])%in%c("Nbre.Variables","Variables.ajoutees")]

# Give us the mean/standard deviation of all indicators over the 10 k-fold cross validation for the Gamma-metric ranking method
GAMMA.MEAN <- data.frame(Variable.ajoutee = RES.GLOB.GAMMA[[1]]$Variables.ajoutees,
                         aggregation(RES.GLOB.GAMMA, "mean", indicateurs))

GAMMA.SD <- data.frame(Variale.ajoutee = RES.GLOB.GAMMA[[1]]$Variables.ajoutees, 
                       aggregation(RES.GLOB.GAMMA, "sd", indicateurs))

# Give us the mean/standard deviation of all indicators over the 10 k-fold cross validation for the MDG ranking method
MDG.MEAN <- data.frame(Variable.ajoutee = RES.GLOB.MDG[[1]]$Variables.ajoutees,
                       aggregation(RES.GLOB.MDG,"mean", indicateurs))
MDG.SD <- data.frame(Variable.ajoutee = RES.GLOB.MDG[[1]]$Variables.ajoutees,
                     aggregation(RES.GLOB.MDG, "sd", indicateurs))

# Give us the mean/standard deviation of all indicators over the 10 k-fold cross validation for the MDA ranking method
MDA.MEAN <- data.frame(Variable.ajoutee = RES.GLOB.MDA[[1]]$Variables.ajoutees,
                       aggregation(RES.GLOB.MDA,"mean", indicateurs))
MDA.SD <- data.frame(Variable.ajoutee = RES.GLOB.MDA[[1]]$Variables.ajoutees,
                     aggregation(RES.GLOB.MDA, "sd", indicateurs))


# Computation of the False Negative Rate and False Positive Rate
GAMMA.MEAN$FNR <- 1 - GAMMA.MEAN$Sensibilite
GAMMA.MEAN$FPR <- 1 - GAMMA.MEAN$Specificite

MDA.MEAN$FNR <- 1 - MDA.MEAN$Sensibilite
MDA.MEAN$FPR <- 1 - MDA.MEAN$Specificite

MDG.MEAN$FNR <- 1 - MDG.MEAN$Sensibilite
MDG.MEAN$FPR <- 1 - MDG.MEAN$Specificite


# Indicators
vars = c("Variable.ajoutee","Mal.Classees","FNR","FPR","AUC","Precision","Rappel","MCC")
vars.pourcent <- c("Mal.Classees","FNR","FPR")
vars.format <- vars[-1]

# TPR and TNR in pourcent
mat <- cbind(Nbr.Var = 1:p,GAMMA.MEAN[,vars],MDG.MEAN[,vars],MDA.MEAN[,vars])
mat[,colnames(mat)%in%vars.pourcent] <- 100*mat[,colnames(mat)%in%vars.pourcent]
mat[,colnames(mat)%in%vars.format] <- round(mat[,colnames(mat)%in%vars.format],3)

stargazer(as.matrix(mat)[,1:9])
stargazer(as.matrix(mat)[,c(1,10:17)])
stargazer(as.matrix(mat)[,c(1,18:25)])


###### Graphics #####
MCC <- data.frame(MCC.mean = c(GAMMA.MEAN$MCC, MDG.MEAN$MCC, MDA.MEAN$MCC),
                  MCC.up = c(GAMMA.MEAN$MCC+GAMMA.SD$MCC, MDG.MEAN$MCC + MDG.SD$MCC, MDA.MEAN$MCC + MDA.SD$MCC),
                  MCC.low = c(GAMMA.MEAN$MCC - GAMMA.SD$MCC, MDG.MEAN$MCC - MDG.SD$MCC, MDA.MEAN$MCC - MDA.SD$MCC),
                  Methode = c(rep("Gamma", p), rep("MDG",p), rep("MDA", p)))

FNR <- 100*(1 - GAMMA.MEAN$Sensibilite)
FPR <- 100*(1 - GAMMA.MEAN$Specificite)
Err <- 100*GAMMA.MEAN$Mal.Classees


max.mcc <- c(which.max(MCC$MCC.mean[MCC$Methode == "Gamma"]), max(MCC$MCC.mean[MCC$Methode == "Gamma"]))
min.fnr <- c(which.min(FNR), min(FNR))
min.fpr <- c(which.min(FPR), min(FPR))
min.err <- c(which.min(Err), min(Err))

par(mar = c(5,5,3,5), cex.axis=1.5)
plot(1:p, FNR, type = "l", col = "blue", xlab = "", ylab = "", lwd = 2)
points(x = min.fnr[1], y = min.fnr[2], pch = 19, col = "blue")
text(x = min.fnr[1]+0.5, y = min.fnr[2]+0.005, labels = round(min.fnr[2],4), col = "blue")
title(xlab = "Number of features", ylab = "Rate(%)", cex.lab = 1.5)
lines(1:p, FPR, type = "l", col = "red", lwd = 2)
points(x = min.fpr[1], y = min.fpr[2], pch = 19, col = "red")
text(x = min.fpr[1], y = min.fpr[2]+0.01, labels = round(min.fpr[2],4), col = "red")
lines(1:p, Err, type = "l", col = "black", lwd = 2)
points(x = min.err[1], y = min.err[2], col = "black", pch = 19)
text(x = min.err[1], y = min.err[2]+0.01, labels = round(min.err[2],4),col = "black")
par(new = TRUE)
plot(x = 1:p, y = MCC$MCC.mean[which(MCC$Methode == "Gamma")],
     type = "l", xaxt = "n", yaxt="n", ylab ="", xlab="", lty =2, col = "black", lwd = 2)
points(x = max.mcc[1], y = max.mcc[2], pch = 19, col = "black")
text(x = max.mcc[1], y = max.mcc[2]+0.0001, labels = round(max.mcc[2],4), col = "black")

axis(side = 4)
mtext("Matthew's Correlation Coefficient", side = 4, line = 3, cex = 1.5) 
legend(x = "topright",
       legend = c("False Negative Rate", "False Positive Rate", "Error Rate", "MCC"),
       col = c("blue","red","black","black"),
       lty = c(1,1,1,2), cex = 1, lwd = 2)

ggplot(data = MCC, aes(x = rep(1:p,3), y = MCC.mean, ymin = MCC.low, ymax = MCC.up, fill = Methode))+
  geom_line(aes(col = Methode), lwd = 2)+
  geom_ribbon(alpha = 0.5)




