#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fGetSpeciesName <- function(){
  fileName = 'Data/SpeciesNames.txt'
  readChar(fileName, file.info(fileName)$size)
  SpNames <- readLines((fileName))
  return(SpNames)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fGetFuncAvgFeatureListMay052015 <- function(TimeSTamp, varNames){
  FeatureList = list()
  for (name in varNames) {
    for (AtTime in TimeSTamp) { 
      FeatureName = paste(name,'at',toString(AtTime),'sec', sep = "")
      FeatureList <- c(FeatureList, FeatureName)            
    }
  }
  FeatureList <- unlist(FeatureList)
  return(FeatureList)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Get ColNos for a given variable #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fGetColNosForVariableFuncAvgMay09 <- function(TimeStamp, varNo){ 
  noOfAvgFuncValsPerVariable <- ncol(TimeStamp)
  ColStartofithVariable <- ((varNo-1)*noOfAvgFuncValsPerVariable) + 1
  ColEndofithVariable <- ColStartofithVariable + noOfAvgFuncValsPerVariable-1
  colNos <- ColStartofithVariable:ColEndofithVariable
  return(colNos)
}
# fGetColNosForVariableFuncAvgMay09(TimeStamp, 1)
#fGetColNosForVariableFuncAvgMay09(TimeStamp, 2)
#fGetColNosForVariableFuncAvgMay09(TimeStamp, 34)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
getOOBnMDGiniFromMultipleRFM <- function(x, y, KFoldNo, TreeCount) {
  nFeatures <- ncol(x)
  
  importanceFlag = TRUE
  #   printOutNo = 250
  printOutNo = 500
  
  OOB.Rates.3Vals.DF <- data.frame(matrix(ncol = 3, 
                                          nrow = KFoldNo))
  colnames(OOB.Rates.3Vals.DF) <- c('OOB', 'ACS', 'CAD' )  
  Mean.Dec.Gini.DF <- matrix(ncol = nFeatures, nrow = KFoldNo) 
  
  for (iter in 1:KFoldNo) {
    OOB.Rates.3Vals <- 0.0
    Mean.Dec.Gini <- 0.0
    x.train <- x
    y.train <- y
    rf.tree <- randomForest(x.train, y.train,
                            importance = importanceFlag,
                            do.trace = printOutNo, 
                            # TreeCount = TreeCount)
                            ntree = TreeCount)
    
    OOB.Rates.3Vals <- getACSCAD.RFM.OOBRatesMay102014(rf.tree)
    Mean.Dec.Gini <- rf.tree$importance[,'MeanDecreaseGini']
    
    Mean.Dec.Gini <- unname(Mean.Dec.Gini, force=F)
    OOB.Rates.3Vals.DF[iter,] <- OOB.Rates.3Vals
    Mean.Dec.Gini.DF[iter,] <- Mean.Dec.Gini
    
  }
  outList <- list('OOBRate'=OOB.Rates.3Vals.DF,'MDGini'=Mean.Dec.Gini.DF)
  return(outList) 
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
getACSCAD.RFM.OOBRatesMay102014 <- function(rf.tree){
  # Returns Last element in 'nCol' th column of OOB object
  print('getACSCAD.RFM.OOBRatesMay102014 function Called')
  print('returning LAST OOB rate, ACS rate, CAD rate')
  
  OOB <- rf.tree$err.rate[,1]
  OOB.ACS <- rf.tree$err.rate[,2]
  OOB.CAD <- rf.tree$err.rate[,3]
  
  outList <- list(meanOOB=last(OOB), 
                  meanOOBACS=last(OOB.ACS), 
                  meanOOBCAD=last(OOB.CAD))
  
  return(outList)
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
get.OOB.MDGini.MultipleRF.1Formula <- function(formula.RF, data.RF, KFoldNo, TreeCount) {
  
  OOB.Rates.3Vals.DF <- data.frame(matrix(ncol = 3, nrow = KFoldNo))
  colnames(OOB.Rates.3Vals.DF) <- c('OOB', 'ACS', 'CAD')
  
  Mean.Dec.Gini.DF <- matrix(ncol = 1, nrow = KFoldNo) 
  
  for (iter in 1:KFoldNo) {
    
    OOB.Rates.3Vals <- 0.0
    Mean.Dec.Gini <- 0.0
    
    rf.tree <- randomForest(formula.RF, data = data.RF, 
                            importance = TRUE, do.trace = FALSE, 
                            ntree = TreeCount)
    
    OOB.Rates.3Vals <- getACSCAD.RFM.OOBRatesMay102014(rf.tree)
    Mean.Dec.Gini <- rf.tree$importance[,'MeanDecreaseGini']
    
    Mean.Dec.Gini <- unname(Mean.Dec.Gini, force=F)
    OOB.Rates.3Vals.DF[iter,] <- OOB.Rates.3Vals
    Mean.Dec.Gini.DF[iter,] <- Mean.Dec.Gini
  }
  outList <- list('OOBRate'=OOB.Rates.3Vals.DF,'MDGini'=Mean.Dec.Gini.DF)
  return(outList) 
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
last <- function(x) { tail(x, n = 1) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
print.mean.sd.Accuracy.from.RF.Result <- function(RF.OOB.MDGini){
  
  OOB.mean <- apply(RF.OOB.MDGini[[1]], 2, mean)
  OOB.sd <- apply(RF.OOB.MDGini[[1]], 2, sd)
  
  print("Mean")
  print(100.0-100.0*OOB.mean)
  print("SD")
  print(100.0*OOB.sd)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
reshapeFuncAvgFeaturesToVariablesCols <- function(Importance,
                                                  FuncValCols, TimeSTamp, 
                                                  SpNames = SpNames){
  
  noOfFuncValsPerVar <- length(TimeStamp)
  variableImportanceDF <- data.frame(matrix(ncol = 34, 
                                            nrow = noOfFuncValsPerVar))
  for (varNo in 1:34){
    colsOfGivenVar <- fGetColNosForVariableFuncAvgMay09(TimeStamp, varNo)
    importanceOfGivenVarSorted <- Importance[colsOfGivenVar]
    
    variableImportanceDF[,varNo] <- importanceOfGivenVarSorted
    names(variableImportanceDF)[varNo] <- SpNames[varNo]   
  }
  return(variableImportanceDF)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



