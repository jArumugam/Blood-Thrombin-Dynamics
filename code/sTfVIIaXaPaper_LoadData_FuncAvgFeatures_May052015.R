foldercaseName <- 'Data/'
caseName <- 'May05'

# LOAD TIME AVERAGED FUNCTION VALUES
fileName <- paste(foldercaseName, 'ACS', caseName, 
                  '_5pM_FunctionAverages.txt', sep="")
print(fileName)
dataACSFuncAvg <- read.table(fileName)

fileName <- paste(foldercaseName, 'CAD', caseName, 
                  '_5pM_FunctionAverages.txt', sep="")
print(fileName)
dataCADFuncAvg <- read.table(fileName)

fileName <- paste(foldercaseName, 'ACSMay05_5pM_TimeStamp.txt', sep="")
print(fileName)
TimeStamp <- read.table(fileName)

#fix(TimeStamp)

#NIntervalsForFuncAVg <- ncol(TimeStamp)
SpNames <- fGetSpeciesName()
FeatureList <- fGetFuncAvgFeatureListMay052015(TimeStamp, SpNames)

names(dataACSFuncAvg) <- FeatureList
names(dataCADFuncAvg) <- FeatureList

#fix(dataACSFuncAvg)
#fix(dataCADFuncAvg)

noOfFuncValsForAllVar <- ncol(dataACSFuncAvg)

FuncValCols <- 1:noOfFuncValsForAllVar

ThrombinFeatureList = c('Thrombin-T2nM', 'Thrombin-MaxL', 'Thrombin-TMaxL', 
                        'Thrombin-MaxR', 'Thrombin-TMaxR', 'Thrombin-AUC')
XaFeatureList = c('fXa-MaxL', 'fXa-TMaxL', 
                  'fXa-MaxR', 'fXa-TMaxR', 'fXa-AUC')

noOfThrombinParams <- 6
noOfXaParams <- 5

ThrombinCols <- (noOfFuncValsForAllVar+1):(noOfFuncValsForAllVar+noOfThrombinParams)
XaCols <- (noOfFuncValsForAllVar+noOfThrombinParams+1):
  (noOfFuncValsForAllVar+noOfThrombinParams+noOfXaParams)

ThrombinCols
XaCols

fileName <- paste(foldercaseName, 'ACS', 'April26', 
                  '_5pM_ActiveThrombinSummaryParameters.txt', sep="")
print(fileName)
dataACS_IIa <- read.table(fileName)
names(dataACS_IIa) <- ThrombinFeatureList
#fix(dataACS_IIa)

fileName <- paste(foldercaseName, 'ACS', 'April26', 
                  '_5pM_XaSummaryParameters.txt', sep="")
print(fileName)
dataACS_Xa <- read.table(fileName)
names(dataACS_Xa) <- XaFeatureList
#fix(dataACS_Xa)


fileName <- paste(foldercaseName, 'CAD', 'April26', 
                  '_5pM_ActiveThrombinSummaryParameters.txt', sep="")
print(fileName)
dataCAD_IIa <- read.table(fileName)
names(dataCAD_IIa) <- ThrombinFeatureList
#fix(dataCAD_IIa)


fileName <- paste(foldercaseName, 'CAD', 'April26', 
                  '_5pM_XaSummaryParameters.txt', sep="")
print(fileName)
dataCAD_Xa <- read.table(fileName)
names(dataCAD_Xa) <- XaFeatureList
#fix(dataCAD_Xa)


dataACS <- data.frame(dataACSFuncAvg, dataACS_IIa)
dataACS <- data.frame(dataACS, dataACS_Xa)
dataCAD <- data.frame(dataCADFuncAvg, dataCAD_IIa)
dataCAD <- data.frame(dataCAD, dataCAD_Xa)

colnames(dataACS)[600:623]

disease = rep("ACS", 200)
dataACS <- data.frame(dataACS, disease)
colnames(dataACS)[ncol(dataACS)] <- 'Disease'

colnames(dataACS)[620:624]

disease = rep("CAD", 200)
dataCAD <- data.frame(dataCAD, disease)
colnames(dataCAD)[ncol(dataCAD)] <- 'Disease'

colnames(dataACS)[620:624]

data.FuncAvg <- rbind(dataACS, dataCAD)
colnames(data.FuncAvg)[610:624]


rm(dataACS, dataCAD,
 dataACSFuncAvg, dataCADFuncAvg,
 dataACS_IIa, dataCAD_IIa,
 dataACS_Xa, dataCAD_Xa)


ThrombinCols[1]
data.FuncAvg[,ThrombinCols[1]]
names(data.FuncAvg)[ThrombinCols]




