# This is the script which processes for the Tf-fVIIa-fXa paper 

rm(list=ls())

# load libraries
library(tree)
library(randomForest)
library(xtable)
library(ggplot2)
# library(zoo) 
library('reshape2')

# load functions
source('sPaperResults_Functions_Oct042015.R')
# source('fAccesories.R')
# source('fPrintTablesTf.R')
# source('fPlotFigures.R')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data
source('sTfVIIaXaPaper_LoadData_FuncAvgFeatures_May052015.R')

# fix(data.FuncAvg)
# str(data.FuncAvg)
# write.table(data.FuncAvg, file = "Data_623features_624DiseaseName", 
#             append = FALSE, sep = " ",
#             eol = "\n", na = "NA", dec = ".", 
#             row.names = TRUE, col.names = TRUE)

# aTmp <- read.table('Data_623features_624DiseaseName')
# colnames(aTmp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit Trees for specific Variables and get mean OOB #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# varNo 
#    No.  Species 
#    10   Tf-fVIIa-fXa
#     9   Tf-fVIIa-fX
#     18  fIXa-fVIIIa-fX
#     7   IIa

varNo = 10

varColNos <- fGetColNosForVariableFuncAvgMay09(TimeStamp, varNo) 
varColNos

kForKFold <- 50
varfuncAvg.RF.OOB.MDGini <- getOOBnMDGiniFromMultipleRFM(x = data.FuncAvg[, varColNos],
                                                   y = data.FuncAvg$Disease,
                                                   KFoldNo = kForKFold, 
                                                   TreeCount = 501)

print.mean.sd.Accuracy.from.RF.Result(varfuncAvg.RF.OOB.MDGini)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#
# Data Tests #
#~~~~~~~~~~~~#
ncol(data.FuncAvg)
names(data.FuncAvg)[ncol(data.FuncAvg)-1]
names(data.FuncAvg)[ncol(data.FuncAvg)]
names(data.FuncAvg)[ncol(data.FuncAvg)-11]
# fix(FuncValCols)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit Trees for funcAvg and get mean OOB #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
last(FuncValCols)
kForKFold <- 50

# funcAvg.RF.OOB.MDGini <- getKFoldOOBnMDGiniFromRFM(x = data.FuncAvg[, FuncValCols],
funcAvg.RF.OOB.MDGini <- getOOBnMDGiniFromMultipleRFM(x = data.FuncAvg[, FuncValCols],
                                                   y = data.FuncAvg$Disease,
                                                   KFoldNo = kForKFold, 
                                                   TreeCount = 501)

print.mean.sd.Accuracy.from.RF.Result(funcAvg.RF.OOB.MDGini)

FuncValColNames <- colnames(data.FuncAvg)[FuncValCols]
length(FuncValColNames)

funcAVg.MDGini.KFold <- funcAvg.RF.OOB.MDGini[[2]]
funcAVg.MDGini.KFold.Mean <- colMeans(funcAVg.MDGini.KFold)
# fix(funcAVg.MDGini.KFold.Mean)
str(funcAVg.MDGini.KFold.Mean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~#
# VariableImportance #
#~~~~~~~~~~~~~~~~~~~~#
varImpFuncAvgKFoldDF <- 
  reshapeFuncAvgFeaturesToVariablesCols(funcAVg.MDGini.KFold.Mean,
                                        FuncValCols, TimeSTamp, 
                                        SpNames = SpNames)
# fix(varImpFuncAvgKFoldDF)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot SVariable Importance #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
TimeStampData <-  unlist(as.list(unname(TimeStamp)))
TimeStampData

TimeStamp.varImpFuncAvgDF <- data.frame(TimeStampData, varImpFuncAvgKFoldDF[,c(6,7,9,10,18)])
# TimeStamp.varImpFuncAvgDF <- data.frame(TimeStampData, varImpFuncAvgKFoldDF[,])
melted.DF = melt(TimeStamp.varImpFuncAvgDF, id.vars="TimeStampData")
# fix(melted.DF)

plt <- qplot(TimeStampData, value, data = melted.DF, 
             geom = c("line","point"), color = variable)
plt <- plt + labs(title = "MDGini for Selected Variables") 
plt <- plt + labs(x = "Time (sec)", y = "MDGini")
plt <- plt + guides(color=guide_legend(title="Species"))
plt <- plt + theme(axis.title = element_text(size=rel(1)), 
                   title = element_text(size=rel(1)), 
                   plot.title = element_text(hjust = 0.1, vjust = 2),
                   legend.text=element_text(size=rel(0.5)), 
                   legend.title=element_blank())

# tiff("plot3.tiff", width=2250, height=1200, res=600)
print(plt)
dev.off()

# plt + theme(axis.title = element_text(size=rel(1.5)))
# ggplot(data=melted.DF, aes(x=TimeStampData, y=value, group=variable)) + geom_line()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# TimeStamp index and time with higest MDGini
# June 07
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Tf-fVIIa-fXa    - 07 - 1300
# Tf-fVIIa-fX     - 07 - 1300
# fIXa-fVIIIa-fX  - 18 - 3500
# IIa             - 10 - 1900
# fXa             - 02 - 0300

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Tf-fVIIa-fXa 3 High Measurement
SelectedCols <- c(18*(10-1)+7, 18*(10-1)+6, 18*(10-1)+8)
# Tf-fVIIa-fX 3 High Measurement
SelectedCols <- c(18*(9-1)+7, 18*(9-1)+6, 18*(9-1)+8)
# fIXa-fVIIIa-fX 3 High Measurement
SelectedCols <- c(18*(18-1)+17, 18*(18-1)+18, 18*(18-1)+16)
# IIa Few Measurement
SelectedCols <- c(18*6+10, 18*6+11, 18*6+12)
# fXa 3 High Measurement
SelectedCols <- c(18*(6-1)+2, 18*(6-1)+1, 18*(6-1)+3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Combinations
SelectedCols <- c(18*9+8, 18*6+10, 18*17+18)
SelectedCols <- c(18*8+8, 18*6+10, 18*17+18)
SelectedCols <- c(18*9+8, 18*8+8, 18*17+18)
SelectedCols <- c(18*9+8, 18*6+10, 18*8+8)

TimeStamp[c(8, 10, 18)]

# SelectedCols <- c(18*9+8, 18*6+8, 18*17+18)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~#
# Selected Few #
#~~~~~~~~~~~~~~#
TimeStamp
SelectedCols

kForKFold <- 50
# selectFew.RF.OOB.MDGini <- getKFoldOOBnMDGiniFromRFM(x = data.FuncAvg[, SelectedCols],
selectFew.RF.OOB.MDGini <- getOOBnMDGiniFromMultipleRFM(x = data.FuncAvg[, SelectedCols],
                                                     y = data.FuncAvg$Disease,
                                                     KFoldNo = kForKFold,
                                                     TreeCount = 501)

# selectFew.OOB.Rates <- selectFew.RF.OOB.MDGini[[1]]
SelectedCols
colnames(data.FuncAvg)[SelectedCols]
print.mean.sd.Accuracy.from.RF.Result(selectFew.RF.OOB.MDGini)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
colnames(data.FuncAvg)[SelectedCols]
#~~~~~~~~~~~~~~~~#
# Single Feature #
#~~~~~~~~~~~~~~~~#
# Tf-fVIIa-fXa    - 07 - 1300 - 18*9 
# Tf-fVIIa-fX     - 07 - 1300 - 18*8
# fIXa-fVIIIa-fX  - 18 - 3500 - 18*17
# IIa             - 10 - 1900 - 18*6
# fXa             - 02 - 0300 - 18*5

singleFeatureCols <- c(18*9 +07, 18*8 +07, 18*17 +18, 18*6 +10, 18*5 +2)
colnames(data.FuncAvg)[singleFeatureCols]

rf.tree.formula <-Disease~Tf.fVIIa.fXaat1500sec
rf.tree.formula <-Disease~Tf.fVIIa.fXat1300sec
rf.tree.formula <-Disease~fIXa.fVIIIa.fXat3500sec
rf.tree.formula <-Disease~fIIaat1900sec
rf.tree.formula <-Disease~fXaat300sec


# rf.tree1.OOB.MDGini <- getKFold.OOB.MDGini.RF.1Formula(formula.RF = rf.tree.formula, 
rf.tree1.OOB.MDGini <- get.OOB.MDGini.MultipleRF.1Formula(formula.RF = rf.tree.formula, 
                                                       data.RF = data.FuncAvg, 
                                                       KFoldNo = 50, TreeCount = 501)

rf.tree.formula
print.mean.sd.Accuracy.from.RF.Result(rf.tree1.OOB.MDGini)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SelectedCols <- c(18*17+18, 18*6+10)
SelectedCols <- c(18*9+7, 18*17+18)
SelectedCols

TwoFeature.DF = data.FuncAvg[,SelectedCols]
TwoFeature.DF = cbind(TwoFeature.DF, data.FuncAvg$Disease)
colnames(TwoFeature.DF)[1:3] <- c('TfVIIaXa','IXVIIIaX', 'Disease')
# fix(TwoFeature.DF)
tree.2Feature = tree(Disease~., 
                     data = TwoFeature.DF)
summary(tree.2Feature) 

randomColOrder <- sample(400)
classNames <- data.FuncAvg[randomColOrder, 'Disease']

# tiff("plot2.tiff", width=2250, height=2250, res=600)
par(mar=c(3.3,3.5,1.1,1.1), mgp=c(2.1, 1, 0))
plot(data.FuncAvg[randomColOrder,SelectedCols[1]], 
     data.FuncAvg[randomColOrder,SelectedCols[2]], 
     pch=1, cex = 1, cex.axis = 0.8, cex.lab=0.95, cex.main=0.95,
     xlab='Tf-fVIIa-fXa (M)',
     ylab='fIXa-fVIIIa-fX (M)',
     col=c("blue", "red")[unclass(classNames)])

partition.tree(tree.2Feature, add=TRUE, lwd=3, cex = 1)
legend("topright", 
       legend=c('ACS','CAD'), 
       col=c('blue','red'), 
       pch=1, cex = 1, pt.cex = 1, text.font = 1, y.intersp = 1.5, 
       bg = 'gray')

# tiff("Plot1.tiff")
# dev.off()





