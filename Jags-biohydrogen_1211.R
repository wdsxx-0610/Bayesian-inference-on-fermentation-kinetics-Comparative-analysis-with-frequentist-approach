#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("/Users/wdsxx0610/Documents/R_directory/A_formation/biohydrogen/Jags-biohydrogen-source_1211.R")
#------------------------------------------------------------------------------- 
# #.............................................................................
# Only one predictor:
rawData = read.csv( file="/Users/wdsxx0610/Documents/R_directory/A_formation/biohydrogen/biohydrogenData.csv" )
xName = "time"
metabolites = c("biomass", "glucose", "hydrogen", "acetate", "lactate")
unit = c(" (g/L)", " (g/L)", " (mol/mol)", " (mM)", " (mM)")
numSavedSteps = 20000
thinSteps = 2
nChains = 4
graphFileType = "png" 
#
pName <- c("Amax", "a", "b", "c", "sigma")

init_biomass <- c(1.358670597,0.317526471,0.005644986,8.459471989,0.083241498)
init_glucose <- c(6.71516631,0.011362886,0.094433556,15.80765327,0.166956094)
init_hydrogen <-c(3.595692089,0.150386599,0.007338925,24.49573055,0.125281006)
init_acetate <- c(53.13713456,0.113432032,0.032782727,30.60206256,2.148719873)
init_lactate <- c(25.01659055,0.117430552,0.032105069,29.21516088,1.5140836)
init_guess <<- data.frame(pName, init_biomass, init_glucose, init_hydrogen, 
                          init_acetate, init_lactate)
#
  i = 4  #1 = biomass, 2 = glucose, 3 = hydrogen, 4 = acetate, 5 = lactate
  xData = rawData[,1]
  
  yName = metabolites[i]
  yData = rawData[,i+1]
  xUnit = " (h)"
  yUnit = unit[i]
#
  init_Amax <<- init_guess[1,i+1]
  init_a <<- init_guess[2,i+1]
  init_b <<- init_guess[3,i+1]
  
  init_c <<- init_guess[4,i+1]
  init_s <<- init_guess[5,i+1] #sigma
  fileNameRoot = paste0("/Users/wdsxx0610/Documents/R_directory/A_formation/biohydrogen/", yName, "/", yName, "_", "123")
  fileNameRoot
  outCSV = paste0(fileNameRoot, "bayesian.csv")
  outSMRY = paste0(fileNameRoot, "bayesian_summary.csv")
#
# Generate the MCMC chain:
mcmcCoda = genMCMC(rawData[,1], rawData[,i+1], 
                   numSavedSteps=numSavedSteps, thinSteps=thinSteps, 
                   nChains=nChains,saveName=fileNameRoot )
mcmcCoda
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:

parameterNames = varnames(mcmcCoda) # get all parameter names
#do not know why not work
#for ( parName in parameterNames ) {
#  diagMCMC( codaObject=mcmcCoda, parName=parName, 
#            saveName=fileNameRoot, saveType=graphFileType )
#}
parameterNames = varnames(mcmcCoda)   # 建议直接这样
print(parameterNames)

for (parName in parameterNames) {
  cat("Now plotting:", parName, "\n")
  diagMCMC(codaObject = mcmcCoda,
           parName     = parName,
           saveName    = fileNameRoot,
           saveType    = graphFileType)
}

#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC(mcmcCoda, saveName=fileNameRoot)
show(summaryInfo)
#summary(mcmcCoda)
#sink(paste0(fileNameRoot,"lactate.txt"), append=TRUE)
#cat("* * * * * * * * * * * * * * *\n")
#cat("* * * Summary(mcmcCoda) * * *\n")
#cat("* * * * * * * * * * * * * * *\n")
summary(mcmcCoda)

#
library(broom.mixed)
#cat("\n* * * * * * * * * * * * * * *")
#cat("\n* * * tidyMCMC(mcmcCoda) * * *")
#cat("\n* * * * * * * * * * * * * * *\n\n")
(tidyMCMC(mcmcCoda, conf.int=TRUE, conf.method="HPDinterval", ess=TRUE))
#sink()
#
#------------------------------------------------------------------------------- 
#prepare data for graphical summary
plotData = cbind(rawData[,1], rawData[,i+1])
yName <- metabolites[i]
colnames(plotData) <- c(xName, yName)
plotData <- data.frame(plotData)
yName = metabolites[i]
#
plotData
plotMCMC( mcmcCoda , data=plotData , xName=xName , yName=yName , 
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot, saveType=graphFileType )
#
#generate MCMC density plots
library(mcmcplots)
#
#主要的问题是在这个部分的某一个参数b或者a因为太稀疏所以没办法画出密度图
#png(paste0(fileNameRoot, "Post_Density.png"))
##  denplot(mcmcCoda, parms = c("Amax", "a","b", "c"))
#dev.off()
#
#generate MCMC trace plots
#
png(paste0(fileNameRoot, "Post_Trace.png"))
  traplot(mcmcCoda, parms = c("Amax", "a","b", "c"))
dev.off()
#
#generate pairs plot
library(ggplot2)
library("GGally")
#
png(paste0(fileNameRoot, "Post_Pairs.png"))
  ggpairs(data.matrix(mcmcCoda), columns=1:length(pName), lower = list(continuous = "smooth"))
#ggsave(paste0(fileNameRoot,"Post_Pairs.png"))
dev.off()
#
#generate residual plot
library(ggplot2)
#
mcmc=as.matrix(mcmcCoda[[nChains]])
coefs = apply(mcmc, 2, median) #use median value for Amax, a, b, c, and sigma
fitted = array(nrow(plotData))
for (j in 1:nrow(plotData)){
  fitted[j] = coefs[1]/(exp(-coefs[2]*(plotData$time[j]-coefs[4])) +
                        exp( coefs[3]*(plotData$time[j]-coefs[4])))
}
resid = plotData[,2] - fitted
ggplot() + 
  geom_point(data=NULL, aes(y=resid, x=fitted))  +
  geom_hline(yintercept=0) +
  labs(title='Residual vs. Fitted Values Plot', x='Fitted Values', y='Residuals')
ggsave(paste0(fileNameRoot, "Post_Residual.png"))
#
#calculate root mean square error
(rmse = sqrt(mean(resid^2)))
#calculate Nash-Sutcliffe model efficiency
meanplotData = mean(plotData[,2])
(nsme = 1 - sum((fitted-plotData[,2])^2)/sum((plotData[,2]-meanplotData)^2))
#
seme <-rbind("rmse", rmse, "nsme", nsme)
write.table(seme, outSMRY, append=TRUE, sep=",")
#generate metabolite profiles (both fitted and observation)
#
library(dplyr)
mcmc=as.matrix(mcmcCoda[[nChains]])
newTime = as.array(x = seq(min(plotData$time,na.rm=TRUE), max(plotData$time,na.rm=TRUE), len=100))
yRep = matrix(, nrow(mcmc), nrow(newTime)) #create yRep matrix object
for (i in 1:nrow(mcmc)){
  for (j in 1:nrow(newTime)){
    yRep[i,j] = mcmc[i,1]/(exp(-mcmc[i,2]*(newTime[j]-mcmc[i,4])) +
                             exp( mcmc[i,3]*(newTime[j]-mcmc[i,4])))  
  }
}
#yRep = matrix(, 3333, nrow(newTime)) #create yRep matrix object
#for (i in 16666:19998){
#  for (j in 1:nrow(newTime)){
#    yRep[i-16665,j] = mcmc[i,1]/(exp(-mcmc[i,2]*(newTime[j]-mcmc[i,4])) +
#                             exp( mcmc[i,3]*(newTime[j]-mcmc[i,4])))  
#  }
#}
#
newMetabolite = cbind(tidyMCMC(as.mcmc(yRep), conf.int=TRUE, conf.method="HPDinterval"))
xmax = max(xData)
ymax = max(newMetabolite$conf.high)
ggplot(newMetabolite, aes(y=estimate, x=newTime)) + 
  geom_line()+
  geom_point(data = plotData, aes(xData, yData))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="blue", alpha=0.3)+
  scale_x_continuous(name="time (h)", breaks=seq(0, ceiling(xmax/10)*10, 5)) + 
  scale_y_continuous(name=paste0(yName," (g/L)", breaks=seq(0, ceiling(ymax/1)*1, 0.2))) +
  guides(x=guide_axis(minor.ticks=TRUE), y=guide_axis(minor.ticks=TRUE)) +
  theme(panel.background=element_rect(fill="white",color="grey50")) +
  annotate("text", x=0.7*max(xData), y=0.8*max(yData), 
           label=paste0("rmse = ", sprintf("%1.4f", rmse))) +
  annotate("text", x=0.7*max(xData), y=0.7*max(yData), 
           label=paste0("nsme = ", sprintf("%1.4f", nsme))) +  
  ggtitle("Bayesian")   
ggsave(paste0(fileNameRoot, "Post_", yName, "Profile.png"))
#
#result output
out1 <- cbind(newTime, newMetabolite[,2:5]) #predicted metabolite profile
HPDinterval <- newMetabolite[,5] - newMetabolite[,4] #conf.high - conf.low
colnames(out1)[1] <- "time"
p1 <- rep(plotData[,1], length = nrow(out1))
p2 <- rep(plotData[,2], length = nrow(out1))
nr = nrow(rawData) + 1

p11 <- replace(p1,nr:nrow(out1),0)
p21 <- replace(p2,nr:nrow(out1),0)
p12 <- cbind (p11, p21)
colnames(p12) <- c("time", yName)
output <- cbind(p12, out1, HPDinterval)
write.csv(output, outCSV)
#
#effect size
mcmc=as.matrix(mcmcCoda[[nChains]])
newTime = as.array(x = c(min(rawData$time,na.rm=TRUE), max(rawData$time,na.rm=TRUE)))
fit = matrix(, nrow(mcmc), nrow(newTime)) #create yRep matrix object
for (i in 1:nrow(mcmc)){
  for (j in 1:nrow(newTime)){
    fit[i,j] = mcmc[i,1]/(exp(-mcmc[i,2]*(newTime[j]-mcmc[i,4])) +
                             exp( mcmc[i,3]*(newTime[j]-mcmc[i,4])))  
  }
}
#raw effect size
(RES = tidyMCMC(as.mcmc(fit[,2]-fit[,1]), conf.int=TRUE, conf.method="HPDinterval"))
#Cohen's D
cohenD = (fit[,2] - fit[,1])/mcmc[,"sigma"]
(cohenDES = tidyMCMC(as.mcmc(cohenD), conf.int=TRUE, conf.method="HPDinterval"))
#percentage change (relative to Group A)
ESp = 100 * (fit[,2] - fit[,1])/fit[,1]
(PES = tidyMCMC(as.mcmc(ESp), conf.int=TRUE, conf.method="HPDinterval"))

#find mode in a range of data
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#example: find mode in Amax
mode(mcmc[,1])

# 提取 Test C 的标准差 (Standard Deviation)
# 假设你已经跑完了 Test C 的代码，且 mcmcCoda 对象还在内存里
TestC_SD <- apply(as.matrix(mcmcCoda), 2, sd)
print(TestC_SD)

