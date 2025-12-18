
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-cofermentation-source1211.r")
#------------------------------------------------------------------------------- 
# #.............................................................................
# Only one predictor:
rawData = read.csv( file="../data/CofermentationData.csv")

head(rawData)
xName = "time"
#Lactate	Acetate	Ethanol	PDO
metabolites = c("Lactate","Acetate", "Ethanol", "PDO")
numSavedSteps=20000 ; thinSteps=2
graphFileType = "png" 
#
pName <- c("Amax", "a","b", "c", "sigma")

init_Acetate <- c(56.2093, 0.175, 0.0046, 25.9546, 1.3968)
init_Ethanol <- c(7.145, 0.6076, 0.0414, 14.9552, 0.3864)
init_Lactate <- c(15.5552, 0.2198, 0.1899, 20.7419, 0.5924)
init_PDO     <- c(55.3459, 0.2979, 0.0047, 32.0571, 1)

init_guess <<- data.frame(pName, init_Lactate,init_Acetate, init_Ethanol ,init_PDO)
#
init_guess
  i = 1  # 1=lactate 2= Acetate, 3 = Ethanol, 4 = PDO
  xData = rawData[,1]
  yName = metabolites[i]
  yData = rawData[,i+1]
  #init_value <<- init_guess[,i+1]
  init_Amax <<- init_guess[1,i+1]
  init_a <<- init_guess[2,i+1]
  init_b <<- init_guess[3,i+1]
  init_c <<- init_guess[4,i+1]
  init_s <<- init_guess[5,i+1]
  fileNameRoot = paste0("./results/", yName, "/", yName, "_","1211")
  fileNameRoot
  outCSV = paste0(fileNameRoot, "bayesian.csv")
  outSMRY = paste0(fileNameRoot, "bayesian_summary.csv")
outCSV
outDir <- "./results/"
  #
# Generate the MCMC chain:
mcmcCoda = genMCMC(rawData[,1], rawData[,i+1], 
                   numSavedSteps=numSavedSteps, thinSteps=thinSteps , 
                   saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
mcmcCoda
parameterNames
graphFileType
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda, parName=parName, 
            saveName=fileNameRoot, saveType=graphFileType )
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
plotMCMC( mcmcCoda , data=plotData , xName=xName , yName=yName , 
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot, saveType=graphFileType )
#
#generate MCMC density plots
library(mcmcplots)
#
png(paste0(fileNameRoot, "Post_Density.png"))
  denplot(mcmcCoda, parms = c("Amax", "a","b", "c"))
dev.off()
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
mcmc=as.matrix(mcmcCoda)
coefs = apply(mcmc, 2, median) #use median value for Amax, a, b, c, and sigma

fitted = array(nrow(plotData))
for (j in 1:nrow(plotData)){
  # Use unified MHSF 4-parameter formula
  t_val = plotData$time[j]
  fitted[j] = coefs["Amax"] / (exp(-coefs["a"]*(t_val - coefs["c"])) + 
                                 exp( coefs["b"]*(t_val - coefs["c"])))
}
head(plotData)
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
mcmc = as.matrix(mcmcCoda)
newTime = seq(min(plotData$time, na.rm=TRUE), max(plotData$time, na.rm=TRUE), length=100)
yRep = matrix(NA, nrow=nrow(mcmc), ncol=length(newTime))

for (k in 1:nrow(mcmc)){
  p <- mcmc[k,] # Extract kth row parameters
  # Unified formula
  yRep[k,] = p["Amax"] / (exp(-p["a"]*(newTime - p["c"])) + 
                            exp( p["b"]*(newTime - p["c"])))
}

#
newMetabolite = cbind(tidyMCMC(as.mcmc(yRep), conf.int=TRUE, conf.method="HPDinterval"))
xmax = max(xData)
ymax = max(newMetabolite$conf.high)

ggplot(newMetabolite, aes(y=estimate, x=newTime)) + 
  geom_line()+
  geom_point(data = plotData, aes(xData, yData))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="blue", alpha=0.3)+
  scale_x_continuous(name="time (h)", breaks=seq(0, ceiling(xmax/10)*10, 10)) + 
  scale_y_continuous(name=paste0(yName," (g/L)"), breaks=seq(0, ceiling(ymax/5)*5, 5)) +
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
colnames(p12) <- c("time", "lactate")
output <- cbind(p12, out1, HPDinterval)
write.csv(output, outCSV)
#
#effect size
mcmc=as.matrix(mcmcCoda)
newTime = as.array(x = c(min(rawData$time,na.rm=TRUE), max(rawData$time,na.rm=TRUE)))
# effect size
# ... (newTime 定义保持不变) ...
fit = matrix(NA, nrow(mcmc), length(newTime)) # 修正 ncol 定义

for (k in 1:nrow(mcmc)){
  p <- mcmc[k,]
  # 统一公式
  fit[k,] = p["Amax"] / (exp(-p["a"]*(newTime - p["c"])) + 
                           exp( p["b"]*(newTime - p["c"])))

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


