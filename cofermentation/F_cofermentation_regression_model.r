#
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
# imports library
library(minpack.lm)
library(ggplot2)
# generate data
dataroots = "../data/"
rawData = read.csv( file=paste0(dataroots,"CofermentationData.csv") )
head(rawData)
xName = "time"
metabolites = c( "Lactate","Acetate", "Ethanol", "PDO")
x = rawData[,1]
y = rawData[,5]
yName = metabolites[4]
fileNameRoot = paste0("./results/", yName, "/", yName, "_")
outCSV = paste0(fileNameRoot, "frequentist.csv")
outSMRY = paste0(fileNameRoot, "frequentist_summary.csv")
# fit the model 
which_level = "confidence"
#which_level = "prediction"

#start_values <- c(Amax=56.2093, a=0.175, b=0.0046, c=25.9546)#Acetate
#start_values <- c(Amax=7.145, a=0.6076, b=0.0414, c=14.9552)#Ethanol
#start_values <- c(Amax=15.5552, a=0.2198, b=0.1899, c=20.7419)#Lactate
start_values <- c(Amax=55.3459*2, a=0.2979, b=0.0047, c=32.0571)#PDO
model <- nls(y ~ Amax / (exp(-a*(x-c)) + exp(b*(x-c))),
		start = start_values,
		algorithm = "port",
		control = nls.control(maxiter = 1000))
s0 <- summary(model)
s1 <- s0$parameters

s2 <- confint(model, level=0.95)
sumy <- cbind(s1, s2)

write.csv (sumy, outSMRY)
#
#calc root mean square error and Nash Sutcliffe model efficiency
pred_y =  predict(model, x)

resid = y - pred_y
rmse = sqrt(mean(resid^2))
mean_y = mean(y)
nsme = 1 - sum((resid)^2)/sum((y-mean_y)^2)
seme <- rbind("rmse", rmse, "nsme", nsme)
write.table(seme, outSMRY, append=TRUE, sep=",",col.names=FALSE)
#
library(investr)
nr = nrow(rawData)

newX <- data.frame(x = seq(0, rawData[nr,1], len=100))
fit <- investr::predFit(model, newX, interval = which_level, level=0.95)
head(fit)
#fit <- investr::predFit(model, newX, interval = "prediction", level=0.95)
newY <- fit[,1]
lwr <- fit[,2]
upr <- fit[,3]
#
# plotting
obsData = data.frame(x, y)
obsData
newData = data.frame(newX, newY)
head(newData)
colnames(newData)[1] <- "newX"
xmax = max(x)
ymax = max(upr)
ggplot(data = newData, aes(newX, newY))+
  geom_line()+
  geom_point(data = obsData, aes(x, y))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill="blue", alpha=0.3)+
  scale_x_continuous(name="time (h)", breaks=seq(0, ceiling(xmax/10)*10, 10)) + 
  scale_y_continuous(name=paste0(yName," (g/L)"), breaks=seq(0, ceiling(ymax/5)*5, 5)) +
  guides(x=guide_axis(minor.ticks=TRUE), y=guide_axis(minor.ticks=TRUE)) +
  theme(panel.background=element_rect(fill="white",color="grey50")) +
  annotate("text", x=0.7*max(x), y=0.8*max(y), 
           label=paste0("rmse = ", sprintf("%1.4f", rmse))) +
  annotate("text", x=0.7*max(x), y=0.7*max(y), 
           label=paste0("nsme = ", sprintf("%1.4f", nsme))) +  
  ggtitle(paste0("Frequentist", " (95% ", which_level," level)")) 
ggsave(paste0(fileNameRoot, "frequentist (", which_level, ")_", yName, "_Profile.png"))
#
#result output
confInterval <- upr - lwr #conf.high - conf.low
out1 <- cbind(newX, newY, lwr, upr, confInterval) #predicted metabolite profile
colnames(out1)[1:2] <- c("time", yName)
p1 <- rep(x, length = nrow(out1))
p2 <- rep(y, length = nrow(out1))
nr = nrow(rawData) + 1
p11 <- replace(p1,nr:nrow(out1),0)
p21 <- replace(p2,nr:nrow(out1),0)
p12 <- cbind (p11, p21)#true point value
colnames(p12) <- c("time", yName)
output <- cbind(p12, out1)
write.csv(output, outCSV)
#
head(newData)
