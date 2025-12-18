#
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
# imports library
library(minpack.lm)
library(ggplot2)
# generate data
rawData = read.csv( file="../data/biohydrogenData.csv" )
head(rawData)
xName = "time"
metabolites = c("biomass", "glucose", "hydrogen", "acetate", "lactate")
i=5
x = rawData[,1]
y = rawData[,i+1]
yName = metabolites[i]
yName
fileNameRoot = paste0("./results/", yName, "/", yName, "_","123")
outCSV = paste0(fileNameRoot, "frequentist.csv")
outSMRY = paste0(fileNameRoot, "frequentist_summary.csv")
# fit the model 
#start_values <- c(Amax=100, a=0.06, b=0.06, c=55)
init_biomass <- c(1.358670597,0.317526471,0.005644986,8.459471989,0.083241498)
init_glucose <- c(6.71516631,0.011362886,0.094433556,15.80765327,0.166956094)
init_hydrogen <-c(3.595692089,0.150386599,0.007338925,24.49573055,0.125281006)
init_acetate <- c(53.13713456,0.113432032,0.032782727,30.60206256,2.148719873)
init_lactate <- c(25.01659055,0.117430552,0.032105069,29.21516088,1.5140836)

#start_values <- c(Amax=1.358670597, a=0.317526471, b=0.005644986, c=8.459471989)#biomass
#start_values <- c(Amax=init_glucose[1]*2, b=init_glucose[3], c=init_glucose[4])
#start_values <- c(Amax=init_hydrogen[1], a=init_hydrogen[2], c=init_hydrogen[4])
#start_values <- c(Amax=init_acetate[1], a=init_acetate[2], b=init_acetate[3], c=init_acetate[4])
start_values <- c(Amax=init_lactate[1], a=init_lactate[2], b=init_lactate[3], c=init_lactate[4])
#!!!How to determine the start_values
#start_values <- c(Amax=100, a=0.06,  c=55)
#model <- nls(y ~ Amax / (1+exp(-a*(x-c))),
#model <- nls(y ~ Amax / (1+exp(b*(x-c))),    
# Does not converge
# Considering the monotonic decreasing characteristic of glucose consumption in biohydrogen fermentation,
# and combined with model identifiability analysis (parameters a,b satisfy the collinearity constraint a ≈ -b,
# leading to singular Fisher information matrix), we adopt the single-sided logistic
# model f(t) = Amax/(1+exp(b(t-c))). This model is numerically stable.
model <- nls(y ~ Amax / (exp(-a*(x-c)) + exp(b*(x-c))),
		start = start_values,
		algorithm = "port",
		control = nls.control(maxiter = 1000))

#library(minpack.lm)
#model <- nlsLM(y ~ Amax / (exp(-a*(x-c)) + exp(b*(x-c))),
#               start = start_values,
  #             control = nls.lm.control(maxiter = 1000))

###text start_value
# orginal data point
plot_data <- data.frame(x = x, y = y)
ggplot(plot_data, aes(x, y)) +
  geom_point(color = "steelblue") +
  #model
  geom_line(
    data = data.frame(
      x = seq(min(x), max(x), length.out = 100),
      y = start_values["Amax"] / (
        exp(-start_values["a"]*(seq(min(x), max(x), length.out = 100) - start_values["c"])) +
          exp(start_values["b"]*(seq(min(x), max(x), length.out = 100) - start_values["c"]))
      )
    ),
    aes(x, y),
    color = "red",
    linewidth = 1
  ) +
  labs(title = "Data vs Initial Model Fit")
###

s0 <- summary(model)
s1 <- s0$parameters
s2 <- confint(model, level=0.95)
results_table <- cbind(s1, s2)

# 将这张表写入 CSV (这会创建文件或覆盖旧文件)
write.csv(results_table, outSMRY)

#calc root mean square error and Nash Sutcliffe model efficiency
pred_y =  predict(model, x)
resid = y - pred_y
rmse = sqrt(mean(resid^2))
mean_y = mean(y)
nsme = 1 - sum((resid)^2)/sum((y-mean_y)^2)
seme <- rbind("rmse", rmse, "nsme", nsme)
write.table(seme, outSMRY, append=TRUE, sep=",")
#
library(investr)
nr = nrow(rawData)
newX <- data.frame(x = seq(0, rawData[nr,1], len=100))
fit <- investr::predFit(model, newX, interval = "confidence", level=0.95)
#fit <- investr::predFit(model, newX, interval = "prediction", level=0.95)
newY <- fit[,1]
lwr <- fit[,2]
upr <- fit[,3]
#
# plotting
obsData = data.frame(x, y)
newData = data.frame(newX, newY)
colnames(newData)[1] <- "newX"
xmax = max(x)
ymax = max(upr)
ggplot(data = newData, aes(newX, newY))+
  geom_line()+
  geom_point(data = obsData, aes(x, y))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill="blue", alpha=0.3)+
  scale_x_continuous(name="time (h)", breaks=seq(0, ceiling(xmax/10)*10, 5)) + 
  scale_y_continuous(name=paste0(yName," (g/L)", breaks=seq(0, ceiling(ymax/1)*1, 0.2))) +
  guides(x=guide_axis(minor.ticks=TRUE), y=guide_axis(minor.ticks=TRUE)) +
  theme(panel.background=element_rect(fill="white",color="grey50")) +
  annotate("text", x=0.7*max(x), y=0.8*max(y), 
           label=paste0("rmse = ", sprintf("%1.4f", rmse))) +
  annotate("text", x=0.7*max(x), y=0.7*max(y), 
           label=paste0("nsme = ", sprintf("%1.4f", nsme))) + 
  ggtitle("Frequentist") 
ggsave(paste0(fileNameRoot, "frequentist_", yName, "Profile.png"))
#
#result output
confInterval <- upr - lwr #conf.high - conf.low
out1 <- cbind(newX, newY, lwr, upr, confInterval) #predicted metabolite profile
colnames(out1)[1:2] <- c("time", "estimate")
p1 <- rep(x, length = nrow(out1))
p2 <- rep(y, length = nrow(out1))
nr = nrow(rawData) + 1
p11 <- replace(p1,nr:nrow(out1),0)
p21 <- replace(p2,nr:nrow(out1),0)
p12 <- cbind (p11, p21)
colnames(p12) <- c("time", yName)
output <- cbind(p12, out1)
write.csv(output, outCSV)

#

