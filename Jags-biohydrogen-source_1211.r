# Jags-Ymet-Xnom1grp-Mnormal.R 
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("/Users/wdsxx0610/Documents/R_directory/A_formation/biohydrogen/DBDA2E-utilities.R")

#===============================================================================

genMCMC = function(x, y, 
                   numSavedSteps=numSavedSteps, thinSteps=thinSteps, 
                   nChains=nChains, saveName=fileNameRoot) { 
  #require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  #x = data[,xName]  #time
  #y = as.matrix(data[,yName],ncol=length(yName))  #all metabolites
  #x = data[,xName]
  #y = data[,yName]
  # Specify data in list, for later shipment to JAGS:
  dataList = list(x = x, y = y, Ntotal = length(y))
  #-----------------------------------------------------------------------------
  # THE MODEL.
  #===============================================================================
  modelString = "
  data {
  
    #(1)prior for biomass
    #inits <- c(0.9859,-1.2918,-5.5327,2.0691,-2.5359)
    #inits_precision <-c(36.2894/4,3.4581/4,1.4058/4,7.5521/4,10.0299/4)
    #(2)prior for glucose
    #inits <- c(2.5767,-4.8622,-2.3767,2.6986,-1.8389)
    #inits_precision <-c(24.011,1.2993,29.7431,8.0794,10.2276)
    #(3)prior for hydrogen
    #inits <- c(1.2689,-1.905,-5.2328,3.1931,-2.1284)
    #inits_precision <-c(46.2536,47.8709,1.5713,92.5573,9.7575)
    #(4)prior for acetate
    #inits <- c(4.6352,-2.2183,-3.7177,3.3933,0.6823)
    #inits_precision <-c(16.2087,11.9845,1.6677,18.021,6.0533)
    #(5)prior for lactate
    inits <- c(3.8568,-2.2254,-3.8923,3.3163,0.3151)
    inits_precision <-c(8.9396,5.9911,1.1023,8.5691,5.017)
    #old value
    #inits <- c(  1.1040, 0.1366, 0.0623, 22.0649, 0.001) #initial guess for biomass
    #inits <- c(10.0000, 0.0001, 0.1859,  0.0001, 0.001)  #initial guess for glucose
    #inits <- c( 3.1429, 0.1598, 0.0001, 52.0000, 0.001) #initial guess for hydrogen
    #inits <- c(58.1180, 5.5815, 0.0001, 52.0000, 0.001) #initial guess for acetate
    #inits <- c(44.1830, 1.9863, 0.0001, 52.0000, 0.001) #initial guess for lactate
    #
    init_Amax = inits[1]
    init_a = inits[2] 
    init_b = inits[3] 
    init_c = inits[4]  
    init_s = inits[5]
  }
  model {
    for ( i in 1:Ntotal ) {
      #y[i] ~ dt(mu[i], sigma, nu)
      y[i] ~ dnorm(mu[i], sigma)
      mu[i] <- Amax/(exp(-a*(x[i]-c))+exp(b*(x[i]-c)))
    }
    Amax ~ dlnorm(init_Amax, inits_precision[1] )
    a ~ dlnorm(init_a, inits_precision[2] )
    b ~ dlnorm(init_b, inits_precision[3])
    c ~ dlnorm(init_c, inits_precision[4])    
    real_sigma ~ dunif(init_s , 100)
    sigma = 1/(real_sigma*real_sigma)
    #nu ~ dexp(1/30.0)
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )

  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  #mu = mean(y) 
  #sigma = sd(y) 
  #initsList = list( mu = mu , sigma = sigma )
  initsList = list(Amax = init_Amax, a = init_a, b = init_b, c = init_c)
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  require(rjags)
  require(runjags)
  #parameters = c( "mu" , "sigma")     # The parameters to be monitored
  #parameters = c("Amax", "a" , "b" , "c", "sigma", "nu")     # The parameters to be monitored
  parameters = c("Amax", "a" , "b" , "c", "sigma")     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 6 
  thinSteps = 2
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                          n.chains=nChains , n.adapt=adaptSteps )
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nIter , thin=thinSteps )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste0(saveName,"Mcmc.Rdata") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function( codaSamples , 
                     compVala=NULL , ropea=NULL , 
                     compValb=NULL , ropeb=NULL , 
                     compValc=NULL , ropec=NULL , 
                     compValSigma=NULL , ropeSigma=NULL , 
                     saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "Amax" = summarizePost( mcmcMat[,"Amax"] , 
                                            compVal=compVala , 
                                            ROPE=ropea ) )
  summaryInfo = rbind( summaryInfo , 
                       "a" = summarizePost( mcmcMat[,"a"] , 
                                               compVal=compVala , 
                                               ROPE=ropea ) )
  summaryInfo = rbind( summaryInfo , 
                       "b" = summarizePost( mcmcMat[,"b"] , 
                                            compVal=compValb , 
                                            ROPE=ropeb ) )
  summaryInfo = rbind( summaryInfo , 
                       "c" = summarizePost( mcmcMat[,"c"] , 
                                            compVal=compValc, 
                                            ROPE=ropec ) )
   summaryInfo = rbind( summaryInfo , 
                       "sigma" = summarizePost( mcmcMat[,"sigma"] , 
                                                compVal=compValSigma , 
                                                ROPE=ropeSigma ) )  
  #summaryInfo = rbind( summaryInfo , 
  #                    "nu" = summarizePost( mcmcMat[,"nu"] , 
  #                                           compVal=NULL , 
  #                                           ROPE=NULL) )
  #summaryInfo = rbind( summaryInfo , 
  #                     "log10(nu)" = summarizePost( log10(mcmcMat[,"nu"]) , 
  #                                                  compVal=NULL , ROPE=NULL ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste0(saveName,"bayesian_summary.csv") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , xName="x" , yName="y" ,
                     showCurve=FALSE , pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  Amax = mcmcMat[,"Amax"]
  a = mcmcMat[,"a"]
  b = mcmcMat[,"b"]
  c = mcmcMat[,"c"]
  sigma = mcmcMat[,"sigma"]
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph(width=7*5/5,height=7*5/5) # 5 parameters to be paired
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      #usr = c(0,1,0,1)
      #par(usr)
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind(Amax,a,b,c,sigma)[plotIdx,] ,
           labels=c( expression(Amax),
                     expression(a),
                     expression(b),
                     expression(c),
                     expression(sigma)), 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,"PostPairs"), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
}

#===============================================================================

