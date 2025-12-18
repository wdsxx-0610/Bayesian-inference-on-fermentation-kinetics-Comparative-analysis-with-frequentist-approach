# Jags-Ymet-Xnom1grp-Mnormal.R 
# Accompanies the book:
#   

source("../utils/DBDA2E-utilities.R")

#===============================================================================

genMCMC = function(x, y, 
                   numSavedSteps=numSavedSteps, 
                   thinSteps=thinSteps ,saveName=fileNameRoot ) { 
  #require(rjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  #x = data[,xName]  #time
  #y = as.matrix(data[,yName],ncol=length(yName))  #all metabolites
  #x = data[,xName]
  #y = data[,yName]
  # Specify data in list, for later shipment to JAGS:
  dataList = list(x = x, y = y, Ntotal = length(y))
  print(dataList)#-----------------------------------------------------------------------------
  # THE MODEL.
  #===============================================================================
  modelString = "
  data {
    #(1)prior for acetate
    #inits <- c(4.0184, -1.7612, -5.7387, 3.2522, -0.012)
    #inits_precision <-c(46.6694/4, 27.2778/4, 1.4055/4, 120.2169/4, 1.4442)
    #(2)prior for ethanol
    #inits <- c(0.4095, -1.0635, -5.5425, 1.8806, -2.3111)
    #inits_precision <- c(0.3211,0.8844,0.2121,0.6064,0.3676)
    #(3)prior for lactate
    inits <- c(2.7348, -1.6316, -1.7737, 3.0193, -0.6804)
    inits_precision <-c(52.1589, 4.2834, 4.4571, 38.8296, 3.1878)
    #(4)prior for PDO
    #inits <- c(3.9939, -1.2502, -6.0002, 3.4651, 0.1206)
    #inits_precision <-c(25.3153, 12.7001, 0.7765, 202.8906, 0.9617)
    ##old one
    #inits = c(26.6305, 0.1129, 0.2156, 25.6914, 0.001) # initial guess for Lactate
    #inits = c(49.4012, 0.1359, 0.0001, 21.3734, 0.001) # initial guess for Acetate
    #inits = c( 6.7617, 0.1598, 0.0001, 26.0465, 0.001) # initial guess for Ethanol
    #inits = c(87.5607, 0.1359, 0.0018, 27.4886, 0.001) # initial guess for PDO


    init_Amax = inits[1] 
    init_a = inits[2] 
    init_b = inits[3] 
    init_c = inits[4]  
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
    c ~ dlnorm(init_c, inits_precision[4] )  
    real_sigma ~ dunif(0, 100)
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
  initsList = list( Amax = init_Amax, a = init_a, b = init_b, c = init_c)
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  require(rjags)
  require(runjags)
  #parameters = c( "mu" , "sigma")     # The parameters to be monitored
  #parameters = c( "Amax", "a" , "b" , "c", "sigma", "nu")     # The parameters to be monitored
  parameters = c( "Amax", "a" , "b" , "c", "sigma")     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 4 
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
                     compValAmax=NULL , ropeAmax=NULL , 
                     compVala=NULL , ropea=NULL , 
                     compValb=NULL , ropeb=NULL , 
                     compValc=NULL , ropec=NULL , 
                     compValSigma=NULL , ropeSigma=NULL , 
                     saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "Amax" = summarizePost( mcmcMat[,"Amax"] , 
                                                compVal=compValAmax , 
                                                ROPE=ropeAmax ) )
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
  #                     "nu" = summarizePost( mcmcMat[,"nu"] , 
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
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( Amax,a,b,c,sigma)[plotIdx,] ,
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

