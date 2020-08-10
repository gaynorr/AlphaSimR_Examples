#Read in results
sc1 = readRDS("Scenario1.rds")
sc2 = readRDS("Scenario2.rds")

#Plot results

#Inbred Mean
plot(-19:20,colMeans(sc2$inbredMean),type="l",col=2,
     main="Inbred Yield",xlab="Year",ylab="Yield (bu/ac)")
lines(-19:20,colMeans(sc1$inbredMean))

#Hybrid Mean
plot(-19:20,colMeans(sc2$hybridMean),type="l",col=2,
     main="Hybrid Yield",xlab="Year",ylab="Yield (bu/ac)")
lines(-19:20,colMeans(sc1$hybridMean))

#Inbred Variance
plot(-19:20,sqrt(colMeans(sc2$inbredVar)),type="l",col=2,
     main="Inbred Standard Deviation",xlab="Year",ylab="Yield (bu/ac)")
lines(-19:20,sqrt(colMeans(sc1$inbredVar)))

#Hybrid Variance
plot(-19:20,sqrt(colMeans(sc2$hybridVar)),type="l",col=2,
     main="Hybrid Standard Deviation",xlab="Year",ylab="Yield (bu/ac)")
lines(-19:20,sqrt(colMeans(sc1$hybridVar)))

#Hybrid Corr
plot(-19:20,sqrt(colMeans(sc2$hybridCorr)),type="l",col=2,
     main="GCA vs per se",xlab="Year",ylab="Correlation")
lines(-19:20,sqrt(colMeans(sc1$hybridCorr)))

#See how well long-term gain compares to published values
#Functions for matching regression lines in Troyer and Wellin (2009)
hybridYearFn = function(yield){ #bu/ac
  y = yield/0.01593 #kg/ha
  year = (y-4434.1)/74+1920
  return(year)
}
inbredYearFn = function(yield){ #bu/ac
  y = yield/0.01593 #kg/ha
  year = (y-1404.1)/48.3+1920
  return(year)
}

round(hybridYearFn(colMeans(sc2$hybridMean)),1)
round(inbredYearFn(colMeans(sc2$inbredMean)),1)
