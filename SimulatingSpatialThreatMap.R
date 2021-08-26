library(ggplot2)
library(rgdal)
library(rgeos)
library(raster)
library(reldist)

PlotSpatialPattern<-function(d)
{
  df <- data.frame(x = at.pts[,1], y = at.pts[,2],z = d)
  ggplot(df, aes(x = x, y = y, col = z)) +
    geom_point() +
    scale_colour_gradient(low="green", high="red")
}


#### Create a simulated autocorrelated surface for a region ####

ecoregion.grid<-readOGR("C:/Users/mikeha/Dropbox/Threat mapping/Regional/GlobalGrid50kmRoundl_Clip_Ecoregion.shp")

afrotropics.i<-which(ecoregion.grid$REALM == "Afrotropic")

afrotropic.grid<-ecoregion.grid[afrotropics.i,]
at.pts<-coordinates(gCentroid(afrotropic.grid,byid = T))
at.pts<-at.pts/50E3


Dd <- as.matrix(dist(at.pts))


## DEGREE OF AUTOCORRELATION
p = 0.05
##
for(p in c(1E-6,0.0001,0.05,0.3))
{
  # weights matrix
  w <- exp(-p * Dd)
  Ww <- chol(w)
  
  # errors
  z <- t(Ww) %*% rnorm(ncol(Ww),0,1)
  z<-scales::rescale(z,to = c(0,1))
  
  # plot
  df <- data.frame(x = at.pts[,1], y = at.pts[,2], z = z)
  
  
  ggplot(df, aes(x = x, y = y, col = z)) +
    geom_point() +
    scale_colour_gradient(low="green", high="red")
  
  
  #z.low.p<-z
  
  
  afrotropic.grid@data[[paste0("ThreatPrevalenceP",p)]]<-z
}


PlotSpatialPattern(afrotropic.grid$`ThreatPrevalenceP1e-06`)


#Save the simulated grid
save(afrotropic.grid, file = "C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Revised/afrotropic.threat.grid.Rdata")

#### Restart from here if have already saved the simulated grid####
load("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/afrotropic.threat.grid.Rdata")




#Read in the shapefile for the afrotropics grid, containing the threat intensity layers
at.pts<-coordinates(gCentroid(afrotropic.grid,byid = T))
at.pts<-at.pts/50E3


load("C:/Users/mikeha/Dropbox/Threat mapping/Regional 2017_3/Afrotropic/mammal_intersection.a.sp.R")
load("C:/Users/mikeha/Dropbox/Threat mapping/Regional 2017_3/Afrotropic/mammal_ranges.R")


access<-raster("F:/Work/Spatial data/accessibility_to_cities_2015_v1.0/accessibility_to_cities_2015_v1.0_moll.tif")
#access.moll<-projectRaster(from = access,crs = CRS(proj4string(afrotropic.grid)))

#hist(log(access))

#Simulate two different modes of sampling for threat assessment.
# 1. Perfect knowledge of the species ranges
# 2. Imperfect knowledge, approximated by the accessibility of the land surface

access.at.pts<-extract(access, coordinates(gCentroid(afrotropic.grid,byid = T)))


z.layers<-colnames(afrotropic.grid@data)[grep("ThreatPrevalence",colnames(afrotropic.grid@data))]

r.a<-r.mammal.a


q.prob<-c(0.25,0.5,0.75)

range.q.threat.intensity<-list()
range.access.weighted.q.threat.intensity<-list()
for(t in z.layers)
{
  z<-afrotropic.grid@data[[t]]
  weighted.q.threat.intensity<-array(NA, dim = c(ncol(r.a),length(q.prob)))
  q.threat.intensity<-array(NA, dim = c(ncol(r.a),length(q.prob)))
  for(sp in 1:ncol(r.a))
  {
    cell.i<-which((r.a[,sp]))
    
    #weighted.q.threat.intensity[sp]<-(log(access.pts[cell.i]) %*% z[cell.i])/sum(log(access.pts[cell.i]))
    #Want low travel times to be weighted higher than high travel times
    weighted.q.threat.intensity[sp,]<-wtd.quantile(z[cell.i][access.at.pts[cell.i] > 0], q=q.prob, na.rm = T, weight=log(1/access.at.pts[cell.i][access.at.pts[cell.i]>0]))
    q.threat.intensity[sp,]<-quantile(z[cell.i],probs = q.prob,na.rm = T)
  }
  range.access.weighted.q.threat.intensity[[t]]<-weighted.q.threat.intensity
  range.q.threat.intensity[[t]]<-q.threat.intensity
}


#Create 2 sets of threat codes for each threat prevalence pattern
#threat.breaks<-list(c(0,0.5,1.0),c(0,0.6,1.0),c(0,0.7,1.0))
draws<-10
ThreatCodes<-data.frame(binomial=r.mammal.ranges$binomial)
for(t in z.layers)
{
  #for(d in seq_len(draws))
  {
    for(q in q.prob)
    {
      
      temp = matrix(unlist(lapply(1:draws,function(x) runif(length(range.access.weighted.q.threat.intensity[[t]][,which(q.prob == q)])) <
            range.access.weighted.q.threat.intensity[[t]][,which(q.prob == q)])),
            ncol = draws,byrow = F)
      ThreatCodes[[paste0(t,".AccessWQ_",q)]]<-rowSums(temp) > draws/2
      
      temp = matrix(unlist(lapply(1:draws,function(x) runif(length(range.q.threat.intensity[[t]][,which(q.prob == q)])) <
                                    range.q.threat.intensity[[t]][,which(q.prob == q)])),
                    ncol = draws,byrow = F)
      ThreatCodes[[paste0(t,".Q_",q)]]<-rowSums(temp) > draws/2
    }
  }
}

write.csv(ThreatCodes,file<-"C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Revised/AftrotropicsThreatCodesConsensus10draws_w_p0.05.csv",row.names = F)

#afrotropic.grid@data = afrotropic.grid@data[,c(1,2,3,6,4)]

save(afrotropic.grid,file = "C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/afrotropic.threat.grid.Rdata")



#### Adding simulated threat using the scope data for birds ####
threat.codes<-read.csv("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/AftrotropicsThreatCodesConsensus10draws_w_p0.55.csv")


load("C:/Users/mikeha/Dropbox/Threat mapping/Regional 2017_3/Afrotropic/bird_intersection.a.sp.R")
load("C:/Users/mikeha/Dropbox/Threat mapping/Regional 2017_3/Afrotropic/bird_ranges.R")

draws<-10

weighted.scope<-list()
unweighted.scope<-list()
weighted.threat<-list()
unweighted.threat<-list()
for(t in z.layers)
{
  z<-afrotropic.grid@data[[t]]

  weighted.scope[[t]]<-rep(NA,ncol(r.bird.a))
  unweighted.scope[[t]]<-rep(NA,ncol(r.bird.a))
  weighted.threat[[t]]<-rep(NA,ncol(r.bird.a))
  unweighted.threat[[t]]<-rep(NA,ncol(r.bird.a))
  for(sp in 1:ncol(r.bird.a))
  {
    cell.i<-which(!is.na(r.bird.a[,sp]))
    #weighted.q.threat.intensity[sp]<-(log(access.pts[cell.i]) %*% z[cell.i])/sum(log(access.pts[cell.i]))
    cell.threat.m<-matrix(unlist(lapply(1:draws,function(x) runif(length(cell.i)) <
                           z[cell.i])),ncol = draws,byrow = F)
    if(sum(access.at.pts[cell.i]>0) > 1)
    {
      weighted.scope.draws<-2500*colSums(cell.threat.m[access.at.pts[cell.i]>0,] * 
        scales::rescale(log(1/access.at.pts[cell.i][access.at.pts[cell.i]>0]),to = c(0,1)))/r.bird.ranges$POLY_AREA[sp]
        #sum(access.at.pts[cell.i]>0)
    } else {
      weighted.scope.draws<-2500*(cell.threat.m[access.at.pts[cell.i]>0,] * scales::rescale(log(1/access.at.pts[cell.i][access.at.pts[cell.i]>0]),to = c(0,1)))/
      r.bird.ranges$POLY_AREA[sp]#sum(access.at.pts[cell.i]>0)
    }
    scope.draws<-colSums(cell.threat.m)/length(cell.i)#r.bird.ranges$POLY_AREA[sp]
    
    weighted.scope[[t]][sp]<-median(weighted.scope.draws, na.rm=T)
    weighted.threat[[t]][sp]<-sum(runif(draws) < 
                            weighted.mean(z[cell.i][(rowSums(cell.threat.m) > draws/2) & (access.at.pts[cell.i]>0)],
                                     weight=log(1/access.at.pts[cell.i][(rowSums(cell.threat.m) > draws/2) & (access.at.pts[cell.i]>0)]), na.rm=TRUE),
                            na.rm=T) > (draws/2)
    
    unweighted.scope[[t]][sp]<-median(scope.draws, na.rm=T)
    unweighted.threat[[t]][sp]<-sum(runif(draws) < median(z[cell.i][rowSums(cell.threat.m) > draws/2], na.rm = T)) > (draws/2)
    
    
  }
  
}



scope.mid.points<-list(
  "Minority (<50%)" = 0.25,
  "Majority (50-90%)" = 0.7,
  "Whole (>90%)" = 0.95,
  "Unknown" = 0.5
)

breaks<-c(-0.1,0.5,0.9,1.0)

weighted.scope.class<-lapply(weighted.scope,function(x) cut(x, breaks = breaks))
unweighted.scope.class<-lapply(unweighted.scope,function(x) cut(x, breaks = breaks))

for(m in names(weighted.scope.class)){
  weighted.scope.class[[m]][!weighted.threat[[m]]]<-"(0.9,1]"
  unweighted.scope.class[[m]][!unweighted.threat[[m]]]<-"(0.9,1]"
}

write.csv(weighted.threat,file="C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Aftrotropics_WeightedThreatCodesConsensus10draws_Scope.csv",row.names = F)
write.csv(unweighted.threat,file="C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Aftrotropics_UnweightedThreatCodesConsensus10draws_Scope.csv",row.names = F)

write.csv(weighted.scope.class,file="C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Aftrotropics_WeightedScopeConsensus10draws_Scope.csv",row.names = F)
write.csv(unweighted.scope.class,file="C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Aftrotropics_UnweightedScopeConsensus10draws_Scope.csv",row.names = F)
