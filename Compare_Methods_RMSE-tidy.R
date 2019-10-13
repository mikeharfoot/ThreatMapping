library(ggplot2)
library(rgdal)
library(rgeos)
library(SpatialPack)
library(RColorBrewer)

PlotSpatialPattern<-function(d,zlims,low.col,high.col)
{
  df <- data.frame(x = at.pts[,1], y = at.pts[,2],z = d)
  ggplot(df, aes(x = x, y = y, col = z)) +
    geom_point(size = 1) +
    scale_colour_gradient(limits = zlims,low=low.col, high=high.col)
}



# Script to compare analytical methods to see which best recreate the simulated spatial pattern of threat data
#
# 1. Threatened species richness	
# 2. Range weighted threatened species richness	
# 3. Proportion of species threatened	
# 4. Logistic regression with inverse range covariate	Test for different range transformations: e.g. log
# 5. Logistic regression weighted by inverse range	Test for different range transformations: e.g. log

#Read in the shapefile for the afrotropics grid, containing the threat intensity layers
load("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/afrotropic.threat.grid.Rdata")
at.pts<-coordinates(gCentroid(afrotropic.grid,byid = T))
at.pts<-at.pts/50E3



#Read the threat code data for each species
threat.codes<-read.csv("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/AftrotropicsThreatCodesConsensus5draws_w_p0.15.csv")

#Read mammal intersection and range data
load("C:/Users/mikeha/Dropbox/Threat mapping/Regional/Rdata/Afrotropic/mammal_intersection.a.sp.R")
load("C:/Users/mikeha/Dropbox/Threat mapping/Regional/Rdata/Afrotropic/mammal_ranges.R")


predicted.threats<-list()

threat.codes.to.model<-colnames(threat.codes)[-1]


#### 1. Threatened species richness   ####
predicted.threats[["threatened.richness"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["threatened.richness"]][[n]]<-apply(r.mammal.a,1,function(x) length(which(!is.na(x) & threat.codes[[n]] > 0)))
}


#### 2. Range weighted threatened species richness	####
predicted.threats[["range.weighted.threat.richness"]]<-list()
for(n in threat.codes.to.model)
{
  proportional.range<-apply(r.mammal.a,1,function(x) x/r.mammal.ranges$RangeArea)
  predicted.threats[["range.weighted.threat.richness"]][[n]]<-apply(proportional.range,2,function(x) sum(x[(!is.na(x) & threat.codes[[n]] > 0)]))
}


#### 3. Proportion of species threatened	####
predicted.threats[["proportion.threatened.richness"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["proportion.threatened.richness"]][[n]]<-predicted.threats[["threatened.richness"]][[n]]/apply(r.mammal.a,1,function(x) length(which(!is.na(x))))
}


new.data<-data.frame(pred = 1E4)

#### 4. Logistic regression with range covariate ####
predicted.threats[["logistic.range.size.covariate"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["logistic.range.size.covariate"]][[n]]<-apply(r.mammal.a,1,function(x) {
    pred<-(r.mammal.ranges$RangeArea[!is.na(x)])
    predict(glm(threat.codes[[n]][!is.na(x)] ~ pred,family = binomial(link = 'logit'),control = list(maxit = 50)),
            newdata = new.data,type = "response")
  })
}

#### 4a. Logistic regression with log range covariate ####
predicted.threats[["logistic.log.range.size.covariate"]]<-list()
for(n in threat.codes.to.model)
{
  
  predicted.threats[["logistic.log.range.size.covariate"]][[n]]<-apply(r.mammal.a,1,function(x) {
    pred<-log(r.mammal.ranges$RangeArea[!is.na(x)])
    predict(glm(threat.codes[[n]][!is.na(x)] ~ pred,family = binomial(link = 'logit'),control = list(maxit = 50)),
            newdata = log(new.data),type = "response")
  })
}

#### 4b. Logistic regression with sqrt range covariate ####
predicted.threats[["logistic.sqrt.range.size.covariate"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["logistic.sqrt.range.size.covariate"]][[n]]<-apply(r.mammal.a,1,function(x) {
    pred<-sqrt(r.mammal.ranges$RangeArea[!is.na(x)])
    predict(glm(threat.codes[[n]][!is.na(x)] ~ pred,family = binomial(link = 'logit'),control = list(maxit = 50)),
            newdata = sqrt(new.data),type = "response")
  })
}


#### 5. Logistic regression weighted by inverse range	 ####
predicted.threats[["logistic.range.size.weight"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["logistic.range.size.weight"]][[n]]<-apply(r.mammal.a,1,function(x) {
    median(predict(glm(threat.codes[[n]][!is.na(x)] ~ 1,weights = 1.0/(r.mammal.ranges$RangeArea[!is.na(x)]),family = binomial(link = 'logit')),type = "response"),na.rm=T)
  })
}


#### 5a. Logistic regression weighted by inverse log range####
predicted.threats[["logistic.log.range.size.weight"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["logistic.log.range.size.weight"]][[n]]<-apply(r.bird.a,1,function(x) {
    median(predict(glm(threat.codes[[n]][!is.na(x)] ~ 1,weights = 1.0/log(r.bird.ranges$POLY_AREA[!is.na(x)]),family = binomial(link = 'logit')),type = "response"),na.rm=T)
  })
}

#### 5b. Logistic regression weighted by inverse cube root range####
predicted.threats[["logistic.cbrt.range.size.weight"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["logistic.cbrt.range.size.weight"]][[n]]<-apply(r.bird.a,1,function(x) {
    median(predict(glm(threat.codes[[n]][!is.na(x)] ~ 1,weights = 1.0/(r.bird.ranges$POLY_AREA[!is.na(x)]^(1/3)),family = binomial(link = 'logit')),type = "response"),na.rm=T)
  })
}


#### 5c. Logistic regression weighted by inverse sqrt range####
predicted.threats[["logistic.sqrt.range.size.weight"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["logistic.sqrt.range.size.weight"]][[n]]<-apply(r.mammal.a,1,function(x) {
    median(predict(glm(threat.codes[[n]][!is.na(x)] ~ 1,weights = 1.0/sqrt(r.mammal.ranges$RangeArea[!is.na(x)]),family = binomial(link = 'logit')),type = "response"),na.rm=T)
  })
}

#### 5c. Logistic regression weighted by 1/2.5 range####
predicted.threats[["logistic.2.5rt.range.size.weight"]]<-list()
for(n in threat.codes.to.model)
{
  predicted.threats[["logistic.2.5rt.range.size.weight"]][[n]]<-apply(r.mammal.a,1,function(x) {
    median(predict(glm(threat.codes[[n]][!is.na(x)] ~ 1,weights = 1.0/(r.mammal.ranges$RangeArea[!is.na(x)]^(1/2.5)),family = binomial(link = 'logit')),type = "response"),na.rm=T)
  })
}



#### For bird scope data only ####
weighted.threat<-read.csv("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Aftrotropics_WeightedThreatCodesConsensus10draws_Scope.csv",stringsAsFactors = F)
unweighted.threat<-read.csv("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Aftrotropics_UnweightedThreatCodesConsensus10draws_Scope.csv",stringsAsFactors = F)

weighted.scope.class<-read.csv(file="C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Aftrotropics_WeightedScopeConsensus10draws_Scope.csv",stringsAsFactors =F)
unweighted.scope.class<-read.csv(file="C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/Aftrotropics_UnweightedScopeConsensus10draws_Scope.csv",stringsAsFactors =F)

scope.mid.points<-list(
  "Minority (<50%)" = 0.25,
  "Majority (50-90%)" = 0.7,
  "Whole (>90%)" = 0.95,
  "Unknown" = 0.5
)

predicted.threats.scope<-list()

#### 6a. Scope instead of range size weight - scope and threat weighted by accessibility####
predicted.threats.scope[["Logistic.scope.weighted"]]<-list()
for(n in names(weighted.threat))
{
  
  predicted.threats.scope[["Logistic.scope.weighted"]][[n]]<-apply(r.bird.a,1,function(x) {
    median(predict(glm(weighted.threat[[n]][!is.na(x)] ~ 1,
                       weights = unlist(lapply(weighted.scope.class[[n]][!is.na(x)], function(c) scope.mid.points[[as.numeric(c)]]))),
                   type = "response"),na.rm=T)
  })
}

#### 6b. Scope instead of range size weight - scope and threat unweighted by accessibility####
predicted.threats.scope[["Logistic.scope.unweighted"]]<-list()
for(n in names(weighted.threat))
{
  predicted.threats.scope[["Logistic.scope.unweighted"]][[n]]<-apply(r.bird.a,1,function(x) {
    median(predict(glm(unweighted.threat[[n]][!is.na(x)] ~ 1,
                       weights = unlist(lapply(unweighted.scope.class[[n]][!is.na(x)], function(c) scope.mid.points[[as.numeric(c)]]))),
                   type = "response"),na.rm=T)
  })
}


#### 6c. Logistic regression weighted by inverse cube root range but using threat code from scope simulation - weighted####
predicted.threats.scope[["Logistic.scope.weighted.cbrt.range.size.weight"]]<-list()
for(n in names(weighted.threat))
{
  predicted.threats.scope[["Logistic.scope.weighted.cbrt.range.size.weight"]][[n]]<-apply(r.bird.a,1,function(x) {
    median(predict(glm(weighted.threat[[n]][!is.na(x)] ~ 1,weights = 1.0/(r.mammal.ranges$RangeArea[!is.na(x)]^(1/3)),family = binomial(link = 'logit')),type = "response"),na.rm=T)
  })
}

#### 6d. Logistic regression weighted by inverse cube root range but using threat code from scope simulation - unweighted####
predicted.threats.scope[["Logistic.scope.unweighted.cbrt.range.size.weight"]]<-list()
for(n in names(weighted.threat))
{
  predicted.threats.scope[["Logistic.scope.unweighted.cbrt.range.size.weight"]][[n]]<-apply(r.bird.a,1,function(x) {
    median(predict(glm(unweighted.threat[[n]][!is.na(x)] ~ 1,weights = 1.0/(r.mammal.ranges$RangeArea[!is.na(x)]^(1/3)),family = binomial(link = 'logit')),type = "response"),na.rm=T)
  })
}



#### Calculate the correlation between threat prevalence and modelled threat #####

threat.map.types<-colnames(afrotropic.grid@data)
threat.map.types<-threat.map.types[grep("ThreatPrevalence",threat.map.types)]

threat.map.types.no.minus<-as.list(gsub(pattern = '-',replacement = '.',x = threat.map.types))
names(threat.map.types.no.minus)<-threat.map.types
threat.code.method<-paste0(c("Q_","AccessWQ_"),c(0.5,0.5,0.75,0.75))

u<-unique(unlist(lapply(threat.codes.to.model[grep(threat.map.types[1],threat.codes.to.model)],function(x) 
  strsplit(x,threat.map.types[1])[[1]][2])))
threat.code.method<-unique(sapply(u[order(u)],function(s) substr(s,start = 2,stop = nchar(s))))
#threat.code.method<-threat.code.method[-grep("0.6",threat.code.method)]


Frmse<-function(x,y)
{
  sqrt(mean((x-y)^2, na.rm=T))
}


# 1. Threatened species richness	
# 2. Range weighted threatened species richness	
# 3. Proportion of species threatened	
# 4. Logistic regression with inverse range covariate	Test for different range transformations: e.g. log
# 5. Logistic regression weighted by inverse range	Test for different range transformations: e.g. log

rmse<-list()
q.rmse<-list()
rmse.probs<-c(0,0.5,1.0)


for(tm in threat.map.types)
{
  tm.no.minus<-threat.map.types.no.minus[[tm]]
  
  rmse[[tm]]<-list()
  q.rmse[[tm]]<-list()
  for(tc in threat.code.method)
  {
    
    rmse[[tm]][[tc]]<-list()
    q.rmse[[tm]][[tc]]<-list()
    
    for(m in names(predicted.threats))
    {
      
      rmse[[tm]][[tc]][[m]]<-Frmse(afrotropic.grid@data[[tm]],predicted.threats[[m]][[paste0(tm,".",tc)]])

      q.rmse[[tm]][[tc]][[m]]<-quantile(unlist(rmse[[tm]][[tc]][[m]]),rmse.probs,na.rm=T)
    }
    
 
  }
  
}


median.rmse<-data.frame(row.names = names(predicted.threats)[-(1:2)])
rank.rmse<-data.frame(row.names = names(predicted.threats)[-(1:2)])
for(tm in threat.map.types[1:3])
{
  rmse.m<-data.frame(do.call(cbind,rmse[[tm]])[-(1:2),])
  median.rmse[[tm]]<-(apply(rmse.m,1,function(x) median(unlist(x), na.rm=T)))
  rank.rmse[[tm]]<-rank(apply(rmse.m,1,function(x) median(unlist(x), na.rm=T)))
}
write.csv(cbind(median.rmse,rowSums(rank.rmse)), file = "C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/MedianRMSE.csv")
rowSums(rank.rmse)


plot.reconstruction.methods<-names(predicted.threats)[c(3,8,9,11)]
plot.threat.code.methods<-list(1:3,4:6)
plot.threat.code.probs<-as.character(c(0.25,0.5,0.75))

plot.data<-list()

for(m in plot.reconstruction.methods)
{
  plot.data[[m]]<-array(NA, dim=c(2,9))
  for(tc in 1:length(plot.threat.code.methods))
  {
    cnt = 1
    tc.i<-plot.threat.code.methods[[tc]]
    for(tm in threat.map.types[c(1,3,2)])
    {
      for(i in 1:length(tc.i))
      {
        plot.data[[m]][tc,cnt]<-rmse[[tm]][[tc.i[i]]][[m]]
        cnt = cnt +1
      }
    }
  }
}
    

cols<-brewer.pal(n = 4,name = "Set1")
x.pos<-as.vector(sapply(1:3,function(x) x+c(-0.15,0,0.15)))

method.names<-c("Proportion of threatened species",
                "Logistic regression weighted by sqrt(R)",
                "Logistic regression weighted by cube root(R)",
                "Logistic regression weighted by ln(R)")

#Plot the summary RMSE for best model methods across different spatial threat maps and different
#simulated threat coding processes
png("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/RMSE_Best_Models.png",
    height = 12, width = 16, units = "cm",res = 600)
par(xpd = F)
plot(NA, xlim=c(0.5,3.5),ylim = c(0.1,0.4),xlab="",ylab = "",axes=F)
abline(v = x.pos,col= "grey",lty="dashed")
for(i in 1:length(plot.reconstruction.methods))
{
  points(jitter(x.pos,amount = 0.04),plot.data[[i]][1,],pch = 16,col = cols[i],cex=0.8)
  points(jitter(x.pos,amount = 0.04),plot.data[[i]][2,],pch = 17,col = cols[i],cex=0.8)
}

axis(side = 1,at = x.pos,labels = rep(c("25%","50%","75%"),3), cex.axis = 0.5,tick = F)
axis(side = 1,at = x.pos[c(2,5,8)],labels = c("Low SAC","Medium SAC","High SAC"), cex.axis = 0.8,tick = F,line = 2)
axis(side = 2,at = seq(0.1,0.4,0.05),las = 2,cex.axis = 0.6)
mtext(side = 2, "Square Root Mean Squared Error",line = 3,cex = 0.8)
par(xpd=T)
legend(x = 1, y = 0.445, legend = c("Perfect knowledge","Imperfect knowledge"),pch=c(16,17),bty = "n",cex=0.6)
legend(x = 2, y = 0.45, legend = method.names,pch=c(16),col = cols,bty="n",cex=0.6)
dev.off()


png("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/ThreatPatternP0.3.png",
    height = 20, width = 20, units = "cm",res = 600)
PlotSpatialPattern(afrotropic.grid@data$ThreatPrevalenceP0.3,c(0,1),"white","red")
dev.off()

png("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/ThreatPatternP0.05.png",
    height = 20, width = 20, units = "cm",res = 600)
PlotSpatialPattern(afrotropic.grid@data$ThreatPrevalenceP0.05,c(0,1),"white","red")
dev.off()

png("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/ThreatPatternP1E-4.png",
    height = 20, width = 20, units = "cm",res = 600)
PlotSpatialPattern(afrotropic.grid@data$`ThreatPrevalenceP1e-04`,c(0,1),"white","red")
dev.off()



#Mean RMSE for each model method
lapply(plot.data,function(x) apply(x,1,median))
lapply(plot.data,median)

# $proportion.threatened.richness
# [1] 0.1521534 0.2095129
# 
# $logistic.log.range.size.weight
# [1] 0.1504516 0.2083537
# 
# $logistic.cbrt.range.size.weight
# [1] 0.1434960 0.1907036
# 
# $logistic.2.5rt.range.size.weight
# [1] 0.1437281 0.1917147



png("C:/Users/mikeha/Dropbox/Threat mapping/SimulatedThreatData/CubrtvsSimulated.png",
    height = 12, width = 12, units = "cm",res = 600)

plot(afrotropic.grid$`ThreatPrevalenceP1e-04`,
     predicted.threats$logistic.cbrt.range.size.weight$`ThreatPrevalenceP1e-04.AccessWQ_0.5`,pch=16,
     xlim = c(0,1),ylim =c(0,1), col = rgb(0,0,0,0.2), cex=0.8, xlab="Simulated threat",ylab="Threat modelled using mammal ranges",
     las=1)
abline(a = 0,b =1,col = "blue")
dev.off()

cor.test(afrotropic.grid$`ThreatPrevalenceP1e-04`,
              predicted.threats$logistic.cbrt.range.size.weight$`ThreatPrevalenceP1e-04.AccessWQ_0.5`)



##### Compare range size versus scope weighted models ####

rmse.scope<-list()
q.rmse.scope<-list()

for(tm in names(weighted.scope))
{
  for(m in names(predicted.threats.scope))
{
    
    rmse.scope[[tm]][[m]]<-Frmse(afrotropic.grid@data[[tm]],predicted.threats.scope[[m]][[tm]])
    q.rmse.scope[[tm]][[m]]<-quantile(unlist(rmse.scope[[tm]][[m]]),rmse.probs,na.rm=T)
    
  }
  
}

plot(afrotropic.grid@data$ThreatPrevalenceP0.3,predicted.threats.scope$Logistic.scope.weighted$ThreatPrevalenceP0.3,pch=16,cex=0.5,ylim=c(0,1),xlim=c(0,1))
points(afrotropic.grid@data$ThreatPrevalenceP0.3,predicted.threats.scope$Logistic.scope.weighted.cbrt.range.size.weight$ThreatPrevalenceP0.3,
       pch=16,cex=0.5,col="red")
abline(a = 0, b = 1,col="blue")

tm<-"ThreatPrevalenceP1e-04"
plot(afrotropic.grid@data[[tm]],predicted.threats.scope$Logistic.scope.weighted[[tm]],pch=16,cex=0.5,ylim=c(0,1),xlim=c(0,1))
points(afrotropic.grid@data[[tm]],predicted.threats.scope$Logistic.scope.weighted.cbrt.range.size.weight[[tm]],
       pch=16,cex=0.5,col="red")
abline(a = 0, b = 1,col="blue")


plot(predicted.threats.scope$Logistic.scope.weighted[[tm]],predicted.threats.scope$Logistic.scope.weighted.cbrt.range.size.weight[[tm]],
     pch=16,cex=0.5,ylim=c(0,1),xlim=c(0,1))


plot(afrotropic.grid@data$`ThreatPrevalenceP1e-04`,
     rowMeans(data.frame(predicted.threats.scope$Logistic.scope.weighted$`ThreatPrevalenceP1e-04`,
                         predicted.threats.scope$Logistic.scope.weighted.cbrt.range.size.weight$`ThreatPrevalenceP1e-04`)),
     pch=16,cex=0.5,ylim=c(0,1),xlim=c(0,1))
abline(a = 0, b = 1,col="blue")



lapply(rmse.scope,function(x) x[c(1,3)])

mammal.predicted.threats.scope<-predicted.threats.scope
mammal.rmse.scope<-rmse.scope


tm<-"ThreatPrevalenceP0.05"
plot(afrotropic.grid@data[[tm]],mammal.predicted.threats.scope$Logistic.scope.weighted[[tm]],pch=16,cex=0.5,ylim=c(0,1),xlim=c(0,1))
points(afrotropic.grid@data[[tm]],mammal.predicted.threats.scope$Logistic.scope.weighted.cbrt.range.size.weight[[tm]],
       pch=16,cex=0.5,col="red")
abline(a = 0, b = 1,col="blue")




plot(afrotropic.grid@data$`ThreatPrevalenceP1e-04`,predicted.threats$logistic.log.range.size.weight$`ThreatPrevalenceP1e-04.Q_0.5`,pch=16,cex=0.5,ylim=c(0,1),xlim=c(0,1))
points(afrotropic.grid@data$`ThreatPrevalenceP1e-04`,predicted.threats$logistic.cbrt.range.size.weight$`ThreatPrevalenceP1e-04.Q_0.75`,pch=16,cex=0.5,col = "red")
abline(a = 0, b = 1,col="blue")
