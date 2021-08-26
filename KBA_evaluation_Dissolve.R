wants <- c("rgeos","rgdal","maptools","plyr","raster","RColorBrewer","viridis","ash","spdep","reshape2")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
has <- wants %in% .packages()
for(p in wants[!has]) library(p,character.only = T)

out.dir = 'C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/'
path = 'C:/Users/mikeha/Dropbox/Threat mapping/'
source(paste0(path,'scripts/ThreatMapFunctions.R'))

wgs1984.proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
moll.proj <- CRS("+proj=moll +datum=WGS84")

grid<-readOGR("C:/Users/mikeha/Dropbox/Threat mapping/Regional","GlobalGrid50kmRoundl_Clip_Ecoregion")

#Read in map of the world
WL <- readShapePoly("C:/Users/mikeha/Dropbox/Threat mapping/Spatial data/ne_10m_land.shp")#download from 'naturalearth.com'
proj4string(WL) <- wgs1984.proj  
WL.moll <- spTransform(WL, moll.proj)



realms = c(
  'Afrotropic',
  'Indomalayan',
  'Nearctic',
  'Neotropic',
  'Palearctic'
)

regions.comb = list(
  Australasia_Oceania = c('Australasia','Oceania')
) 

taxa = c('mammal','bird','amph')

ThreatCodes = GetThreatCodesPhases()

phase = 'phase1'
th = ThreatCodes[[phase]][4]

models<-c("cub.rt.R_W_Pred_","PropThSp")#,"sqrtR_W_Pred_"
model.var.names<-c("p.r.cube.rt","p.th")

path = 'C:/Users/mikeha/Dropbox/Threat mapping/'


g.q<-list()


for(i in 1:length(models))
{
  print(models[i])
  for(th.i in 1:length(ThreatCodes[[phase]]))
  {
    th <- ThreatCodes[[phase]][th.i]
    print(th)
    for(t in taxa)
    {
      print(t)
      p = c(0,0.025,0.25,0.5,0.75,0.975,1.0)
      g.q[[model.var.names[i]]][[t]][[names(th)]]<-GetPredictionQuantiles(path = path,
                                                             model = models[i],model.name = model.var.names[i],t = t,th = th,
                                                             realms = realms,regions.comb = regions.comb,p = p)
    }
  }
}


#Remove P(Th) predictions where N species < 10
for(i in 1:length(models))
{
  print(models[i])
  for(th.i in 1:length(ThreatCodes[[phase]]))
  {
    th <- ThreatCodes[[phase]][th.i]
    print(th)
    for(t in taxa)
    {
      print(t)
      lt.10sp<-which(g.q[[model.var.names[i]]][[t]][[names(th)]]$global.sr < 10)
      g.q[[model.var.names[i]]][[t]][[names(th)]]$globalq[lt.10sp]<-NA
      
    }
  }
}


color.schemes<-list(
  "Logging" = brewer.pal(8,"Greens"),
  "Agriculture" = rev(brewer.pal(11,"BrBG")[1:6]),
  "Hunting" = brewer.pal(8,"Reds"),
  "Pollution" = brewer.pal(8,"Blues"),
  "Invasives" = brewer.pal(8,"Purples"),
  "Climate change" = brewer.pal(8,"Oranges")
)


primary.threats<-names(color.schemes) 


### Print correlation plots from the two models ####

png("C:/Users/mikeha/Dropbox/Threat mapping/KBA/ProportionThreatened vs Inv cube root weighted.png",unit="cm",height = 17, width=15,
    res = 600)
cor.results<-list()
taxa.label <- list("Mammals", "Birds","Amphibians")
names(taxa.label) = taxa
layout(matrix(1:28, ncol=4,byrow=F), widths = c(0.1,0.3,0.3,0.3),heights = c(rep(0.15,6),0.1))
par(mar = c(0.25,0.25,0.25,0.25))
for(i in 1:3){
  plot(NA, axes=F, xlab="",ylab= "", xlim=c(0,1),ylim = c(0,1))
  mtext(primary.threats[i],side = 2, line = -2, cex=0.6)
}
mtext("Threat probability from inverse cube root weighted proportion model", side= 2,line = -1, cex=0.8)
for(i in 4:7){
  plot(NA, axes=F, xlab="",ylab= "", xlim=c(0,1),ylim = c(0,1))
  mtext(primary.threats[i],side = 2, line = -2, cex=0.6)
}
for(t in taxa)
{
  
  for(th in primary.threats)
  {
    plot(g.q$p.th[[t]][[th]]$globalq,
         g.q$p.r.cube.rt[[t]][[th]]$globalq,pch=16,cex=0.2,col=rgb(0,0,0),
         xlim=c(0,1), ylim = c(0,1), xlab="",ylab="",xaxt='n',yaxt=ifelse(t == taxa[1],'s','n'), las = 1)
    
    if(th == primary.threats[1]) text(x = 0.1,y = 0.9, taxa.label[t],cex = 0.8)
    axis(1,labels = F)
    axis(2,labels=F)
    abline(a = 0, b = 1, col="blue")
    cor.results[[t]][[th]]<-cor.test(g.q$p.th[[t]][[th]]$globalq,
                                     g.q$p.r.cube.rt[[t]][[th]]$globalq)
    
  }
  axis(1,tick=F)
  if(t == taxa[2]) mtext("Threat probability from proportion of threatened species model", side= 1,line = 2.5, cex=0.8)
  plot.new()
}

dev.off()

do.call(rbind,lapply(cor.results, function(c.t) lapply(c.t, function(c.t.th) c.t.th$estimate^2)))
do.call(rbind,lapply(cor.results, function(c.t) lapply(c.t, function(c.t.th) c.t.th$p.value)))


#### Calculate threat probabilities in PAS ####
dsn<-"C:/Users/mikeha/Documents/ArcGIS/Projects/ThreatMapIntersection/ThreatMapIntersection.gdb"
ogrListLayers(dsn)
kba<-readOGR(dsn=dsn,layer="KBAGlobal2019_September02_moll")

kba.data<-read.csv(paste0(path,"KBA/Threats data downloaded 10 Jan 2020.csv"),stringsAsFactors = F)
kba.threat.match<-read.csv(paste0(path,"KBA/Threat_categories_KBA_match.csv"),stringsAsFactors = F)

#Add threat codes to the kba.data
kba.data$ThreatLevelCode<-kba.threat.match$ThreatLevel2.1[match(kba.data$ThreatLevel2,kba.threat.match$ThreatLevel2)]

#Exclude future threat
kba.data.now<-kba.data[kba.data$TimingScore == 3,]

#transform from long format to wide
kba.impact.scores.wide <- dcast(kba.data.now, SiteID + Country + NationalName + InternationalName ~ ThreatLevelCode, value.var="ImpactScore",fun.aggregate = max)
kba.scope.scores.wide <- dcast(kba.data.now, SiteID + Country + NationalName + InternationalName ~ ThreatLevelCode, value.var="ScopeScore",fun.aggregate = max)
kba.severity.scores.wide <- dcast(kba.data.now, SiteID + Country + NationalName + InternationalName ~ ThreatLevelCode, value.var="SeverityScore",fun.aggregate = max)

for( c in grep("^[[:digit:]]+", colnames(kba.impact.scores.wide)))
{
  kba.impact.scores.wide[is.infinite(kba.impact.scores.wide[,c]),c]<-NA
  kba.scope.scores.wide[is.infinite(kba.scope.scores.wide[,c]),c]<-NA
  kba.severity.scores.wide[is.infinite(kba.severity.scores.wide[,c]),c]<-NA
}


kba.severity.score.levels<-list(
  "no" = 0,
  "slow" = 1,
  "moderate" = 2,
  "very rapid" = 3
)


#kba<-spTransform(kba,moll.proj)

grid.kba<-gIntersects(kba,grid,byid = T,returnDense = F)

kba.areas<-gArea(kba,T)

#length(which(kba.areas > 0.5*(50E3)^2))


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



probs<-c(0,0.025,0.25,0.5,0.75,0.975,1.0)

#mett.pa.realm<-apply(small.grid.mett.wdpa,2,function(x) getmode(small.grid.data$REALM[x]))

kba.p.th<-list()
for(i in 1:length(models))
{
  print(models[i])
  for(th in primary.threats)
  {
    print(th)
    
    #loop over taxanomic group
    for(t in taxa)
    {
      print(t)
      #find the data column
      
      p.th<-g.q[[model.var.names[i]]][[t]][[th]]$globalq
      
      #find the median P(th) for each mett pa
      kba.p.th[[model.var.names[i]]][[t]][[th]]<-sapply(grid.kba,function(x) 
        quantile(unlist(p.th[match(x,g.q[[model.var.names[i]]][[t]][[th]]$global.i)]),probs = probs, na.rm=T))
        
      kba.p.th[[model.var.names[i]]][[t]][[th]][is.infinite(kba.p.th[[model.var.names[i]]][[th]][[t]])]<-NA
    }
  }
}



Translate.KBAid.PTh<-function(kba.ids,pth)
{
  
  #first match the kba.ids to kba.grid intersection rows and find the grid cell 
  grid.ids.to.use<-unique(unlist(sapply(grid.kba[match(kba.ids,kba$SitRecID)], function(x)
         match(x,g.q[[m]][[t]][[th]]$global.i))))
  
  r<-unlist(pth[grid.ids.to.use])
  if(sum(is.na(r)) > 0) r<-r[-which(is.na(r))]
  r
  
}



use.quant<-which(probs == 0.5)


#Now aggregate p(th) by mett pa category 
w.gt<-list() #holds wilcox.test results
w.2sided<-list()
q.p.th.kba.level<-list()
q.p.th.kba.group<-list()
mn.p.th.kba.level<-list()
sd.p.th.kba.level<-list()
count.kbas.level<-list()
p.th.kba.level<-list()
p.th.kba.group<-list()
kendall.pth.kba<-list()
for(th in primary.threats)
{
  print(th)
  #first find the maximum category for each mett pa
  #find the mett columns for this threat
  kba.th.cols<-grep(paste0("^",ThreatCodes[[phase]][[th]]),colnames(kba.severity.scores.wide),value = T)
  if(length(kba.th.cols) > 1)
  {
    # agg.mett.th.levels<-apply(mett.w.wdpa.coded[,mett.th.cols],1,function(x) {
    #   u<-unique(x)
    #   if(length(u) > 1)
    #     NA
    #   else
    #     u
    #   })
    
    # If there is more than one threat code - find the max threat code for each kba across the threats
    agg.kba.th.levels<-apply(kba.severity.scores.wide[,kba.th.cols],1,function(x) min(x,na.rm=T) )
  } else
  {
    #else take the threat codes for the single threat
    agg.kba.th.levels<-kba.severity.scores.wide[,kba.th.cols]
  }
  
  #List the kbas falling into each threat code
  kba.IDs.by.level<-lapply(kba.severity.score.levels, #seq(-3,0,1)
                            function(m) kba.severity.scores.wide$SiteID[which(agg.kba.th.levels %in% m)])
  names(kba.IDs.by.level)<-names(kba.severity.score.levels)#rev(names(kba.impact.score.levels))
  
  kba.areas.by.level<-do.call(rbind,
                              lapply(kba.IDs.by.level,function(c) 
                                quantile(kba.areas[match(c,kba$SitRecID)],probs=probs, na.rm=T)))
  
  count.kbas.level[[th]]<-lapply(kba.IDs.by.level, length)
  
  for(m in model.var.names)
  {
    print(m)
    for(t in taxa)
    {
      print(t)
      
    
      q.p.th.kba.level[[m]][[th]][[t]]<-do.call(rbind,
            lapply(kba.IDs.by.level,function(c) 
              quantile(Translate.KBAid.PTh(c,
                                  g.q[[m]][[t]][[th]]$globalq), probs = probs, na.rm = T)))
      

      mn.p.th.kba.level[[m]][[th]][[t]]<-do.call(rbind,
                                                 lapply(kba.IDs.by.level,function(c) 
                                                   mean(Translate.KBAid.PTh(c,
                                                                                g.q[[m]][[t]][[th]]$globalq), na.rm = T)))
      
      sd.p.th.kba.level[[m]][[th]][[t]]<-do.call(rbind,
                                                 lapply(kba.IDs.by.level,function(c) 
                                                   sd(Translate.KBAid.PTh(c,
                                                                            g.q[[m]][[t]][[th]]$globalq), na.rm = T)))
      
      p.th.kba.level[[m]][[t]][[th]]<-lapply(kba.IDs.by.level,function(c) 
        (Translate.KBAid.PTh(c,g.q[[m]][[t]][[th]]$globalq)))
        
      
      p.th.kba.group[[m]][[t]][[th]][["low"]]<-Translate.KBAid.PTh(c(kba.IDs.by.level$no,kba.IDs.by.level$slow),
                                                                     g.q[[m]][[t]][[th]]$globalq)
      p.th.kba.group[[m]][[t]][[th]][["high"]]<-Translate.KBAid.PTh(c(kba.IDs.by.level$moderate,kba.IDs.by.level$`very rapid`),
                                                                   g.q[[m]][[t]][[th]]$globalq)
      
      q.p.th.kba.group[[m]][[th]][[t]]<-do.call(rbind,
                                                lapply(p.th.kba.group[[m]][[t]][[th]],function(c) 
                                                  quantile(unlist(c), probs = probs, na.rm = T)))
        
      kendall.pth.kba[[m]][[t]][[th]]<-cor.test(unlist(p.th.kba.level[[m]][[t]][[th]]),
               unlist(sapply(1:4,function(i) rep(i,lapply(p.th.kba.level[[m]][[t]][[th]],length)[[i]]))),
               method = "kendall")
      
      low<-c(p.th.kba.level[[m]][[t]][[th]]$no,p.th.kba.level[[m]][[t]][[th]]$slow)
      high<-c(p.th.kba.level[[m]][[t]][[th]]$moderate,p.th.kba.level[[m]][[t]][[th]]$`very rapid`)
      
      w.gt[[m]][[th]][[t]]<-wilcox.test(x = high,
                                   y = low,
                                   exact = F,
                                   alternative = "great")
      w.2sided[[m]][[th]][[t]]<-wilcox.test(x = high,
                                       y = low,
                                       exact = F,alternative = "two.sided")
      
    }
  }
  
}


get.signif<-function(d)
{
  code = ""
  if(d < 0.01)
    code = "***"
  else if(d < 0.05)
    code = "**"
  else if(d < 0.1)
    code = "*"
  code
}


cols<-list("red","green","blue")
names(cols)<-taxa

x.pos<-list(-0.25,0,0.25)
x.pos.adj<-c(-0.045,-0.015,0.015,0.045)
names(x.pos)<-taxa

png("C:/Users/mikeha/Dropbox/Threat mapping/KBA/KendallRankCorrelation_Dissolve_PTh.png",unit="cm",height = 15, width=15,
    res = 600)
par(mar=c(6,5,1,2))
plot(NA, xlim = c(0.75,length(primary.threats)+0.2),ylim=c(-0.1,0.3),axes=F,xlab="",ylab="")
for(i in 1:length(primary.threats))
{
  th = primary.threats[i]
  print(th)
  y.cubrt<-do.call(rbind,lapply(kendall.pth.kba$p.r.cube.rt, function(d) c(d[[th]]$estimate,d[[th]]$p.value)))
  y.p.th<-do.call(rbind,lapply(kendall.pth.kba$p.th, function(d) c(d[[th]]$estimate,d[[th]]$p.value)))
  
  points(i+unlist(x.pos),
         y.cubrt[,1],
         pch=16,col = unlist(cols))
  points(i+unlist(x.pos),
         y.p.th[,1],
         pch=17,col = unlist(cols))

  print(unlist(lapply(q.p.th.kba.level$p.r.cube.rt[[th]], function(d) d[4,4]/d[1,4]))/
          unlist(lapply(q.p.th.kba.level$p.th[[th]], function(d) d[4,4]/d[1,4])))
  signif.p.th<-sapply(y.cubrt[,2],get.signif)
  signif.cubrt<-sapply(y.p.th[,2],get.signif)
  text(i+unlist(x.pos),y.p.th[,1],signif.p.th,pos = 2,col= "grey", cex=0.8)
  text(i+unlist(x.pos),y.cubrt[,1],signif.cubrt,pos = 2,col= "black", cex=0.8)
  
}
abline(h = 0,col = "grey")
abline(v=seq(1.5,length(primary.threats)),lty="dashed")
axis(side = 1,at=seq(1,length(primary.threats)),labels = primary.threats,las=2,tick = F,cex.axis=0.8)
axis(side =2,las=1,cex.axis=0.8)

mtext(side = 2,text = "Kendall tau",cex=0.8,line = 3)
legend("topright",legend = c("Mammals","Birds","Amphibians"),pch=rep(16,3), col = unlist(cols),bty='n')
dev.off()

quantile(unlist(lapply(kendall.pth.kba$p.r.cube.rt,function(d) lapply(d, function(th) th$estimate))),probs=c(0,0.5,1))
quantile(unlist(lapply(kendall.pth.kba$p.th,function(d) lapply(d, function(th) th$estimate))),probs=c(0,0.5,1))






x.pos.adj<-c(-0.09,-0.03,0.03,0.09)
names(x.pos)<-taxa

png("C:/Users/mikeha/Dropbox/Threat mapping/KBA/PImpactInKBAClassesDissolve.png",unit="cm",height = 15, width=15,
    res = 600)
par(mar=c(6,5,1,2))
plot(NA, xlim = c(0.75,length(primary.threats)+0.2),ylim=c(0,0.6),axes=F,xlab="",ylab="")
for(i in 1:length(primary.threats))
{
  th = primary.threats[i]
  print(th)
  y.cubrt<-do.call(rbind,lapply(q.p.th.kba.level$p.r.cube.rt[[th]], function(d) d[,4]))
  y.cubrt.lq<-do.call(rbind,lapply(q.p.th.kba.level$p.r.cube.rt[[th]], function(d) d[,3]))
  y.cubrt.uq<-do.call(rbind,lapply(q.p.th.kba.level$p.r.cube.rt[[th]], function(d) d[,5]))
  
  segments(x0 = matrix(sapply(i+unlist(x.pos), function(j) j + x.pos.adj), nrow=1),
           x1 = matrix(sapply(i+unlist(x.pos), function(j) j + x.pos.adj), nrow=1),
           y0 = matrix(t(y.cubrt.lq), nrow=1), y1 = matrix(t(y.cubrt.uq), nrow=1),
           col = sapply(cols,function(c) rep(c,4)))
  
  points(x = matrix(sapply(i+unlist(x.pos), function(j) j + x.pos.adj), nrow=1),
       y = matrix(t(y.cubrt), nrow=1),pch = 21:24, cex=0.7,
       col = sapply(cols,function(c) rep(c,4)), bg=sapply(cols,function(c) rep(c,4)))
  
  
  
  
}

abline(h = 1.0,col = "grey")
abline(v=seq(1.5,length(primary.threats)),lty="dashed")
axis(side = 1,at=seq(1,length(primary.threats)),labels = primary.threats,las=2,tick = F,cex.axis=0.8)
axis(side =2,las=1,cex.axis=0.8)
mtext(side = 2,text = "Impact probabilities in KBA impact severity classes ",cex=0.8,line = 3)
legend("top",legend = c("Mammals","Birds","Amphibians"),pch=rep(16,3), col = unlist(cols),bty='n')
legend("topright",legend = rev(c("No","Slow","Moderate","Very rapid")),pch=rev(21:24), col = "grey",pt.bg="grey",bty='n')

dev.off()



q.p.th.kba.level$p.r.cube.rt$Logging$mammal
q.p.th.kba.level$p.th$Logging$mammal



df.kba.counts<-as.data.frame(do.call("rbind",count.kbas.level))
for(c in colnames(df.kba.counts)) df.kba.counts[[c]]<-unlist(df.kba.counts[[c]])
write.csv(df.kba.counts,file = "C:/Users/mikeha/Dropbox/Threat mapping/KBA/KBACounts.csv")



sapply(primary.threats,function(th) lapply(w.2sided$p.r.cube.rt[[th]], function(x) x$p.value))




cols<-list("red","green","blue")
names(cols)<-taxa

x.pos<-list(-0.15,0,0.15)
x.pos.adj<-c(-0.045,-0.015,0.015,0.045)
names(x.pos)<-taxa

png("C:/Users/mikeha/Dropbox/Threat mapping/KBA/EvaluationRatioThreatTaxa_comparison_signif_Dissolve.png",unit="cm",height = 15, width=15,
    res = 600)
par(mar=c(6,5,1,2))
plot(NA, xlim = c(0.75,length(primary.threats)),ylim=c(0.5,10),log="y",axes=F,xlab="",ylab="")
for(i in 1:length(primary.threats))
{
  th = primary.threats[i]
  print(th)
  y.cubrt<-unlist(lapply(q.p.th.kba.level$p.r.cube.rt[[th]], function(d) d[4,4]/d[1,4]))
  y.p.th<-unlist(lapply(q.p.th.kba.level$p.th[[th]], function(d) d[4,4]/d[1,4]))
  
  points(i+unlist(x.pos),
         y.cubrt,
         pch=16,col = unlist(cols))
  # points(i+unlist(x.pos),
  #        y.p.th,
  #        pch=17,col = unlist(cols))
  
  print(unlist(lapply(q.p.th.kba.level$p.r.cube.rt[[th]], function(d) d[4,4]/d[1,4]))/
    unlist(lapply(q.p.th.kba.level$p.th[[th]], function(d) d[4,4]/d[1,4])))
  signif.p.th<-sapply(unlist(lapply(w.2sided$p.th[[th]], function(x) x$p.value)),get.signif)
  signif.cubrt<-sapply(unlist(lapply(w.2sided$p.r.cube.rt[[th]], function(x) x$p.value)),get.signif)
  # text(i+unlist(x.pos),y.p.th,signif.p.th,pos = 2,col= "grey", cex=0.8)
  text(i+unlist(x.pos),y.cubrt,signif.cubrt,pos = 2,col= "black", cex=0.8)
  #points(i,mn.th,pch=17,col = "black")
}

med.cubrt<-quantile(sapply(primary.threats,function(th) unlist(lapply(q.p.th.kba.group$p.r.cube.rt[[th]],
                                                                    function(d) (d[2,4]-d[1,4])/d[1,4]))),probs=c(0,0.025,0.25,0.5,0.75,0.975,1.0))
med.p.th<-quantile(sapply(primary.threats,function(th) unlist(lapply(q.p.th.kba.group$p.th[[th]],
                                                                      function(d) (d[2,4]-d[1,4])/d[1,4]))),probs=c(0,0.025,0.25,0.5,0.75,0.975,1.0),na.rm=T)


#med.p.th<-median(sapply(primary.threats,function(th) unlist(lapply(q.p.th.kba.level$p.th[[th]], function(d) d[4,4]/d[1,4]))))


abline(h = 1.0,col = "grey")
abline(v=seq(1.5,length(primary.threats)),lty="dashed")
axis(side = 1,at=seq(1,length(primary.threats)),labels = primary.threats,las=2,tick = F,cex.axis=0.8)
axis(side =2,las=1,cex.axis=0.8)
mtext(side = 2,text = "Ratio of median threat probability\nfor Very rapid or moderate vs no or slow severity KBAs",cex=0.8,line = 3)
legend("bottomleft",legend = c("Mammals","Birds","Amphibians"),pch=rep(16,3), col = unlist(cols),bty='n')
# legend("bottom",legend = c("Inv. cube root weight","Proportion"),pch=c(16,17), col = "grey",bty='n')
dev.off()

cubrt<-function(x) { sign(x)*(abs(x)^(1/3)) }


png("C:/Users/mikeha/Dropbox/Threat mapping/KBA/EvaluationPercentThreatTaxa_comparison_signif_Dissolve.png",unit="cm",height = 17, width=17,
    res = 600)
par(mar=c(6,5,1,3))
plot(NA, xlim = c(0.75,length(primary.threats)),ylim=cubrt(c(-0.5,6)),axes=F,xlab="",ylab="")
for(i in 1:length(primary.threats))
{
  th = primary.threats[i]
  print(th)
  #y.cubrt<-cubrt(unlist(lapply(q.p.th.kba.level$p.r.cube.rt[[th]], function(d) (d[4,4]-d[1,4])/d[1,4])))
  y.cubrt<-cubrt(unlist(lapply(q.p.th.kba.group$p.r.cube.rt[[th]], function(d) (d[2,4]-d[1,4])/d[1,4])))
  y.p.th<-cubrt(unlist(lapply(q.p.th.kba.group$p.th[[th]], function(d) (d[2,4]-d[1,4])/d[1,4])))
  
  points(i+unlist(x.pos),
         y.cubrt,
         pch=16,col = unlist(cols))
  # points(i+unlist(x.pos),
  #        y.p.th,
  #        pch=17,col = unlist(cols))

  print(unlist(lapply(q.p.th.kba.level$p.r.cube.rt[[th]], function(d) d[4,4]/d[1,4]))/
          unlist(lapply(q.p.th.kba.level$p.th[[th]], function(d) d[4,4]/d[1,4])))
  signif.p.th<-sapply(unlist(lapply(w.2sided$p.th[[th]], function(x) x$p.value)),get.signif)
  signif.cubrt<-sapply(unlist(lapply(w.2sided$p.r.cube.rt[[th]], function(x) x$p.value)),get.signif)
  # text(i+unlist(x.pos),y.p.th,signif.p.th,pos = 2,col= "grey", cex=0.8)
  text(i+unlist(x.pos),y.cubrt,signif.cubrt,pos = 2,col= "black", cex=0.8)
  #points(i,mn.th,pch=17,col = "black")
}


abline(h = 0.0,col = "grey")
abline(v=seq(1.5,length(primary.threats)),lty="dashed")
axis(side = 1,at=seq(1,length(primary.threats)),labels = primary.threats,las=2,tick = F,cex.axis=0.8)
axis(side =2,las=1,cex.axis=0.8, at = cubrt(c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,3,4,5)), labels = 100*c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,3,4,5))
mtext(side = 2,text = "Difference in median threat probability\nfor Very rapid or moderate vs no or slow severity KBAs (%)",cex=0.8,line = 3)
#arrows(x0=4.15,x1=4.15,y0 = cubrt(4),y1 = cubrt(6),code = 2,length=0.1,lwd=2,col=cols$amph)
#arrows(x0=6.15,x1=6.15,y0 = cubrt(4),y1 = cubrt(6),code = 2,length=0.1,lwd=2,col=cols$amph)
legend("topleft",legend = c("Mammals","Birds","Amphibians"),pch=rep(16,3), col = unlist(cols),bty='n',cex = 0.6)
# legend(x = 1.5,y = cubrt(6.5),legend = c("Inv. cube root weight","Proportion"),pch=c(16,17), col = "grey",bty='n',cex = 0.6)
dev.off()









png("C:/Users/mikeha/Dropbox/Threat mapping/KBA/EvaluationAbsThreatTaxa_propth.png",unit="cm",height = 15, width=17,
    res = 600)
par(mar=c(6,5,1,2))
plot(NA, xlim = c(0.75,length(primary.threats)+1),ylim=c(0,0.5),log="",axes=F,xlab="",ylab="")
for(i in 1:length(primary.threats))
{
  th = primary.threats[i]
  print(th)
  for(t in taxa){
    points(i+unlist(x.pos[[t]]) + x.pos.adj,
           unlist(q.p.th.kba.level$p.th[[th]][[t]][,4]),
           pch=c(15:18),col = cols[[t]],cex=0.75)
  }
  # signif<-sapply(unlist(lapply(w.2sided[[th]], function(x) x$p.value)),get.signif)
  # text(i+unlist(x.pos),unlist(sum.th),signif,pos = 3)
  #points(i,mn.th,pch=17,col = "black")
}
abline(v=seq(1.5,length(primary.threats)),lty="dashed")
axis(side = 1,at=seq(1,length(primary.threats)),labels = primary.threats,las=2,tick = F,cex.axis=0.8)
axis(side =2,las=1,cex.axis=0.8)
mtext(side = 2,text = "Ratio of mean threat probability in PAs with High vs Low METT threat score",cex=0.8,line = 3)
legend(x = 5.3,y = 0.5,legend = c("Mammals","Amphibians","Birds"),pch=rep(16,3), col = unlist(cols),bty='n')
legend(x = 5.3,y = 0.4,legend = rev(names(kba.severity.score.levels)),pch=rev((15:18)), col = "grey",bty='n')
dev.off()



png("C:/Users/mikeha/Dropbox/Threat mapping/METT validation/High_Low_ratios.png",unit="cm",height = 15, width=15,
    res = 600)
par(mar=c(6,5,1,2))
plot(NA, xlim = c(0.75,length(primary.threats)),ylim=c(0.2,5),log="y",axes=F,xlab="",ylab="")
for(i in 1:length(primary.threats))
{
  th = primary.threats[i]
  print(th)
  sum.th<-lapply(q.p.th.mett.level[[th]], function(d) d[2,4]/d[1,4])
  mn.th<-mean(unlist(sum.th))
  print(sum.th)
  points(i+unlist(x.pos),unlist(sum.th),pch=16,col = unlist(cols))
  signif<-sapply(unlist(lapply(w.2sided[[th]], function(x) x$p.value)),get.signif)
  text(i+unlist(x.pos),unlist(sum.th),signif,pos = 3)
  #points(i,mn.th,pch=17,col = "black")
}
abline(h = 1.0,col = "grey")
abline(v=seq(1.5,length(primary.threats)),lty="dashed")
axis(side = 1,at=seq(1,length(primary.threats)),labels = primary.threats,las=2,tick = F,cex.axis=0.8)
axis(side =2,las=1,cex.axis=0.8)
mtext(side = 2,text = "Ratio of mean threat probability in PAs with High vs Low METT threat score",cex=0.8,line = 3)
legend("bottomleft",legend = c("Mammals","Amphibians","Birds"),pch=rep(16,3), col = unlist(cols),bty='n')
dev.off()




