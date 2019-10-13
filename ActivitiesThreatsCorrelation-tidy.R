wants <- c("rgeos","rgdal","maptools","plyr","raster","RColorBrewer","betareg")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
has <- wants %in% .packages()
for(p in wants[!has]) library(p,character.only = T)

out.dir = 'C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/'
path = 'C:/Users/mikeha/Dropbox/Threat mapping/'
source(paste0(path,'scripts/ThreatMapFunctions.R'))
taxa = c('mammal','bird','amph')

wgs1984.proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
moll.proj <- CRS("+proj=moll +datum=WGS84")
grid.wgs84<-readOGR("C:/Users/mikeha/Dropbox/Threat mapping/Regional","GlobalGrid50kmRoundl_Clip_Ecoregion_wgs84")
grid<-spTransform(grid.wgs84,moll.proj)

r.grid<-raster(extent(grid),nrows = (bbox(grid)[2,2]-bbox(grid)[2,1])/500,
               ncols = (bbox(grid)[1,2]-bbox(grid)[1,1])/500)




#Read in map of the world
WL <- readShapePoly("C:/Users/mikeha/Dropbox/Threat mapping/Spatial data/ne_10m_land.shp")#download from 'naturalearth.com'
proj4string(WL) <- wgs1984.proj  
WL.moll <- spTransform(WL, moll.proj)




# Read in the Hansen deforestation layer
f.loss<-read.csv("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/Deforestation/cover2000_and_loss_per_cell.csv",stringsAsFactors = F)
f.loss<-f.loss[,-which(colnames(f.loss) == ".geo")]
f.loss<-f.loss[order(f.loss$Id),]
gc()
f.loss$RelativeLoss<-f.loss$area_loss/f.loss$area_forest2000

cell.areas<-gArea(grid,byid = T)
f.loss$ProportionCellLost<-f.loss$area_loss/cell.areas
#f.loss$ProportionCellLost[f.loss$area_forest2000 == 0]<-NA
f.loss$ProportionCellCovered<-f.loss$area_forest2000/cell.areas
#f.loss$ProportionCellCovered[f.loss$area_forest2000 == 0]<-NA
#plot.order<-sapply(f.loss$Id, which(grid@data$Id == x))


values(r.grid)<-f.loss$area_loss

#image(as.matrix(f.loss))

PlotSingleMap<-function(d,cell.i,n.breaks,lower.cut,range,axis.range,fname)
{
  t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  t.labels = ''
  plot.norm = normalise_set_b(d, t.breaks,t.labels)
  
  features.to.plot = which(plot.norm[[1]] > 0)
  
  png(fname,width = 17,height = 12, units = 'cm',res = 600)
  
  par(oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  par(xpd=NA)
  par(cex=1)
  
  layout(matrix(c(1,2),nrow = 2), heights = c(0.9,0.1))
  
  plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
  #plot(bound,lwd=0.5,col = NA,border = "black", add=T)
  
  cols <- (colorRampPalette((brewer.pal(8,"Reds")))(n.breaks))
  palette(cols)
  plot(grid[features.to.plot,],col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
  
  
  par(mar=c(2,3,0,3))
  
  x.lim = axis.range
  y.lim = c(0,1)
  
  width = (axis.range[2]-axis.range[1])/n.breaks
  x_step = seq(x.lim[1], x.lim[2], width)
  
  plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
  
  #rect(0, 0, 1, 1, col="light grey", border=NA)
  
  sapply(seq(1,n.breaks),
         function(i) {
           rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
           #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
         })
  #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
  
  x.values = seq(x.lim[1],x.lim[2],length.out = 6)
  axis(1,padj=-2, cex.axis=0.5,
       at = x.values)
  
  
  dev.off()
}

PlotSingleMap(d = (f.loss$ProportionCellCovered),
              cell.i =seq_len(dim(grid)[1]),
              n.breaks = 100,lower.cut = 0,range = c(0,1),
              axis.range = c(0,1),
              fname = paste0("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/Deforestation/PcellForestCover.png"))



GetPredictionQuantiles<-function(path,model,t,th,realms, regions.comb,p)
{
  if(exists("global.q")) rm(global.q)
  global.i = NULL
  global.sr = NULL
  #Loop over realms
  for(r in realms)
  {
    print(r)
    #load predictions
    load(file = paste0(path,'Occupancy modelling/Logistic/',r,'/',model,t,'_',names(th)))
    #if(model == "logR_W_Pred_")
    #{
    #  p.r = p.r.cube.rt
    #}
    # else if(model == "sqrtR_W_Pred_")
    # {
    #   p.r = p.r.sqrt
    # }
    
    #q.p.r = apply(p.r,1,function(x) quantile(x,probs = p, na.rm=T))
    q.p.r<-p.r.cube.rt#apply(p.r,1,function(x) max(x,na.rm=T))
    q.p.r[is.infinite(q.p.r)]<-NA
    #global.sr = c(global.sr,apply(p.r,1,function(x) sum(!is.na(x))))
    
    
    load(paste0(path,"Regional 2017_3/",r,"/",t,"_intersection.a.sp.R"))
    #load(paste0(path,"Regional 2017_3/Rdata/",r,"/",t,"_ranges.R"))
    
    if(t == 'mammal')
    {
      inter.a = r.mammal.a
    } else if(t == 'bird')
    {
      inter.a = r.birds.a
    } else if(t == 'amph')
    {
      inter.a = r.amphibian.a
    }
    
    
    
    global.sr = c(global.sr,apply(inter.a,1, function(x) sum(x)))
    
    #add these data to the global set
    if(exists("global.q"))
    {
      #global.q = rbind(global.q,(q.p.r))
      global.q<-c(global.q,q.p.r)
    } else
    {
      #global.q = t(q.p.r)
      global.q<-q.p.r
    }
    
    
    
    #find the cell indices of these quantiles
    r.i = which(grid@data$REALM == r) 
    global.i = c(global.i,r.i)
  }
  
  for(r in 1:length(regions.comb))
  {
    region = names(regions.comb)[r]
    print(region)
    #load predictions
    load(file = paste0(path,'Occupancy modelling/Logistic/',region,'/',model,t,'_',names(th)))
    # 
    # if(model == "logR_W_Pred_")
    # {
    #   p.r = p.r.log
    # }
    # else if(model == "sqrtR_W_Pred_")
    # {
    #   p.r = p.r.sqrt
    # }
    
    q.p.r<-p.r.cube.rt#apply(p.r,1,function(x) max(x,na.rm=T))
    q.p.r[is.infinite(q.p.r)]<-NA
    #global.sr = c(global.sr,apply(p.r,1,function(x) sum(!is.na(x))))
    
    
    load(paste0(path,"Regional 2017_3/",region,"/",t,"_intersection.a.sp.R"))
    #load(paste0(path,"Regional 2017_3/Rdata/",r,"/",t,"_ranges.R"))
    
    if(t == 'mammal')
    {
      inter.a = r.mammal.a
    } else if(t == 'bird')
    {
      inter.a = r.birds.a
    } else if(t == 'amph')
    {
      inter.a = r.amphibian.a
    }
    
    
    
    global.sr = c(global.sr,apply(inter.a,1, function(x) sum(x)))
    
    #global.sr = c(global.sr,apply(p.r,1,function(x) sum(!is.na(x))))
    
    #add these data to the global set
    #global.q = rbind(global.q,t(q.p.r))
    global.q<-c(global.q,q.p.r)
    
    r.i = which(grid@data$REALM == regions.comb[[r]][1] | grid@data$REALM == regions.comb[[r]][2])
    global.i = c(global.i,r.i)
  }
  
  
  list(
    global.i = global.i,
    global.sr = global.sr,
    globalq = global.q
    )

}


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

taxa = c('mammal','bird','amph','All')

ThreatCodes = GetThreatCodesPhases()

phase = 'phase1'
th = ThreatCodes[[phase]][4]
models<-c("cub.rt.R_W_Pred_")#,"sqrtR_W_Pred_"

path = 'C:/Users/mikeha/Dropbox/Threat mapping/'
model = models[1]
t = taxa[1]


p = c(0,0.025,0.25,0.5,0.75,0.975,1.0)
g.q<-GetPredictionQuantiles(path = path,
                       model = model,t = t,th = th,
                       realms = realms,regions.comb = regions.comb,p = p)

plot(f.loss$area_loss[g.q$global.i],g.q$globalq,pch=16,cex=0.5,col = rgb(0,0,0,0.05),
     xlab = "GFW Proportion of cell deforested",ylab = paste0("Predicted cell threat, ",t))
cor.test(f.loss$ProportionCellLost[g.q$global.i],g.q$globalq[,4])
abline(coef(lm(g.q$globalq[,4]~f.loss$area_forest2000[g.q$global.i])), col="red")

model.data<-data.frame(predicted.threat = g.q$globalq,
                       forest.cover = f.loss$ProportionCellCovered[g.q$global.i],
                       forest.loss = f.loss$ProportionCellLost[g.q$global.i],
                       rel.forest.loss = f.loss$RelativeLoss[g.q$global.i],
                       sr = g.q$global.sr)


model.data<-na.omit(model.data)

#use.rows<-which(model.data$forest.cover > 0 & model.data$sr > 10) # selecting rows where there is forest and species
# use.rows<-which(model.data$rel.forest.loss > 0.2)
# model.data<-model.data[use.rows,]

png("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/MammalPreds-pLoss.png",height = 15, width = 15,
    res = 600, units = "cm")
plot(model.data$forest.loss,model.data$predicted.threat,pch=16,cex=0.5,col = rgb(0,0,0,0.1),
     xlab = "Proportional loss of forest in 2000",ylab="Predicted threat occurrence probability")
abline(a=0, b=1, col = "red", lty="dashed")
dev.off()

#Fit a logistic regression model to test the effect of ploss on PTh
m1<-glm(predicted.threat ~ forest.loss, data = model.data, family = binomial)



png("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/MammalPreds-pCover.png",height = 15, width = 15,
    res = 600, units = "cm")
plot(model.data$forest.cover,model.data$predicted.threat,pch=16,cex=0.5,col = rgb(0,0,0,0.1),
     xlab = "Proportional forest cover in 2000",ylab="Predicted threat occurrence probability")
abline(a=0, b=1, col = "red", lty="dashed")
dev.off()

plot(model.data$forest.cover,model.data$forest.loss,pch=16,cex=0.5,col = rgb(0,0,0,0.1))


cor.test(model.data$forest.loss,model.data$predicted.threat)

plot(model.data$rel.forest.loss,model.data$predicted.threat,pch=16,cex=0.5,col = rgb(0,0,0,0.1))
cor.test(model.data$rel.forest.loss,model.data$predicted.threat)
plot(model.data$forest.cover,model.data$predicted.threat,pch=16,cex=0.5,col = rgb(0,0,0,0.1))


lms = list()
lms[[names(th)]][[model]][[t]]<-glm(sqrt(predicted.threat) ~ forest.cover*forest.loss,
                                    weights =  model.data$sr, data = model.data,family = binomial(link = "logit"))

summary(lms[[names(th)]][[model]][[t]])


