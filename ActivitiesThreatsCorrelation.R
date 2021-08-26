wants <- c("rgeos","rgdal","maptools","plyr","raster","RColorBrewer","betareg","lmtest")
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


#Read biome information
grid.biome.fgdb<-"C:/Users/mikeha/Documents/ArcGIS/Projects/ThreatMapIntersection/ThreatMapIntersection.gdb"
fc_list <- ogrListLayers(grid.biome.fgdb)
print(fc_list)

grid.biomes<-readOGR(dsn=grid.biome.fgdb,layer="GlobalGrid50km_Ecoregion_Biomes")

#Forested biomes
forest.biomes<-c(1,#. Tropical and subtropical moist broadleaf forests
                2,# Tropical and subtropical dry broadleaf forests
                3,# Tropical and subtropical coniferous forests
                4,# Temperate broadleaf and mixed forests
                5,# Temperate conifer forests
                6,# Boreal forests or taiga
                14)# Mangroves


forest.biome.cells<-which(grid.biomes$BIOME_NUM %in% forest.biomes)


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

PlotSingleMapAll<-function(d,cell.i,n.breaks,lower.cut,range,axis.range,fname)
{
  t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  t.labels = ''
  plot.norm = normalise_set_b(d[], t.breaks,t.labels)
  
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
#th = ThreatCodes[[phase]][4]

th = ThreatCodes$phase1[1]
th[[1]]<-th[[1]][-1] # remove residential to leave logging and agriculture

models<-c("cub.rt.R_W_Pred_","PropThSp")#,"sqrtR_W_Pred_"
model.var.names<-c("p.r.cube.rt","p.th")

path = 'C:/Users/mikeha/Dropbox/Threat mapping/'


g.q<-list()


for(i in 1:length(models))
{
  print(i)
  for(t in taxa)
  {
    print(t)
    p = c(0,0.025,0.25,0.5,0.75,0.975,1.0)
    g.q[[model.var.names[i]]][[t]]<-GetPredictionQuantiles(path = path,
                           model = models[i],model.name = model.var.names[i],t = t,th = th,
                           realms = realms,regions.comb = regions.comb,p = p)
  }
}


grid.temp<-grid.biomes



m = 1
for(m in 1:length(models))
{
  #t = taxa[2]
  for(t in taxa)
  {

    grid.temp[[paste0(m,".",t)]]<-NA
    grid.temp[[paste0(m,".",t)]][g.q[[m]][[t]]$global.i]<-(g.q[[m]][[t]]$globalq-f.loss$ProportionCellLost[g.q[[m]][[t]]$global.i])
    
    use.i<-intersect(g.q[[m]][[t]]$global.i,forest.biome.cells)
    g.i.use.i<-which(g.q[[m]][[t]]$global.i %in% forest.biome.cells)
    
    # plot(f.loss$area_loss[use.i],g.q[[t]]$globalq[g.i.use.i],pch=16,cex=0.5,col = rgb(0,0,0,0.05),
    #      xlab = "GFW Proportion of cell deforested",ylab = paste0("Predicted cell threat, ",t))
    # 
    # cor.test(f.loss$ProportionCellLost[g.q$global.i],g.q$globalq[,4])
    # abline(coef(lm(g.q$globalq[,4]~f.loss$area_forest2000[g.q$global.i])), col="red")
    
    model.data<-data.frame(predicted.threat = g.q[[m]][[t]]$globalq[g.i.use.i],
                           forest.cover = f.loss$ProportionCellCovered[use.i],
                           forest.loss = f.loss$ProportionCellLost[use.i],
                           rel.forest.loss = f.loss$RelativeLoss[use.i])
    
    
    PlotSingleMap(d = grid.temp[[paste0(m,".",t)]],
                  cell.i = g.q[[m]][[t]]$global.i[g.i.use.i],
                  n.breaks = 100,
                  lower.cut = -1,
                  range = c(-1,1),axis.range = c(-1,1),
                  fname = paste0("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/",t,"ResidualsForestLoss.png"),
                  color.scheme = brewer.pal(8,"BrBG"),
                  sqrt = F)
    
    grid.temp[[paste0(m,".",t)]][g.q[[m]][[t]]$global.i]<-(g.q[[m]][[t]]$globalq-f.loss$ProportionCellCovered[g.q[[m]][[t]]$global.i])
    
    PlotSingleMap(d = grid.temp[[paste0(m,".",t)]],
                  cell.i = g.q[[m]][[t]]$global.i[g.i.use.i],
                  n.breaks = 100,
                  lower.cut = -1,
                  range = c(-1,1),axis.range = c(-1,1),
                  fname = paste0("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/",t,"ResidualsForestCOver.png"),
                  color.scheme = brewer.pal(8,"BrBG"),
                  sqrt = F)
    
    
    # #Plot map of residuals
    # PlotSingleMap(d = (model.data$predicted.threat-model.data$forest.loss),
    #               cell.i = g.q[[m]][[t]]$global.i[g.i.use.i],
    #               n.breaks = 100,lower.cut = -1,range = c(-1,1),
    #               axis.range = c(-1,1),
    #               fname = paste0("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/",t,"ResidualsForestLoss.png"))
    # 
    # PlotSingleMap(d = (model.data$predicted.threat-model.data$forest.cover),
    #               cell.i = g.q[[m]][[t]]$global.i[g.i.use.i],
    #               n.breaks = 100,lower.cut = -1,range = c(-1,1),
    #               axis.range = c(-1,1),
    #               fname = 
    #                 paste0("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/",t,"ResidualsForestCOver.png"))
    
    
    
    
    model.data<-na.omit(model.data)
    
    #use.rows<-which(model.data$forest.cover > 0 & model.data$sr > 10) # selecting rows where there is forest and species
    # use.rows<-which(model.data$rel.forest.loss > 0.2)
    # model.data<-model.data[use.rows,]
    
    png(paste0("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/",t,"Preds-pLoss_forest_biomes_",models[m],".png"),height = 15, width = 15,
        res = 600, units = "cm")
    plot(model.data$forest.loss,model.data$predicted.threat,pch=16,cex=0.5,col = rgb(0,0,0,0.1),ylim = c(0,1.0),
         xlab = "Proportional loss of forest in 2000",ylab="Predicted threat occurrence probability")
    abline(a=0, b=1, col = "red", lty="dashed")
    
    
    #Fit a logistic regression model to test the effect of ploss on PTh
    m1<-glm(predicted.threat ~ forest.loss, data = model.data, family = binomial)
    lines(seq(0,1,0.05),
           predict(m1,newdata = data.frame(forest.loss = seq(0,1,0.05)) ,type="response"),col="blue")
    dev.off()
    
    png(paste0("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/",t,"Preds-pCover_forest_biomes_",models[m],".png"),height = 15, width = 15,
        res = 600, units = "cm")
    plot(model.data$forest.cover,model.data$predicted.threat,pch=16,cex=0.5,col = rgb(0,0,0,0.1),ylim = c(0,1.0),
         xlab = "Proportional forest cover in 2000",ylab="Predicted threat occurrence probability")
    abline(a=0, b=1, col = "red", lty="dashed")
    dev.off()
  }
}


Frmse<-function(x,y)
{
  sqrt(mean((x-y)^2, na.rm=T))
}

Fmae<-function(x,y)
{
  mean(abs(x-y), na.rm=T)
}

rmse<-list()
mae<-list()
m1<-list()
m2<-list()
m3.full<-list()
m3.interaction<-list()

taxa.labels<-list(
  mammal = "Mammals",
  bird = "Birds",
  amph = "Amphibians"
)

png(paste0("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/Preds-pLoss_forest_biomes_Cubrt.png"),height = 15, width = 15,
    res = 600, units = "cm")
par(mar = c(1,5,0.5,0.5))
layout(matrix(1:4,ncol=1), heights = c(3,3,3,1))

#Plot comparison between Cubrt model and proportional threat model
for(t in taxa)
{
  m = 1
  grid.temp[[paste0(m,".",t)]]<-NA
  grid.temp[[paste0(m,".",t)]][g.q[[m]][[t]]$global.i]<-(g.q[[m]][[t]]$globalq-f.loss$ProportionCellLost[g.q[[m]][[t]]$global.i])
  
  use.i<-intersect(g.q[[m]][[t]]$global.i,forest.biome.cells)
  g.i.use.i<-which(g.q[[m]][[t]]$global.i %in% forest.biome.cells)
  
  model.data.cbrt<-data.frame(predicted.threat = g.q[[1]][[t]]$globalq[g.i.use.i],
                         P.forest.cover = f.loss$ProportionCellCovered[use.i],
                         P.forest.loss = f.loss$ProportionCellLost[use.i],
                         forest.loss = f.loss$area_loss[use.i],
                         rel.forest.loss = f.loss$RelativeLoss[use.i])
  
  m = 2
  grid.temp[[paste0(m,".",t)]]<-NA
  grid.temp[[paste0(m,".",t)]][g.q[[m]][[t]]$global.i]<-(g.q[[m]][[t]]$globalq-f.loss$ProportionCellLost[g.q[[m]][[t]]$global.i])
  
  use.i<-intersect(g.q[[m]][[t]]$global.i,forest.biome.cells)
  g.i.use.i<-which(g.q[[m]][[t]]$global.i %in% forest.biome.cells)
  
  model.data.prop.threat<-data.frame(predicted.threat = g.q[[2]][[t]]$globalq[g.i.use.i],
                              P.forest.cover = f.loss$ProportionCellCovered[use.i],
                              P.forest.loss = f.loss$ProportionCellLost[use.i],
                              forest.loss = f.loss$area_loss[use.i],
                              rel.forest.loss = f.loss$RelativeLoss[use.i])
  
  
  model.data.cbrt<-na.omit(model.data.cbrt)
  model.data.prop.threat<-na.omit(model.data.prop.threat)
  
  rmse[[t]]<-data.frame("cbrt" = Frmse(x = model.data.cbrt$predicted.threat,model.data.cbrt$rel.forest.loss),
                        "prop.th" = Frmse(x = model.data.prop.threat$predicted.threat,model.data.prop.threat$rel.forest.loss))
  mae[[t]]<-data.frame("cbrt" = Fmae(x = model.data.cbrt$predicted.threat,model.data.cbrt$rel.forest.loss),
                       "prop.th" = Fmae(x = model.data.prop.threat$predicted.threat,model.data.prop.threat$rel.forest.loss))
  
  
  
  plot(model.data.cbrt$rel.forest.loss,model.data.cbrt$predicted.threat,pch=16,cex=0.5,col = rgb(230,85,13,0.1*255,maxColorValue = 255),ylim = c(0,1.0),xlim=c(0,1),
       xlab = "",ylab="",las=1,axes = F)
  #points(model.data.prop.threat$rel.forest.loss,model.data.prop.threat$predicted.threat,pch=16,cex=0.5,col = rgb(49,130,189,0.05*255,maxColorValue = 255))
  abline(a=0, b=1, col = "grey", lty="dashed")
  

  axis(side = 2, las = 1)
  if(t == "bird") mtext(side = 2,text = "Predicted threat occurrence probability",line = 3.5,cex = 0.8)
  mtext(side = 2, text= taxa.labels[[t]], line = 2.5, cex = 0.8)
  
  if(t == "amph")
  {
    axis(side = 1)
    mtext(side = 1,text = "Proportional loss of forest in 2000",line = 2,cex = 0.8)
  } else {
    axis(side = 1,labels = F)
  }
  
  #Fit a logistic regression model to test the effect of ploss on PTh
  m1[[t]]<-glm(predicted.threat ~ rel.forest.loss, data = model.data.cbrt, family = binomial)
  lines(seq(0,1,0.05),
        predict(m1[[t]],newdata = data.frame(rel.forest.loss = seq(0,1,0.05)) ,type="response"),col="black",lwd=1.5)#rgb(230,85,13,maxColorValue = 255)
  m2[[t]]<-glm(predicted.threat ~ rel.forest.loss, data = model.data.prop.threat, family = binomial)
  # lines(seq(0,1,0.05),
  #       predict(m2[[t]],newdata = data.frame(rel.forest.loss = seq(0,1,0.05)) ,type="response"),col = rgb(8,81,156,maxColorValue = 255),lwd=1.5)
  # 
  model.data.cbrt$Model<-rep("Cube root",nrow(model.data.cbrt))
  model.data.prop.threat$Model<-rep("Prop threatened",nrow(model.data.prop.threat))
  model.data.all<-rbind(model.data.cbrt,model.data.prop.threat)
  
  
  m3.full[[t]]<-glm(predicted.threat ~ rel.forest.loss, data = model.data.all, family = binomial)
  m3.interaction[[t]]<-glm(predicted.threat ~ rel.forest.loss * Model, data = model.data.all, family = binomial)

}
dev.off()


lapply(m1, summary)




plot(model.data$forest.cover,model.data$forest.loss,pch=16,cex=0.5,col = rgb(0,0,0,0.1))


cor.test(model.data$forest.loss,model.data$predicted.threat)

plot(model.data$rel.forest.loss,model.data$predicted.threat,pch=16,cex=0.5,col = rgb(0,0,0,0.1))
cor.test(model.data$rel.forest.loss,model.data$predicted.threat)
plot(model.data$forest.cover,model.data$predicted.threat,pch=16,cex=0.5,col = rgb(0,0,0,0.1))


lms = list()
lms[[names(th)]][[model]][[t]]<-glm(sqrt(predicted.threat) ~ forest.cover*forest.loss,
                                    weights =  model.data$sr, data = model.data,family = binomial(link = "logit"))

summary(lms[[names(th)]][[model]][[t]])


plot(model.data$predicted.threat,(predict(lms[[names(th)]][[model]][[t]], type = "response"))^2,
     pch=16,cex =0.5, col = rgb(0,0,0,0.1), xlim = c(0,1), ylim = c(0,1))
abline(a=0,b = 1,col="grey")
cor.test(model.data$predicted.threat,(predict(lms[[names(th)]][[model]][[t]], type = "response")))

plot(model.data$forest.loss,residuals(lms[[names(th)]][[model]][[t]]),pch=16,cex=0.5,
     col = rgb(0,0,0,0.1))



betaregs<-list()
betaregs[[names(th)]][[model]][[t]]<-betareg(predicted.threat ~ forest.cover+forest.loss, 
                                             data = model.data,na.action = "na.exclude")
plot(model.data$predicted.threat,(predict(betaregs[[names(th)]][[model]][[t]], type = "response")),
     pch=16,cex =0.5, col = rgb(1,0,0,0.1), xlim = c(0,1), ylim = c(0,1))
abline(a=0,b = 1,col="grey")
cor.test(model.data$predicted.threat,(predict(betaregs[[names(th)]][[model]][[t]], type = "response")))
