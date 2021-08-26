wants <- c("rgeos","rgdal","maptools","plyr","raster","RColorBrewer","viridis","ash","spdep")
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


normalise_set_b <- function(v,b, labels)
{
  #na.v.i = which(v == -1.0)
  #v = log(v)
  v[is.infinite(v)] = NA
  
  v.norm = as.integer(cut(v,breaks=b,include.lowest = T))
  
  v.norm[is.na(v.norm)] = -1.0
  #v.norm = v.norm-1.0
  
  return(list(v.norm,labels))
}

#Read in map of the world
WL <- readShapePoly("C:/Users/mikeha/Dropbox/Threat mapping/Spatial data/ne_10m_land.shp")#download from 'naturalearth.com'
proj4string(WL) <- wgs1984.proj  
WL.moll <- spTransform(WL, moll.proj)

WL.moll<-crop(WL.moll,extent(c(-17702327,17880979,-7.0E6,8.6E6)))
WL.moll.thin<-crop(WL.moll, extent(-1.5E7,1.75E7,-7E6,8.6E6))
WL.moll<-WL.moll.thin


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

taxa = c('mammal','amph','bird')

ThreatCodes = GetThreatCodesPhases()

phase = 'phase1'

p = c(0,0.025,0.25,0.5,0.75,0.975,1.0)

grid_DDThreatCodes_Uncert<-array(NA, dim = c(nrow(grid@data),length(taxa),
                                             length(ThreatCodes[[phase]]),length(p)))
grid_DDCells_Uncert<-array(NA, dim = c(nrow(grid@data),length(taxa),
                                       length(ThreatCodes[[phase]]),length(p)))


#th = ThreatCodes[[phase]][1]

models<-c("cub.rt.R_W_Pred_")#,"sqrtR_W_Pred_"

# Extract predictions for mammals and amphibians ######
#for(model in models)
model<-models[1]
{
  
  for(t in taxa)
  {
    print(t)
    
    for(th.i in 1:length(ThreatCodes[[phase]]))
    {
      th <- ThreatCodes[[phase]][th.i]
    
      print(names(th))
      
      if(exists("global.q")) rm(global.q)
      global.i = NULL
      global.sr = NULL
      
      DD_ThreatCodes_preds<-NULL
      DD_CellCodes_preds<-NULL
      
      #Loop over realms
      for(r in realms)
      {
        print(r)
        q.p.r<-NULL
        
        
        load(paste0(path,"Regional 2017_3/",r,"/",t,"_intersection.a.sp.R"))
        
        load(file = paste0(path,'Occupancy modelling/Logistic/',r,'/',model,t,'_',names(th)))
        q.p.r<-p.r.cube.rt
        q.p.r[is.infinite(q.p.r)]<-NA
        if(t == 'mammal')
        {
          inter.a = r.mammal.a
        } else if(t == 'bird')
        {
          inter.a = r.bird.a
        } else if(t == 'amph')
        {
          inter.a = r.amphibian.a
        }
      
        
        
        global.sr = c(global.sr,apply(inter.a,1, function(x) sum(x, na.rm=T)))
        
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
        r.i<-if(t == "bird") as.integer(rownames(r.bird.a)) else (which(grid@data$REALM == r)) 
        
        #test where the indices and predictions misalign
        if(length(r.i) != length(q.p.r)) print("ERROR - misaligned")
        
        global.i = c(global.i,r.i)
        
        #Now read in the uncertainty matrices: into object called predictions
        load(file = paste0(path,'Occupancy modelling/Logistic/',r,'/cub.rt.R_W_Pred_Matrix_DDThreatCodes',t,'_',names(th)))
        DD_ThreatCodes_preds<-rbind(DD_ThreatCodes_preds,t(apply(predictions,1,quantile,probs=p, na.rm=T)))
        load(file = paste0(path,'Occupancy modelling/Logistic/',r,'/cub.rt.R_W_Pred_Matrix_DD',t,'_',names(th)))
        DD_CellCodes_preds<-rbind(DD_CellCodes_preds,t(apply(predictions,1,quantile,probs=p, na.rm=T)))
        
      }
      
      for(r in 1:length(regions.comb))
      {
        region = names(regions.comb)[r]
        print(region)
        
        q.p.r<-NULL
        
        
        load(paste0(path,"Regional 2017_3/",region,"/",t,"_intersection.a.sp.R"))
        
        if(t == 'mammal')
        {
          load(file = paste0(path,'Occupancy modelling/Logistic/',region,'/',model,t,'_',names(th)))
          q.p.r<-p.r.cube.rt#apply(p.r,1,function(x) max(x,na.rm=T))
          q.p.r[is.infinite(q.p.r)]<-NA
          inter.a = r.mammal.a
        } else if(t == 'bird')
        {
          load(file = paste0(path,'Occupancy modelling/Logistic/',region,'/',model,t,'_',names(th)))
          q.p.r<-p.r.cube.rt
          q.p.r[is.infinite(q.p.r)]<-NA
          inter.a = r.bird.a
        } else if(t == 'amph')
        {
          load(file = paste0(path,'Occupancy modelling/Logistic/',region,'/',model,t,'_',names(th)))
          q.p.r<-p.r.cube.rt#apply(p.r,1,function(x) max(x,na.rm=T))
          q.p.r[is.infinite(q.p.r)]<-NA
          inter.a = r.amphibian.a
        }
        
        global.sr = c(global.sr,apply(inter.a,1, function(x) sum(x, na.rm=T)))
        
        
        #add these data to the global set
        #global.q = rbind(global.q,t(q.p.r))
        global.q<-c(global.q,q.p.r)
        
        r.i<-if(t == "bird") as.integer(rownames(r.bird.a)) else (which(grid@data$REALM == regions.comb[[r]][1] | grid@data$REALM == regions.comb[[r]][2]))
        global.i = c(global.i,r.i)
        
        #Now read in the uncertainty matrices: into object called predictions
        load(file = paste0(path,'Occupancy modelling/Logistic/',region,'/cub.rt.R_W_Pred_Matrix_DDThreatCodes',t,'_',names(th)))
        DD_ThreatCodes_preds<-rbind(DD_ThreatCodes_preds,t(apply(predictions,1,quantile,probs=p, na.rm=T)))
        load(file = paste0(path,'Occupancy modelling/Logistic/',region,'/cub.rt.R_W_Pred_Matrix_DD',t,'_',names(th)))
        DD_CellCodes_preds<-rbind(DD_CellCodes_preds,t(apply(predictions,1,quantile,probs=p, na.rm=T)))
        
      }
      
      
      use.cells<-which(global.sr > 10)
      
      
      grid@data[[paste0(model,t,'_',names(th))]]<-NA
      grid@data[[paste0(model,t,'_',names(th))]][global.i]<-global.q
      
      
      grid@data[[paste0(model,t,'_',names(th),'gt10sp')]]<-NA
      grid@data[[paste0(model,t,'_',names(th),'gt10sp')]][global.i[use.cells]]<-global.q[use.cells]
      
      if(th.i == 1)
      {
        grid@data[[paste0(t,'_',"SR")]]<-NA
        grid@data[[paste0(t,'_',"SR")]][global.i]<-unlist(global.sr)
        
        grid@data[[paste0(t,'_',"SRgt10sp")]]<-NA
        grid@data[[paste0(t,'_',"SRgt10sp")]][global.i[use.cells]]<-global.sr[use.cells]
        
      }
      
      grid_DDThreatCodes_Uncert[global.i,which(taxa == t),th.i,] <- DD_ThreatCodes_preds
      grid_DDCells_Uncert[global.i,which(taxa == t),th.i,] <- DD_CellCodes_preds
      
      # for(j in seq(ncol(DD_ThreatCodes_preds))){
      #   grid_DDThreatCodes_Uncert[[paste0(model,t,'_',names(th),j)]]<-NA
      #   grid_DDThreatCodes_Uncert[[paste0(model,t,'_',names(th),j)]][global.i]<-DD_ThreatCodes_preds[,j]
      # }
      
      # for(j in seq(ncol(DD_CellCodes_preds))){      
      #   grid_DDCells_Uncert[[paste0(model,t,'_',names(th),j)]]<-NA
      #   grid_DDCells_Uncert[[paste0(model,t,'_',names(th),j)]][global.i]<-DD_CellCodes_preds[,j]
      # }
      
    
    }
  }
}


colnames(grid@data)[3:ncol(grid@data)] = sub(pattern = models[1],
                                             replacement = "cubrt_",
                                             x = colnames(grid@data)[3:ncol(grid@data)])

#Write out the global probability of impact shapefile
writeOGR(obj = grid,dsn = "C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic",
         layer = "Grid_mammal_amph_threat_predictions",
         driver = "ESRI Shapefile",overwrite_layer = T)


#Restart point here - read in preconstructed grid
grid<-readOGR("C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic","Grid_mammal_amph_threat_predictions")


# plotting threats using threat specific colours -----------------------------


color.schemes<-list(
  "Logging" = brewer.pal(8,"Greens"),
  "Agriculture" = rev(brewer.pal(11,"BrBG")[1:6]),
  "Hunting" = brewer.pal(8,"Reds"),
  "Pollution" = brewer.pal(8,"Blues"),
  "Invasives" = brewer.pal(8,"Purples"),
  "Climate change" = brewer.pal(8,"Oranges")
)


primary.threats<-names(color.schemes) 

for(t in taxa)
{
  print(t)
  output.dir<-paste0(out.dir,t,"/")
  if(!dir.exists(output.dir)) dir.create(output.dir)
  
  sr.col.i<-grep(paste0(t,"_SR"),colnames(grid@data))[1]
  global.i<-which(!is.na(grid@data[[sr.col.i]]))
  # PlotSingleMap(d = grid@data[[grep(paste0(t,"_SR"),colnames(grid@data))[1]]],
  #               cell.i = global.i,
  #               n.breaks = 100,
  #               lower.cut = 0,
  #               range = range(grid@data[[grep(paste0(t,"_SR"),colnames(grid@data))[1]]],na.rm=T),
  #               axis.range = range(grid@data[[grep(paste0(t,"_SR"),colnames(grid@data))[1]]],na.rm=T),
  #               fname = paste0(output.dir,t,"_SR.png"),color.scheme = "Reds")
  # 
  
  gt10.col.i<-grep(paste0(t,"_SRgt10sp"),colnames(grid@data))
  gt10sp.i<-which(!is.na(grid@data[[gt10.col.i]]))
  # PlotSingleMap(d = grid@data[[grep(paste0(t,"_SRgt10sp"),colnames(grid@data))]],
  #               cell.i = gt10sp.i,
  #               n.breaks = 100,
  #               lower.cut = 0,
  #               range = range(grid@data[[grep(paste0(t,"_SRgt10sp"),colnames(grid@data))]],na.rm=T),
  #               axis.range = range(grid@data[[grep(paste0(t,"_SRgt10sp"),colnames(grid@data))]],na.rm=T),
  #               fname = paste0(output.dir,t,"_SR.png"),color.scheme = "Reds")
  
  for(th in primary.threats)
  {
    print(th)
    PlotSingleMap(d = grid@data[[grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data))]],
                  cell.i = gt10sp.i,
                  n.breaks = 100,
                  lower.cut = 0,
                  range = c(0,1),
                  axis.range = c(0,1),
                  fname = paste0(output.dir,t,"_",th,"gt10sp.png"),
                  color.scheme = color.schemes[[th]])
    
    PlotSingleMap(d = grid@data[[grep(paste0(t,"_",th),colnames(grid@data))[1]]],
                  cell.i = global.i,
                  n.breaks = 100,
                  lower.cut = 0,
                  range = c(0,1),
                  axis.range = c(0,1),
                  fname = paste0(output.dir,t,"_",th,".png"),
                  color.scheme = color.schemes[[th]])
    
  }
  
  #Plot the probability of threat from any threat
  # th = "AnyThreat"
  # print(th)
  # PlotSingleMap(d = grid@data[[grep(paste0(t,"_",th,"gt10sp"),colnames(grid@data))]],
  #               cell.i = gt10sp.i,
  #               n.breaks = 100,
  #               lower.cut = 0,
  #               range = c(0,1),
  #               axis.range = c(0,1),
  #               fname = paste0(output.dir,t,"_",th,"gt10sp.png"),
  #               color.scheme = color.schemes$Hunting)
  # 
  # PlotSingleMap(d = grid@data[[grep(paste0(t,"_",th),colnames(grid@data))[1]]],
  #               cell.i = global.i,
  #               n.breaks = 100,
  #               lower.cut = 0,
  #               range = c(0,1),
  #               axis.range = c(0,1),
  #               fname = paste0(output.dir,t,"_",th,".png"),
  #               color.scheme = color.schemes$Hunting)
  
}




#Write out shapefiles for rasterization and zonal stats with METT PAs.
use.cols<-(grep("gt10sp",colnames(grid@data), value = T))
use.cols<-use.cols[-grep("SR",use.cols)]

for(c in use.cols)
  writeOGR(obj = grid[,c],
           dsn = "C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic",
           layer = c,
           driver = "ESRI Shapefile",overwrite_layer = T)



#Construct ecoregion summary statistics for each taxa


x.at<-as.vector(sapply(1:7, function(x) seq(6)+((x-1)*7)))
xlim = c(0,max(x.at))
for(t in taxa)
{
  output.dir<-paste0(out.dir,t,"/")
  ecor.t.th.quants = lapply(primary.threats, function(th) aggregate(grid@data[grid@data$REALM != "Antarctica",c(grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data),value = T))],
              by = list(realm = grid@data$REALM[grid@data$REALM != "Antarctica"]),FUN = function(x) quantile(x,probs=c(0.025,0.5,0.975), na.rm=T)))
  
  df.ecor.t.th.quants<-do.call("cbind",lapply(ecor.t.th.quants,function(x) x[,2]))
  colnames(df.ecor.t.th.quants)<-paste(sapply(primary.threats,function(s) rep(s,3)),colnames(ecor.t.th.quants[[1]][,2]))
  df.ecor.t.th.quants<-cbind(ecor.t.th.quants[[1]][1],df.ecor.t.th.quants)
  
  # cols<-sapply(primary.threats,function(th)
  #   if(color.schemes[[th]] == "Browns") browns[4] else (brewer.pal(3,color.schemes[[th]])[3]))
  #   
  # 
  # png(paste0(output.dir,t,"_","ecoregion_P(Th)_summaries.png"),width = 17,height = 12, units = 'cm',res = 600)
  # 
  # plot(x.at,as.vector(sapply(seq(7),function(x) as.numeric(df.ecor.t.th.quants[x,grep("50",colnames(df.ecor.t.th.quants))]))),
  #      pch=16,ylim =c(0,1),axes=F,xlab="",ylab="",col=cols)
  # segments(x0 = x.at,x1=x.at,
  #   y0=as.vector(sapply(seq(7),function(x) as.numeric(df.ecor.t.th.quants[x,grep("2.5",colnames(df.ecor.t.th.quants))]))),
  #   y1=as.vector(sapply(seq(7),function(x) as.numeric(df.ecor.t.th.quants[x,grep("97.5",colnames(df.ecor.t.th.quants))]))),col=cols)
  # abline(v = seq(7)*7, col = "dark grey")
  # axis(side = 1,at = (seq(7)*7)-3.5,labels = df.ecor.t.th.quants[,1],las = 2,tick = F,line = -1)
  # axis(side = 2,las=1)
  # dev.off()
}



# Plotting uncertainties -------------


for(t in taxa)
{
  print(c(min(grid_DDThreatCodes_Uncert[,t.i,,6]-grid_DDThreatCodes_Uncert[,t.i,,2], na.rm=T),
             (max(grid_DDThreatCodes_Uncert[,t.i,,6]-grid_DDThreatCodes_Uncert[,t.i,,2], na.rm=T))))
}

#plot Multipanel figure of by threat taxa uncertainty maps
for(t in taxa)
{
  print(t)
  output.dir<-paste0(out.dir,t,"/")
  t.i<-which(taxa == t)
  
  t.range<-c(min(grid_DDThreatCodes_Uncert[,t.i,,6]-grid_DDThreatCodes_Uncert[,t.i,,2], na.rm=T),
             ceiling(max(grid_DDThreatCodes_Uncert[,t.i,,6]-grid_DDThreatCodes_Uncert[,t.i,,2], na.rm=T)))
  
  t.range=sqrt(c(0,0.7))
  
  #png(paste0(paste0(output.dir,"MultipanelUncertaintyThreatCodes_SRgt10sp"),".png"),height = 20, width = 20,units = "cm",res = 600)
  pdf(paste0(paste0(output.dir,"MultipanelUncertaintyThreatCodes_SRgt10sp"),".pdf"),height = 21/2.54, width = 18/2.54)
  layout(matrix(1:12,ncol = 2,byrow = F),heights = rep(c(0.85,0.15)/3,3))
  #layout(matrix(1:6,ncol = 2,byrow = F))
  i = 1
  for(th in primary.threats)
  {
    print(th)
    th.i = which(names(ThreatCodes[[phase]]) == th)
    
    gt10.col.i<-grep(paste0(t,"_SRgt10sp"),colnames(grid@data))
    gt10sp.i<-which(!is.na(grid@data[[gt10.col.i]]))
    
    
    PlotSingleMapInMultiPanel(d = sqrt(grid_DDThreatCodes_Uncert[,t.i,th.i,6]-grid_DDThreatCodes_Uncert[,t.i,th.i,2]),
                              cell.i = gt10sp.i,
                              n.breaks = 100,
                              lower.cut = 0,
                              range = t.range,
                              axis.range = t.range,
                              color.scheme = color.schemes[[th]],i = i,sqrt = T)
    i = i+1
  }
  dev.off()
  
  # t.range<-c(min(grid_DDCells_Uncert[,t.i,,6]-grid_DDCells_Uncert[,t.i,,2], na.rm=T),
  #            ceiling(max(grid_DDCells_Uncert[,t.i,,6]-grid_DDCells_Uncert[,t.i,,2], na.rm=T)))
  # t.range=sqrt(c(0,0.7))
  # 
  # png(paste0(paste0(output.dir,"MultipanelUncertaintyCellCodes_SRgt10sp"),".png"),height = 20, width = 20,units = "cm",res = 600)
  # layout(matrix(1:12,ncol = 2,byrow = F),heights = rep(c(0.85,0.15)/3,3))
  # #layout(matrix(1:6,ncol = 2,byrow = F))
  # i = 1
  # for(th in primary.threats)
  # {
  #   print(th)
  #   th.i = which(names(ThreatCodes[[phase]]) == th)
  #   
  #   gt10.col.i<-grep(paste0(t,"_SRgt10sp"),colnames(grid@data))
  #   gt10sp.i<-which(!is.na(grid@data[[gt10.col.i]]))
  #   
  #   
  #   PlotSingleMapInMultiPanel(d = sqrt(grid_DDCells_Uncert[,t.i,th.i,6]-grid_DDCells_Uncert[,t.i,th.i,2]),
  #                             cell.i = gt10sp.i,
  #                             n.breaks = 100,
  #                             lower.cut = 0,
  #                             range = t.range,
  #                             axis.range = t.range,
  #                             color.scheme = color.schemes[[th]],i = i,sqrt = T)
  #   i = i+1
  # }
  # dev.off()
  # 
  # png(paste0(paste0(output.dir,"MultipanelPThvsUncertainty"),".png"),height = 20, width = 20,units = "cm",res = 600)
  # layout(matrix(1:6, ncol=2))
  # for(th in primary.threats)
  # {
  #   print(th)
  #   th.i = which(names(ThreatCodes[[phase]]) == th)
  #   
  #   gt10.col.i<-grep(paste0(t,"_SRgt10sp"),colnames(grid@data))
  #   gt10sp.i<-which(!is.na(grid@data[[gt10.col.i]]))
  #   
  #   colname<-grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data), value = T)
  #   
  #   
  #   plot(grid@data[[colname]][gt10sp.i],
  #        grid_DDThreatCodes_Uncert[gt10sp.i,t.i,th.i,6]-grid_DDThreatCodes_Uncert[gt10sp.i,t.i,th.i,2],
  #        pch=16,cex=0.5,col=rgb(0,0,0,0.1), las = 1,ylab="95% range",xlab=paste0("PTh(",th,")"))
  # 
  #   
  #   
  # }
  # dev.off() 
  #   
  
}







# plotting composite maps -----------------------------

max.threat<-list()
second.threat<-list()
cumulative.threat<-list()
mean.threat<-list()
median.threat<-list()
sum.threat<-list()
for(t in taxa)
{
  
  print(t)
  gt10.col.i<-grep(paste0(t,"_SRgt10sp$"),colnames(grid@data))
  gt10sp.i<-which(!is.na(grid@data[[gt10.col.i]]))
  
  primary.threat.columns<-unlist(lapply(primary.threats,function(th) grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data))))
  
  #use.colours<-unlist(lapply(primary.threats,function(th) brewer.pal(3,color.schemes[[th]])[3]))
  
  # max.threat[[t]]<-unlist(apply(grid@data[,primary.threat.columns],1,function(x)
  #   ifelse(sum(is.na(x)) == length(x),yes = NA, no = which(!is.na(x))[which.max(x[!is.na(x)])])))
  # 
  # second.threat[[t]]<-unlist(apply(grid@data[,primary.threat.columns],1,function(x)
  # {
  #   i.no.na <- which(!is.na(x))
  #   x.no.na<-x[i.no.na]
  #   i.2nd<-NA
  #   if(length(x.no.na) > 1)
  #   {
  #     i.max<-which.max(x.no.na)
  #     x.no.na.no.max<-x.no.na[-i.max]
  #     i.2nd<-which.max(x.no.na.no.max[(x.no.na.no.max > (0.9*x.no.na[i.max]))])  
  #   }
  #   i.2nd
  # }))
  # 
  
  order.threats<-t(apply(grid@data[,primary.threat.columns],1,order, decreasing = T))
  sum.not.na<-apply(grid@data[,primary.threat.columns],1, function(x) sum(!is.na(x)))
  
  # use.second.threat<-which(grid@data[,primary.threat.columns[order.threats[,2]]] > 
  #           0.99*grid@data[,primary.threat.columns[order.threats[,1]]])
  
  max.threat[[t]]<-rep(NA, nrow(grid@data))
  max.threat[[t]][sum.not.na > 0]<-order.threats[sum.not.na > 0,1]
  # second.threat[[t]]<-rep(NA, nrow(grid@data))
  # second.threat[[t]][sum.not.na > 1]<-order.threats[sum.not.na > 1,2]
  # 
  
  # threat.medians<-apply(grid@data[,primary.threat.columns],2,median, na.rm=T)
  # cumulative.threat[[t]]<-rowSums(t(apply(grid@data[,primary.threat.columns],1,function(x) x >= threat.medians)))
  # na.i<-which(unlist(apply(grid@data[,primary.threat.columns],1,function(x) sum(is.na(x)))) == length(primary.threat.columns))
  # cumulative.threat[[t]][na.i]<-NA
  
  mean.threat[[t]]<-apply(grid@data[,primary.threat.columns],1,mean,na.rm=T)
  # median.threat[[t]]<-apply(grid@data[,primary.threat.columns],1,median,na.rm=T)
  # sum.threat[[t]]<-apply(grid@data[,primary.threat.columns],1,sum,na.rm=T)
  
  output.dir<-paste0(out.dir,t,"/")
  
  # PlotCompositeMap(d = max.threat[[t]],cell.i = which(!is.na(max.threat[[t]])),
  #                  n.breaks = 100,lower.cut = 0,range = c(0,1),
  #                  axis.range = c(0,1),
  #                  fname = paste0(output.dir,t,"_MaxThreatCompositeContinuous.png"),
  #                  cols = color.schemes)
  
  # PlotCompositeMap(d = max.threat[[t]],cell.i = which(!is.na(max.threat[[t]])),
  #                  n.breaks = 2,lower.cut = 0,range = c(0,1),
  #                  axis.range = c(0,1),
  #                  fname = paste0(output.dir,t,"_MaxThreatCompositeFlat.png"),
  #                  cols = color.schemes)
  
  # PlotCompositeMap(d = second.threat[[t]],cell.i = which(!is.na(second.threat[[t]])),
  #                  n.breaks = 100,lower.cut = 0,range = c(0,1),
  #                  axis.range = c(0,1),
  #                  fname = paste0(output.dir,t,"_2ndThreatCompositeContinuous.png"),
  #                  cols = color.schemes)
  # 
  
  # PlotCompositeMap(d = max.threat[[t]],cell.i = which(!is.na(max.threat[[t]])),
  #                  n.breaks = 100,lower.cut = 0,range = c(0,1),
  #                  axis.range = c(0,1),
  #                  fname = paste0(output.dir,t,"_MaxThreatCompositeContinuous2ndTh.png"),
  #                  cols = color.schemes,
  #                  d2 = second.threat[[t]],
  #                  cell.i.2 = which(!is.na(second.threat[[t]])))
  # 
  
  # PlotSingleMap(d = cumulative.threat[[t]],cell.i = which(!is.na(cumulative.threat[[t]])),
  #                  n.breaks = length(primary.threats)+2,lower.cut = 0,range = c(0,length(primary.threats)+1),
  #                  axis.range = c(0,length(primary.threats)+1),
  #                  fname = paste0(output.dir,t,"_CumulativeThreatComposite.png"),
  #                  color.scheme = "Reds")#brewer.pal(length(primary.threats)+2,"Reds")
  
  # PlotSingleMap(d = sqrt(mean.threat[[t]]),cell.i = which(!is.na(mean.threat[[t]])),
  #               n.breaks = 100,lower.cut = 0,
  #               range = c(0,1),#range(mean.threat[[t]],na.rm = T)
  #               axis.range = c(0,1),
  #               fname = paste0(output.dir,t,"_MeanThreatCompositeMagma.png"),
  #               color.scheme = (magma(n = 8,alpha = 0, direction = -1)),sqrt = T)
  # 
  
  # PlotSingleMap(d = median.threat[[t]],cell.i = which(!is.na(median.threat[[t]])),
  #               n.breaks = 100,lower.cut = 0,range = c(0,1),
  #               axis.range = c(0,1),
  #               fname = paste0(output.dir,t,"_MedianThreatComposite.png"),
  #               color.scheme =  "Reds")
  
  
  # PlotSingleMap(d = sum.threat[[t]],cell.i = which(!is.na(sum.threat[[t]])),
  #               n.breaks = 100,lower.cut = 0,range = range(sum.threat[[t]],na.rm=T),
  #               axis.range = range(sum.threat[[t]],na.rm=T),
  #               fname = paste0(output.dir,t,"_SumThreatComposite.png"),
  #               color.scheme =  "Reds")
  
  
  #Plot any threat
  
  
}




max.threat.coverage<-lapply(taxa,function(t) data.frame(table(max.threat[[t]])/sum(!is.na(max.threat[[t]]))))
names(max.threat.coverage)<-taxa
for(t in taxa) max.threat.coverage[[t]]$primary.threat<-primary.threats[as.integer(levels(max.threat.coverage[[t]]$Var1))]
write.csv(do.call(rbind,max.threat.coverage),
          file = "C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/Max.Threat.Coverage.csv",row.names = T)

mean.threat.quantiles<-lapply(mean.threat, function(m) quantile(m,probs = c(0.025,0.5,0.975), na.rm=T))
mean.threat.mean.sd<-lapply(mean.threat, function(m) c(mean(m, na.rm=T),sd(m,na.rm = T)))



breaks<-seq(0,1,0.02)

# png("C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/ThreatPresenceHistograms.png",
#     width = 12,height = 17, units = 'cm',res = 600)
# layout(matrix(1:14,nrow = 7),heights = c(rep(0.155,6),0.07),widths = c(0.5,0.5))
# par(mar = c(0,3,0.25,0))
#for(i in 1:7) plot.new()
threat.data<-list()
for(th in primary.threats)
{
  primary.taxa.columns<-unlist(lapply(taxa[c(2,1,3)],function(t) grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data))))
  
  threat.data[[th]]<-as.matrix(grid@data[,primary.taxa.columns])
}

mean.threat.across.taxa<-unlist(lapply(threat.data, function(d) mean(apply(d,2,median, na.rm=T), na.rm=T)))
mean.threat.across.taxa<-mean.threat.across.taxa[order(-mean.threat.across.taxa)]

library(vioplot)
pdf("C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/ThreatTaxaVioplots.pdf",
    width = 18/2.54,height = 9/2.54)
# png("C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/ThreatTaxaVioplots.png",
#     width = 20,height = 9, units = 'cm',res = 600)
layout(matrix(1:7,ncol = 7),widths = c(0.05,rep(0.95/6,6)))
par(mar = c(0,0,0,0))
plot(NA, xlim = c(0,4),ylim = c(0,1),las=1,axes = F, xlab="",ylab = "", bty='n')
par(mar = c(6,0,0.5,0))
for(th in names(mean.threat.across.taxa))
{
  plot(NA, xlim = c(0,4),ylim = c(0,1),las=1,axes = F, xlab="",ylab = "", bty='n')
  abline(h = seq(0,1,0.25), lty = "dashed", col = "grey")
  
  if(th == names(mean.threat.across.taxa)[1]) axis(side = 2,at = seq(0,1,0.25),las = 1)
  sapply(1:3, function(i)
    vioplot(na.omit(threat.data[[th]][,i]),ylim = c(0,1),add = T, at = i, col = color.schemes[[th]][5],drawRect = T)
  
  )
  axis(side =1, at = 1:3, labels = taxa[c(2,1,3)],cex.axis = 1, las = 2)
  mtext(text = th,line = -2,cex=0.7, font = 2)
}

dev.off()


#Calculate statistics about the proportion of calculable area with gt a given threat probability

for(th in names(mean.threat.across.taxa))
{
  print(paste(th,sapply(1:3, function(i)
    sum(na.omit(threat.data[[th]][,i]) >= 0.5) * 2500/1E6 ))#/length(na.omit(threat.data[[th]][,i]))
  )
}

for(th in names(mean.threat.across.taxa))
{
  print(paste(th,sapply(1:3, function(i)
    (sum(na.omit(threat.data[[th]][,i]) >= 0.5)) /length(na.omit(threat.data[[th]][,i])))))
  
}





#### Plot Figs 1 - 3 ####
#plot Multipanel figure of taxa risk maps
for(t in taxa)
{
  
  print(t)
  output.dir<-paste0(out.dir,t,"/")
  
  t.i<-which(taxa == t)
  
  t.range<-c(0,1)
  
  pdf(paste0(paste0(output.dir,"MultipanelPThreat_SRgt10sp_MaxThreat"),".pdf"),height = 21/2.54, width = 18/2.54)
  #png(paste0(paste0(output.dir,"MultipanelPThreat_SRgt10sp_MaxThreat"),".png"),height = 20, width = 20,units = "cm",res = 600)

  #layout(matrix(c(1:3,7,4:6,8),ncol = 2,byrow = F),heights = rep(0.25,4))#Used for the 1st submission figures
  layout(matrix(c(1,2,7,8,3,4,9,10,5,6,11,12,13,13,14,15),ncol = 4,byrow = T),
         heights = rep(0.25,4),widths = c(0.4,0.1,0.4,0.1))#Layout with uncertainty figure  
    
  
  i = 1
  for(th in primary.threats)
  {
    print(th)
    colname<-grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data), value = T)
    
    gt10.col.i<-grep(paste0(t,"_SRgt10sp"),colnames(grid@data))
    gt10sp.i<-which(!is.na(grid@data[[gt10.col.i]]))
    
    
    PlotSingleMapInMultiPanel(d = sqrt(grid@data[[colname]]),
                              cell.i = gt10sp.i,
                              n.breaks = 100,
                              lower.cut = 0,
                              range = t.range,
                              axis.range = t.range,
                              color.scheme = color.schemes[[th]],i = i,sqrt = T,
                              axis = T,horiz = F)
    mtext(side = 1,text = th,line = 0, adj = 1.02,cex = 1)
    i = i+1
  }
  
  PlotCompositeMapInMultiPanel(d = max.threat[[t]],cell.i = which(!is.na(max.threat[[t]])),
                   n.breaks = 2,lower.cut = 0,range = c(0,1),
                   axis.range = c(0,1),
                   cols = color.schemes)
  mtext(text = letters[i],side = 3,adj = 0.1,line = -2, cex = 1)
  i = i + 1
  
  PlotSingleMapInMultiPanel(d = sqrt(apply(grid_DDCells_Uncert[,t.i,,6]-grid_DDCells_Uncert[,t.i,,2],1,median,na.rm=T)),
                            cell.i = gt10sp.i,
                            n.breaks = 100,
                            lower.cut = 0,
                            range = sqrt(c(0,0.4)),
                            axis.range = c(0,0.4),
                            color.scheme = viridis(8),i = i,sqrt = T,
                            axis = T, horiz=F)
  mtext(side = 1,text = "Uncertainty",line = 0, adj = 1.02, cex = 1)
  i = i+1
  
  # PlotAxisInMultiPanel(threats.to.use = rev(primary.threats),
  #                      n.breaks = 100,lower.cut = 0,range = t.range,
  #                      axis.range = t.range,cols = color.schemes,sqrt = T)
  # #mtext(text = letters[i],side = 3,adj = 0.1,line = -2)
  
  dev.off()
  
}

####







# Calculate the SR weighted threat occurence layer -----
for(t in taxa)
{
  print(t)
  output.dir<-paste0(out.dir,t,"/")
  if(!dir.exists(output.dir)) dir.create(output.dir)
  
  sr.col.i<-grep(paste0(t,"_SR"),colnames(grid@data))[1]
  global.i<-which(!is.na(grid@data[[sr.col.i]]))
  
  
  gt10.col.i<-grep(paste0(t,"_SRgt10sp"),colnames(grid@data))
  gt10sp.i<-which(!is.na(grid@data[[gt10.col.i]]))
  
  
  range<-c(0,
           ceiling(max(
             sapply(primary.threats, 
                    function(th) 
                      grid@data[,grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data), value = T)]*grid@data[,gt10.col.i]), na.rm=T)))
  
  
  for(th in primary.threats)
  {
    print(th)
    colname<-grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data), value = T)
    
    grid@data[[paste0(colname,"_SRgt10sp")]]<-
      grid@data[,colname]*grid@data[,gt10.col.i]
   
    # PlotSingleMap(d = grid@data[[paste0(colname,"_SRgt10sp")]],
    #               cell.i = gt10sp.i,
    #               n.breaks = 100,
    #               lower.cut = 0,
    #               range = range,
    #               axis.range = range,
    #               fname = paste0(paste0(output.dir,colname,"_SRgt10sp"),".png"),
    #               color.scheme = color.schemes[[th]]) 
  }

  
  th = "AnyThreat"
  print(th)
  colname<-grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data), value = T)
  
  grid@data[[paste0(colname,"_SRgt10sp")]]<-
    grid@data[,colname]*grid@data[,gt10.col.i]
  
  # PlotSingleMap(d = grid@data[[paste0(colname,"_SRgt10sp")]],
  #               cell.i = gt10sp.i,
  #               n.breaks = 100,
  #               lower.cut = 0,
  #               range = range(grid@data[[paste0(colname,"_SRgt10sp")]],na.rm = T),
  #               axis.range = range(grid@data[[paste0(colname,"_SRgt10sp")]], na.rm=T),
  #               fname = paste0(paste0(output.dir,colname,"_SRgt10sp"),".png"),
  #               color.scheme = "YlOrRd")
  
  
}



#plot Multipanel figure of taxa risk maps
for(t in taxa)
{
  print(t)
  output.dir<-paste0(out.dir,t,"/")
  
  t.range<-c(0,
           ceiling(max(
             sapply(primary.threats, 
                    function(th) 
                      grid@data[,grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)]), na.rm=T)))
  
  pdf(paste0(paste0(output.dir,"MultipanelThreat_SRgt10sp"),".pdf"),height = 21/2.54, width = 18/2.54)
  #png(paste0(paste0(output.dir,"MultipanelThreat_SRgt10sp"),".png"),height = 20, width = 20,units = "cm",res = 600)
  layout(matrix(1:12,ncol = 2,byrow = F),heights = rep(c(0.85,0.15)/3,3))
  #layout(matrix(1:6,ncol = 2,byrow = F))
  i = 1
  for(th in primary.threats)
  {
    print(th)
  
    colname<-grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)
  
    gt10.col.i<-grep(paste0(t,"_SRgt10sp"),colnames(grid@data))
    gt10sp.i<-which(!is.na(grid@data[[gt10.col.i]]))
    
      
    PlotSingleMapInMultiPanel(d = grid@data[[colname]],
                  cell.i = gt10sp.i,
                  n.breaks = 100,
                  lower.cut = 0,
                  range = t.range,
                  axis.range = t.range,
                  color.scheme = color.schemes[[th]],i = i)
    i = i+1
  }
  dev.off()
  
}



#Write out shape files of the risk maps
# for(c in grep("gt10sp_SRgt10sp",colnames(grid@data), value = T))
#   writeOGR(obj = grid[,c],
#            dsn = "C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic",
#            layer = c,
#            driver = "ESRI Shapefile",overwrite_layer = T)



##### Plot Fig 4 - cumulative risk map #####

for(th in primary.threats)
{
  print(th)
  colname<-grep(paste0("_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)
  grid@data[paste0("All_",th,"_SRRisk")]<-rowSums(grid@data[,colname], na.rm=T)
}

range<-c(0,
         ceiling(max(
           sapply(primary.threats, 
                  function(th) 
                    grid@data[,grep(paste0("All_",th,"_SRRisk"),colnames(grid@data), value = T)]), na.rm=T)))


cmyk <- function(C,M,Y,K) {
  
  C <- C / 100.0
  M <- M / 100.0
  Y <- Y / 100.0
  K <- K / 100.0
  
  n.c <- (C * (1-K) + K)
  n.m <- (M * (1-K) + K)  
  n.y <- (Y * (1-K) + K)
  
  r.col <- ceiling(255 * (1-n.c))
  g.col <- ceiling(255 * (1-n.m))
  b.col <- ceiling(255 * (1-n.y))
  
  return(col2rgb(sprintf("#%02s%02s%02s",
                         as.hexmode(r.col), 
                         as.hexmode(g.col), 
                         as.hexmode(b.col))))
  
}

# png("C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/CMYK_MultipanelRiskPriorities_summary-rev.png",
#     height = 20, width = 20,units = "cm",res = 600)
pdf("C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/CMYK_MultipanelRiskPriorities_summary-rev.pdf",
    height = 21/2.54, width = 18/2.54)
{
layout(matrix(c(rep(c(1:3,7,8),2),
                c(4:6,7,8),
                c(4:6,9,10)),ncol = 4,byrow = F),heights = c(rep(0.21,3),0.3,0.07))
#Previous summary layout
#layout(matrix(c(1:3,7:8,4:6,9:10),ncol = 2,byrow = F),heights = c(rep(0.21,3),0.31,0.06))

i = 1

par(oma=c(0,0,0,0))
par(mar=c(0,0,0,0))
par(xpd=NA)
par(cex=1)

hotspot.count<-rep(0,nrow(grid@data))

for(th in primary.threats)
{
  print(th)

  
  q.t<-lapply(taxa,function(t) quantile(grid@data[,grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)],
                                        probs= c(0.9),na.rm=T))
  names(q.t)<-taxa
  hotspots<-lapply(taxa,function(t)
    grid@data[,grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)] >= q.t[[t]])
  
  hotspots<-do.call(cbind, hotspots)
  # p.taxa<-do.call(cbind, p.taxa)
  # p.max<-do.call(cbind,p.max)
  # use.i<-which(apply(p.max,1, function(x) sum(is.na(x))) < 3)
  # 
  # p.max[use.i,][is.na(p.max[use.i,])]<-0
  # 
  use.i<-which(apply(hotspots,1, function(x) sum(x, na.rm=T)) > 0 )
  hotspots[use.i,][is.na(hotspots[use.i,])]<-0
  
  hotspot.count<-hotspot.count+rowSums(hotspots, na.rm=T)
  
  cmyk.t<-cmyk(80*hotspots[use.i,1],100*hotspots[use.i,2],80*hotspots[use.i,3],20)/255
  
  plot(WL.moll,lwd=0.05,col = "lightgrey", border = "lightgrey")
  plot(grid[use.i,],
       col = rgb(cmyk.t[1,],cmyk.t[2,],cmyk.t[3,]),border = NA,add=T)

  mtext(side = 1,text = th,line = -2, adj = 0.6)
  mtext(text = letters[i],cex=1.1,side = 3,adj = 0.1,line = -1.5)
  i = i+1
  
}

max.count<-max(hotspot.count)
use.i<-which(hotspot.count > 0)
plot(WL.moll,lwd=0.05,col = "lightgrey", border = "lightgrey")
plot(grid[use.i,],
     col = rev(plasma(n = max.count))[hotspot.count[use.i]],border = NA,add=T)
mtext(text = letters[i],cex=1.1,side = 3,adj = 0.1,line = -1.5)
i = i+1
t.breaks = seq(1,max.count,length.out = max.count)#c(0,seq(lower.cut,1,length.out = n.breaks))
par(mar=c(2,6,0,6))

x.lim = c(1,max.count)
y.lim = c(0,1)

plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
y = 1
sapply(seq(1,max.count),
       function(i) {
         rect((t.breaks[i]), y-1, (t.breaks[i+1]), y, col=rev(plasma(n = max.count))[i], border=NA)
         #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
       })

axis(1,padj=-1, at = t.breaks, cex.axis=0.7)

dev.off()
}

#Calculate the proportion of calculable area that is covered by a hotspot

#Find cells that have more that 10sp of mammals, birds or amphibians
modellable.cells<-sum(rowSums(grid@data[,c("mammal_SRgt10sp","bird_SRgt10sp","amph_SRgt10sp")], na.rm = T) > 10)
sum(hotspot.count > 0)/modellable.cells

#proportion of the terrestrial surface
sum(hotspot.count > 0)/sum(table(grid@data$REALM))


  
  
for(th in primary.threats)
{
  print(th)
  
  # PlotSingleMap(d = grid@data[[grep(paste0(t,"_SRgt10sp"),colnames(grid@data))]],
  #               cell.i = gt10sp.i,
  #               n.breaks = 100,
  #               lower.cut = 0,
  #               range = range(grid@data[[grep(paste0(t,"_SRgt10sp"),colnames(grid@data))]],na.rm=T),
  #               axis.range = range(grid@data[[grep(paste0(t,"_SRgt10sp"),colnames(grid@data))]],na.rm=T),
  #               fname = paste0(output.dir,t,"_SR.png"),color.scheme = "Reds")
  
  
  # PlotSingleMap(d = sqrt(grid@data[,paste0("All_",th,"_SRRisk")]),
  #               cell.i = which(grid@data[,paste0("All_",th,"_SRRisk")] > 0.1),
  #               n.breaks = 100,
  #               lower.cut = 0,
  #               range = sqrt(range),
  #               axis.range = (range),
  #               fname = paste0("C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/All",th,"SRRisk_sqrt.png"),
  #               color.scheme = color.schemes[[th]], sqrt = T)
  
  # p.taxa<-lapply(taxa,function(t)
  #   grid@data[,grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)]/
  #     grid@data[,paste0("All_",th,"_SRRisk")])
  # 
  # p.max<-lapply(taxa,function(t)
  #   grid@data[,grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)]/
  #     max(grid@data[,grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)], na.rm = T))
  
  #paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data) #Risk map data
  
  q.t<-lapply(taxa,function(t) quantile(grid@data[,grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)],
                                        probs= c(0.9),na.rm=T))
  names(q.t)<-taxa
  hotspots<-lapply(taxa,function(t)
    grid@data[,grep(paste0(t,"_",th,"gt10sp_SRgt10sp$"),colnames(grid@data), value = T)] >= q.t[[t]])
  
  hotspots<-do.call(cbind, hotspots)
  # p.taxa<-do.call(cbind, p.taxa)
  # p.max<-do.call(cbind,p.max)
  # use.i<-which(apply(p.max,1, function(x) sum(is.na(x))) < 3)
  # 
  # p.max[use.i,][is.na(p.max[use.i,])]<-0
  # 
  use.i<-which(apply(hotspots,1, function(x) sum(x, na.rm=T)) > 0 )
  hotspots[use.i,][is.na(hotspots[use.i,])]<-0
  
  
  png(paste0("C:/Users/mikeha/Dropbox/Threat mapping/Occupancy modelling/RegionalLogistic/CMYK_",th,"_LikelihoodHotspots.png"),
      width = 17,height = 12, units = 'cm',res = 600)
  
  par(oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  par(xpd=NA)
  par(cex=1)
  
  cmyk.t<-cmyk(80*hotspots[use.i,1],100*hotspots[use.i,2],80*hotspots[use.i,3],20)/255
  
  layout(matrix(c(1,2),nrow = 2), heights = c(0.9,0.1))
  plot(WL.moll,lwd=0.05,col = "lightgrey", border = "lightgrey")
  plot(grid[use.i,],
       col = rgb(cmyk.t[1,],cmyk.t[2,],cmyk.t[3,]),border = NA,add=T)

  dev.off()
  
}


lapply(primary.threats,function(th) quantile(grid@data[,paste0("All_",th,"_SRRisk")], probs=c(0.025,0.5,0.975), na.rm=T))

for(th in primary.threats)
{
  temp<-grid@data[,paste0("All_",th,"_SRRisk")]
  #print(c(th, (quantile(temp[temp > 0], probs=c(0.025,0.5,0.975), na.rm=T))))
  print(c(th, (mean(temp[temp > 0], na.rm=T)),(sd(temp[temp > 0], na.rm=T))))
}


# Fig 5.
# Plot Heat map of threat vs HFI ####


#Extract the median HFI for each grid cell
#rHFI<-raster("C:/Users/mikeha/Work/Spatial data/Dryadv3/Maps/HFP2009.tif")
#Pre calculated zonal statistics
rHFI.grid.data<-read.csv("C:/Users/mikeha/Dropbox/Threat mapping/SpatialHumanActivityData/HFI/ZonalSt_shp2_TableToExcel.csv")

#grid.hfi<-extract(rHFI,grid,fun= median,na.rm=T,df=T)


#Add HFI data to the grid


grid@data$HFI2009<-rHFI.grid.data$MEAN[match(row.names(grid@data), rHFI.grid.data$FID)]

max.HFI = 50
grid@data$HFI2009.scaled<-grid@data$HFI2009/max.HFI

library(MASS)


HFI.threats<-list(
  'Land use' = c("Logging","Agriculture"),
  Hunting = c("Hunting"),#,"Climate change","Pollution","Invasives"),
  Climate = "Climate change",
  Invasives = "Invasives",
  Pollution = "Pollution"
)

HFI.taxa<-list(
  'Land use' = c("mammal","amph","bird"),
  Hunting = c("mammal","amph","bird"),
  Climate = "bird",
  Invasives = c("mammal","amph","bird"),
  Pollution = c("mammal","amph","bird")
)


HFI.means<-lapply(1:length(HFI.threats),function(i)
{
  c<-unlist(sapply(HFI.threats[[i]],function(th) 
    lapply(HFI.taxa[[i]], function(t) 
      grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data), value = T))))
  if(length(c) > 1)
  {
    rowMeans(grid@data[,c],na.rm=T)
  } else
  {
    grid@data[,c]
  }
  })#,1,function(x) mean(x,na.rm=T)))




lapply(taxa, function(t) lapply(primary.threats, function(th) median(grid@data[,grep(paste0(t,"_",th,"gt10sp$"),colnames(grid@data), value = T)], na.rm=T)))
primary.threats


gt10sp.i<-lapply(HFI.means, function(d) which(apply(cbind(grid@data$HFI2009.scaled,d),1,function(x) sum(is.na(x))) == 0))

#Using Kernel Density Estimator
h1 <- hist(grid@data$HFI2009.scaled, breaks=100, plot=F)
h2 <- lapply(HFI.means,function(d) hist(d, breaks=100, plot=F))
taxa.top <- max(unlist(lapply(h2,function(x) max(x$counts))))

k <- lapply(1:length(HFI.means),function(t) {
  #gt10sp.i<-which(apply(cbind(grid@data$HFI2009,mean.threat[[t]]),1,function(x) sum(is.na(x))) == 0)
  kde2d(grid@data$HFI2009.scaled[gt10sp.i[[t]]], HFI.means[[t]][gt10sp.i[[t]]], n=500, lims = c(0,1,0,1))
  })


for(t in 1:length(k)) k[[t]]$z<-sqrt(k[[t]]$z)

z.range<-c(1E-4,max(unlist(lapply(k, function(x) max(x$z)))))
z.range =c(1E-4,9.05)



r<-rev(viridis(50))
#r <- rev(brewer.pal(11,'Spectral'))#viridis(n = 32)

# margins
oldpar <- par()
output.dir<-paste0(out.dir,t,"/")

lims<-c(0,1.0)


#Calculate linear relationship, accounting for spatial autocorrelation
coords<-coordinates(gCentroid(grid,byid = T))/1E3

#which grid cells to consider
use.cells<-lapply(HFI.means,function(x) which(!is.na(x)))



#png(paste0(out.dir,"Multipanel_ThreatvsHFI_Heatmap_.png"), width = 20, height = 20, units = "cm", res = 600)
pdf(paste0(out.dir,"Multipanel_ThreatvsHFI_Heatmap_.pdf"), width = 18/2.54, height = 21/2.54)
{
layout(matrix(
  c(rep(1,3),2,3:10), byrow=F,ncol=3),heights = c(rep(0.31,3),0.07),widths = c(0.1,0.35,0.55))

par(mar=c(6,4,6,1))
plot(NA, axes=F, xlab="",ylab="",xlim = c(0,1), ylim = c(0,z.range[2]^2))

sqrt<-T

t.breaks<-seq(0,z.range[2],length.out = length(r))
sapply(seq(1,length(r)),
       function(i) {
         rect( 0,(t.breaks[i])^2,1, (t.breaks[i+1])^2, col=r[i], border=NA)
         #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
       })
axis(2, cex.axis=1.0,
     at = pretty(t.breaks^2),las =1)#,padj=-2
mtext(side = 2, text = "Pixel count",line = 2.5, cex=0.9)

par(mar=c(1,4.5,1,0))
plot.new()

par(mar=c(1,4,1,0))

for(t in 1:3)
{

  #png(paste0(out.dir,"ThreatvsHFI_Heatmap_",names(HFI.threats)[t],".png"), width = 15, height = 15, units = "cm", res = 600)
  {
    image(k[[t]], col=r, zlim = z.range,ylim = lims,xlim=lims,las = 1,axes=F) #plot the image
    axis(side = 2,las=1)
    axis(side = 1,labels = F)
    abline(a = 0, b = 1, col = "grey",lwd=2)
    mtext(side=2,bquote(paste('P'['Th']*'(', .(names(HFI.threats)[t]),')')),line = 2.75,cex = 0.9)
    mtext(side = 3,letters[t],cex=1.1,line = -0.5,adj = -0.3)
  }
  #dev.off()

}
axis(side = 1)
mtext(side = 1, "HFI",cex=0.9,line = 2.75)

plot.new()

outlier.distance<-lapply(HFI.means,function(d) d-grid@data$HFI2009.scaled)


lapply(1:3, function(i)
  
  PlotSingleMapInMultiPanel(d = outlier.distance[[i]],sqrt = F,
                cell.i = gt10sp.i[[i]],i = 3+i,axis = F,axis.only = F,
                n.breaks = 100,
                lower.cut = 0,
                range = c(-1,1),
                axis.range = c(-1,1),
                color.scheme = brewer.pal(n = 9,"RdBu"))
)

#par(mar=c(3,5,0,5))
PlotSingleMapInMultiPanel(d = outlier.distance[[1]],sqrt = F,
                          cell.i = gt10sp.i[[1]],i = 3+i,axis = T,axis.only = T,
                          n.breaks = 100,
                          lower.cut = 0,
                          range = c(-1,1),
                          axis.range = c(-1,1),
                          color.scheme = brewer.pal(n = 9,"RdBu"))

mtext(side = 1,bquote(paste('P'['Th']*' - HFI')),line = 2,cex = 0.9)

dev.off()
}


png(paste0(out.dir,"ThreatvsHFI_Heatmap_legend.png"), width = 4, height = 15, units = "cm", res = 600)
{

  plot(NA, axes=F, xlab="",ylab="",xlim = c(0,1), ylim = c(0,z.range[2]^2))

sqrt<-T
if(sqrt)
{
  t.breaks<-seq(0,z.range[2],length.out = length(r))
  sapply(seq(1,length(r)),
         function(i) {
           rect( 0,(t.breaks[i])^2,1, (t.breaks[i+1])^2, col=r[i], border=NA)
           #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
         })
  
} else {
  sapply(seq(1,n.breaks),
         function(i) {
           rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
           #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
         })
  
}

axis(2, cex.axis=1.0,
     at = pretty(t.breaks^2),las =1)#,padj=-2

dev.off()
}





