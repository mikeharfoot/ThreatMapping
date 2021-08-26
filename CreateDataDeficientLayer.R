wants <- c("rgeos","rgdal","maptools","plyr","raster","RColorBrewer","reshape2","xlsx")
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
  'Palearctic',
  'Australasia_Oceania'
)



#Read the attribute table from each taxa
attributes.dir<-"C:/Users/mikeha/Documents/ArcGIS/Projects/ThreatMapIntersection/"
table.str<-"_TableToExcel.xls"



taxa = c('mammal','amph','bird')

RL.taxa<-list(
  mammal = "MAMMALIA",
  amph = "AMPHIBIA",
  bird = "AVES"
)

att.tables<-list()
dd.sp<-list()

for(t in names(RL.taxa))
{
  att.tables[[t]]<-read.xlsx(list.files(path = attributes.dir,pattern = paste0(RL.taxa[[t]],table.str),full.names = T),sheetIndex = 1,stringsAsFactors = F)
  att.tables[[t]]<-att.tables[[t]][(att.tables[[t]]$presence == 1),]
  dd.sp[[t]]<-att.tables[[t]]$binomial[att.tables[[t]]$code == "DD"]
}

lapply(att.tables,function(d) table(d$code))



sp.r<-list()
dd.sp.r<-list()
prop.sp.dd<-list()

grid@data$Prop.sp.dd<-NA

#Loop over realms
for(r in realms)
{
  print(r)
  for (t in taxa)
  {
    print(t)
    load(paste0(path,"Regional 2017_3/",r,"/",t,"_intersection.a.sp.R"))
    
    if(t == 'mammal')
    {
      inter.a = r.mammal.a
    } else if(t == 'bird')
    {
      inter.a = r.bird.a
      inter.a[inter.a == 1]<-TRUE
      inter.a[is.na(inter.a)]<-FALSE
      for(c in 1:ncol(inter.a)) inter.a[,c]<-as.logical(unlist(inter.a[,c]))
    } else if(t == 'amph')
    {
      inter.a = r.amphibian.a
    }
    
    #count data deficient species per cell
    spc<-rowSums(inter.a,na.rm=T)
    if(is.null(sp.r[[r]]))
      sp.r[[r]] <-spc
    else
    {
      if(length(spc) != length(sp.r[[r]])) print(paste(t,r,length(sp.r[[r]]), length(spc)))
      sp.r[[r]] <- sp.r[[r]]+spc
    }
    
    dd.inter.a<-inter.a[,which(colnames(inter.a) %in% dd.sp[[t]])]
    dd.spc<-rowSums(dd.inter.a,na.rm=T)
    if(is.null(dd.sp.r[[r]]))
      dd.sp.r[[r]]<-dd.spc
    else{
      
      if(length(dd.spc) != length(dd.sp.r[[r]])) print(paste("DD",t,r,length(dd.sp.r[[r]]), length(dd.spc)))
      dd.sp.r[[r]]<-dd.sp.r[[r]]+dd.spc
    }
    
    
    prop.dd<-dd.spc/spc
    
    save(spc,file = paste0(path,"Regional 2017_3/",r,"/",t,"_SpOccCount.R"))
    save(dd.spc,file = paste0(path,"Regional 2017_3/",r,"/",t,"_DDSpOccCount.R"))
    save(prop.dd,file = paste0(path,"Regional 2017_3/",r,"/",t,"_ProportionDD.R"))
  }
  
  
  prop.sp.dd[[r]]<-dd.sp.r[[r]]/sp.r[[r]]

  if(r == "Australasia_Oceania")
  {
    r.split<-strsplit("Australasia_Oceania","_")[[1]]
    grid@data$Prop.sp.dd[which((grid$REALM == r.split[1]) | (grid$REALM == r.split[2]))]<-prop.sp.dd[[r]]
  } else
  {
    grid@data$Prop.sp.dd[which(grid$REALM == r)]<-prop.sp.dd[[r]]
  }
}






print("Saving predictions")
save(sp.r,file = paste0(path,"Regional 2017_3/All_SpOccCount.R"))
save(dd.sp.r,file = paste0(path,"Regional 2017_3/All_DDSpOccCount.R"))
save(prop.sp.dd,file = paste0(path,"Regional 2017_3/All_ProportionDD.R"))


global.i = NULL

for(r in realms){
  if(r == "Australasia_Oceania")
  {
    r.split<-strsplit("Australasia_Oceania","_")[[1]]
    global.i<-c(global.i,which((grid$REALM == r.split[1]) | (grid$REALM == r.split[2])))
  } else
  {
    global.i<-c(global.i,which(grid$REALM == r))
  }
}


PlotSingleMap(d = grid@data$Prop.sp.dd,
              cell.i = global.i,
              n.breaks = 100,
              lower.cut = 0,
              range = c(0,0.07),
              axis.range = c(0,0.07),
              fname = "C:/Users/mikeha/Dropbox/Threat mapping/Uncertainty/DataDeficients.pdf",
              color.scheme = brewer.pal(8,"Reds"))
