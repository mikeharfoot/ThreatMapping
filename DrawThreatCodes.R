library(reshape2)

path = 'C:/Users/mikeha/Dropbox/Threat mapping/'
source(paste0(path,'scripts/ThreatMapFunctions.R'))

GetSpThreatCodes<-function(threat.species.info,threat.codes, threats)
{
  taxa.th.sp.ids = list()
  
  for(th in names(threat.codes))
  {
    th.sp.ids = list()
    for(code in threat.codes[[th]])
    {
      code.col = grep(paste("^X",code,sep=""),colnames(threats))
      if(length(code.col) > 0) th.sp.ids = c(th.sp.ids,threats$Species.ID[threats[,code.col] > 0])
      
    }
    th.sp.ids = unique(unlist(th.sp.ids))
    taxa.th.sp.ids[[th]] = th.sp.ids
  }
  taxa.th.sp.ids
  
  th.sp.binomials = list()
  for(th in names(threat.codes))
  {
    th.sp.binomials[[th]] = threat.species.info$binomial[threat.species.info$Species.ID %in% taxa.th.sp.ids[[th]]]
  }
  
  th.sp.binomials
}



i<-read.csv("C:/Users/mikeha/Work/Spatial data/Red List/2017_3/Intersections/50km/test_result2018.csv", stringsAsFactors = F)
i['count']<-1


#Read the mammal area-intersections
load("C:/Users/mikeha/Work/Spatial data/Red List/2017_3/MAMMALS/RangeIntersections/intersects_master.u.R")
mammal.a.array = intersects

#Read the mammal range areas
mammal.range.sizes<-read.csv("C:/Users/mikeha/Work/Spatial data/Red List/2017_3/MAMMALS/RangeIntersections/mammal_range_area_master.u.csv")


# Read the Amphibians data --------


#Read the area-intersections
load("C:/Users/mikeha/Work/Spatial data/Red List/2017_3/AMPHIBIANS/RangeIntersections/intersects_master.u.R")
amphibian.a.array<-intersects

#Read the range areas
amphibian.range.sizes<-read.csv("C:/Users/mikeha/Work/Spatial data/Red List/2017_3/AMPHIBIANS/RangeIntersections/AMPHIBIANS_range_area_master.u.csv")




#read grid.dbf to find cell ids in each region
grid<-foreign::read.dbf("C:/Users/mikeha/Dropbox/Threat mapping/Regional/GlobalGrid50kmRoundl_Clip_Ecoregion.dbf",as.is = T)

#read bird range data to get sp_id and binomial data
aves<-read.table(file = "ArcGIS/Projects/ThreatMapIntersection/AVES_extent_no_passage.txt",sep = ',',header =T, stringsAsFactors = F)

load(file = "C:/Users/mikeha/Dropbox/Threat mapping/Regional 2017_3/All_ProportionDD.R")


realms = c(
  'Afrotropic',
  'Indomalayan',
  'Nearctic',
  'Neotropic',
  'Palearctic',
  'Australasia_Oceania'
)


global.i = NULL

grid$Prop.sp.dd<-NA

for(r in realms){
  
  if(r == "Australasia_Oceania")
  {
    r.split<-strsplit("Australasia_Oceania","_")[[1]]
    grid$Prop.sp.dd[which((grid$REALM == r.split[1]) | (grid$REALM == r.split[2]))]<-prop.sp.dd[[r]]
  } else
  {
    grid$Prop.sp.dd[which(grid$REALM == r)]<-prop.sp.dd[[r]]
  }
  
  # if(r == "Australia_Oceania")
  # {
  #   r.split<-strsplit("Australasia_Oceania","_")[[1]]
  #   global.i<-c(global.i,which((grid$REALM == r.split[1]) | (grid$REALM == r.split[2])))
  # } else
  # {
  #   global.i<-c(global.i,which(grid$REALM == r))
  # }
}




sp<-unique(i$ID_NO)

probs<-c(0,0.025,0.25,0.5,0.75,0.975,1.0)

#DD stats over bird species
q.sp.dd.aves<-sapply(sp, function(s) {
  quantile(grid$Prop.sp.dd[i$WDPAID[i$ID_NO == s]],probs = probs,na.rm=T) 
})

df.sp.th<-data.frame(binomial = aves$binomial[match(sp,aves$sp_id)], t(q.sp.dd),stringsAsFactors = F)


#DD stats for mammal species
q.sp.dd.mammals<-apply(mammal.a.array, 2, function(s){
  quantile(grid$Prop.sp.dd[s],probs = probs, na.rm=T)
})
df.sp.th<-rbind(df.sp.th,data.frame(binomial = mammal.range.sizes$binomial,t(q.sp.dd.mammals)))


#DD stats for mammal species
q.sp.dd.amphibians<-apply(amphibian.a.array, 2, function(s){
  quantile(grid$Prop.sp.dd[s],probs = probs, na.rm=T)
})
df.sp.th<-rbind(df.sp.th,data.frame(binomial = amphibian.range.sizes$binomial,t(q.sp.dd.amphibians)))



#for each threat, draw a set of threat codes for the species using the dd data
ThreatCodes = GetThreatCodesPhases()


threats = read.csv(paste0(path,"Birds_Mammals_Amphibians_Threats_Cleaned.csv"), stringsAsFactors = F)
threat.species.info = read.csv(paste0(path,"Birds_Mammals_Amphibians_Threats_Cleaned_Species_Info.csv"), stringsAsFactors = F)
#Add a column for the binomial name
threat.species.info$binomial = paste(threat.species.info$Genus,
                                     threat.species.info$Species)



# #cut threats to birds only
# aves.species.id<-threat.species.info$Species.ID[match(aves$binomial,threat.species.info$binomial)]
# threat.aves<-threats[match(aves.species.id,threats$Species.ID),]
# threat.species.info.aves<-


phase = 'phase1'


orig.threat.codes<-GetSpThreatCodes(threat.species.info,threats = threats,
               threat.codes = ThreatCodes[[phase]])


n.draws<-100

p.u.zero<-0.005

drawn.threat.codes<-list()

for(th.i in 1:length(ThreatCodes[[phase]]))
{
  th <- ThreatCodes[[phase]][th.i]
  
  
  df.sp.th[[paste0("th.o.",names(th))]]<-0
  df.sp.th[[paste0("th.o.",names(th))]][match(orig.threat.codes[[names(th)]],df.sp.th$binomial)]<-1
  
  threat.code.draw<-sapply(seq_len(n.draws), function(d){
    th.d<-df.sp.th[[paste0("th.o.",names(th))]]
    to.draw.new.codes<-runif(n = nrow(df.sp.th)) < (p.u.zero + ((1-p.u.zero)*df.sp.th$X50.))
    th.d[which(to.draw.new.codes)] <-1-th.d[which(to.draw.new.codes)]
    th.d
  })
  
  for(c in seq_len(n.draws)) df.sp.th[[paste0("Draw.",c)]]<-threat.code.draw[,c]
  
  drawn.threat.codes[[names(th)]]<-df.sp.th
}


#Save the threat codes
save(drawn.threat.codes,file = paste0(path,"ThreatCodesDrawnWithUncertainty.csv"))

