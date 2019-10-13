#### Regional model #####

path = 'C:/Users/mikeha/Dropbox/Threat mapping/'
source(paste0(path,'scripts/ThreatMapFunctions.R'))
#taxa = c('mammal','bird','amph','All')
taxa = c('mammal','amph','bird')

realms = c(
  'Afrotropic',
  'Indomalayan',
  'Nearctic',
  'Neotropic',
  'Palearctic',
  'Australasia_Oceania'
)


ThreatCodes = GetThreatCodesPhases()


threats = read.csv(paste0(path,"Birds_Mammals_Amphibians_Threats_Cleaned.csv"), stringsAsFactors = F)
threat.species.info = read.csv(paste0(path,"Birds_Mammals_Amphibians_Threats_Cleaned_Species_Info.csv"), stringsAsFactors = F)
#Add a column for the binomial name
threat.species.info$binomial = paste(threat.species.info$Genus,
                                     threat.species.info$Species)

phase = 'phase1'


#Loop over taxonomic group
for (t in taxa)
{
  
  #Loop over threats
  for(th.i in 1:length(ThreatCodes[[phase]]))
  #Run for logging only first
  #th = ThreatCodes[[phase]][3]
  {
    th <- ThreatCodes[[phase]][th.i]
    #Loop over realms
    for(r in realms)
    {
      print(r)
      if(t == "All")
      {
        
        load(paste0(path,"Regional 2017_3/",r,"/","mammal","_intersection.a.sp.R"))
        load(paste0(path,"Regional 2017_3/",r,"/","mammal","_ranges.R"))
        inter.a = r.mammal.a
        ranges = r.mammal.ranges
        load(paste0(path,"Regional 2017_3/",r,"/","bird","_intersection.a.sp.R"))
        load(paste0(path,"Regional 2017_3/",r,"/","bird","_ranges.R"))
        inter.a<-cbind(inter.a,r.bird.a)
        ranges<-rbind(ranges,r.bird.ranges)
        load(paste0(path,"Regional 2017_3/",r,"/","amph","_intersection.a.sp.R"))
        load(paste0(path,"Regional 2017_3/",r,"/","amph","_ranges.R"))
        inter.a<-cbind(inter.a,r.amphibian.a)
        ranges<-rbind(ranges,r.amphibian.ranges)
        
      }
      else
      {
        load(paste0(path,"Regional 2017_3/",r,"/",t,"_intersection.a.sp.R"))
        load(paste0(path,"Regional 2017_3/",r,"/",t,"_ranges.R"))
      
        if(t == 'mammal')
        {
          inter.a = r.mammal.a
          ranges = r.mammal.ranges
        } else if(t == 'bird')
        {
          inter.a = r.bird.a
          inter.a[inter.a == 1]<-TRUE
          inter.a[is.na(inter.a)]<-FALSE
          for(c in 1:ncol(inter.a)) inter.a[,c]<-as.logical(unlist(inter.a[,c]))
          ranges = r.bird.ranges
          colnames(ranges)[2]<-"RangeArea"
        } else if(t == 'amph')
        {
          inter.a = r.amphibian.a
          ranges = r.amphibian.ranges
        }
      }
      r.th.sp.bin = ByCellThreatenedSpeciesBinomials(threat.species.info = threat.species.info,
                                                     range.sizes = ranges,
                                                     threats = threats,
                                                     threat.codes = ThreatCodes[[phase]])
      
      
      occ.matrices = ConstructOccupancyMatricesSmall(r.array = inter.a,
                                                     r.ranges = ranges,
                                                     th.sp.binomials = r.th.sp.bin[[names(th)]])
      
      
      nsite = nrow(occ.matrices$threat.code.matrix)
      
      r.models = vector(mode = 'list',length = nsite)
      
      
      print("Fitting models")
      for(i in 1:nsite)
      {
        if(i%%floor(0.1*nsite) ==0) print(i)
        if(sum(!is.na(occ.matrices$threat.code.matrix[i,])) > 0)
        {
          r.models[[i]] = glm(occ.matrices$threat.code.matrix[i,] ~ 1,
                              weights = 1.0/(ranges$RangeArea^(1/3)), family = binomial(link = 'logit'))
        }
      }
      
      
      print("Saving models")
      if(!dir.exists(paste0(path,'Occupancy modelling/Logistic/',r))) dir.create(paste0(path,'Occupancy modelling/Logistic/',r))
      

      print("Predicting threat")
      p.r.cube.rt = rep(NA, nsite)
      
      for(i in seq_len(nsite))
      {
        if(i%%100 ==0)print(i)
        
        if(sum(!is.na(occ.matrices$threat.code.matrix[i,])) > 0)
        {
          p.r.cube.rt[i] = predict(r.models[[i]],type='response')[1]
          
        }
      }
      
      print("Saving predictions")
      save(p.r.cube.rt,file = paste0(path,'Occupancy modelling/Logistic/',r,'/cub.rt.R_W_Pred_',t,'_',names(th)))
      
    }
  }
}


