browns<-c(rgb(0.950,0.855,0.808),
          rgb(0.850,0.688,0.595),
          rgb(0.800,0.608,0.480),
          rgb(0.600,0.379,0.210),
          rgb(0.400,0.187,0.000),
          rgb(0.2,0.1,0))


DrawThreatCode<-function(u,u.c,tc)
{
  u.m<-1-u.c
  d.tc<-tc
  use.i<-which(!is.na(tc))
  DrawCode<-(runif(length(use.i)) < (u.c + (u.m*u)))
  d.tc[use.i[DrawCode]]<-1-tc[use.i[DrawCode]]
  d.tc
}



GetPredictionQuantiles<-function(path,model,model.name,t,th,realms, regions.comb,p)
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
    q.p.r<-get(model.name)#apply(p.r,1,function(x) max(x,na.rm=T))
    q.p.r[is.infinite(q.p.r)]<-NA
    #global.sr = c(global.sr,apply(p.r,1,function(x) sum(!is.na(x))))
    
    
    load(paste0(path,"Regional 2017_3/",r,"/",t,"_intersection.a.sp.R"))
    #load(paste0(path,"Regional 2017_3/Rdata/",r,"/",t,"_ranges.R"))
    
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
    
    q.p.r<-get(model.name)#apply(p.r,1,function(x) max(x,na.rm=T))
    q.p.r[is.infinite(q.p.r)]<-NA
    #global.sr = c(global.sr,apply(p.r,1,function(x) sum(!is.na(x))))
    
    
    load(paste0(path,"Regional 2017_3/",region,"/",t,"_intersection.a.sp.R"))
    #load(paste0(path,"Regional 2017_3/Rdata/",r,"/",t,"_ranges.R"))
    
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



PlotAxisInMultiPanel<-function(threats.to.use,n.breaks,lower.cut,range,axis.range,cols,sqrt = F)
{
  t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  par(mar=c(2,8,0,8))
  
  x.lim = axis.range
  y.lim = c(0,length(threats.to.use))
  
  width = 1/n.breaks
  x_step = seq(x.lim[1], x.lim[2], width)
  
  x.values = seq(x.lim[1],x.lim[2],length.out = 6)
  
  plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
  
  #rect(0, 0, 1, 1, col="light grey", border=NA)
  
  
  #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
  
  
  for(y in seq_len(length(threats.to.use)))
  {
    th = threats.to.use[y]
    th.cols <- colorRampPalette(cols[[th]])(n.breaks) 
    if(sqrt)
    {
      sapply(seq(1,n.breaks),
             function(i) {
               rect((t.breaks[i])^2, y-1, (t.breaks[i+1])^2, y, col=th.cols[i], border=NA)
               #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
             })
      
    }
    else
    {
      sapply(seq(1,n.breaks),
             function(i) {
               rect((t.breaks[i]), y-1, (t.breaks[i+1]), y, col=th.cols[i], border=NA)
               #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
             })
      
    }
  }
  
  axis(2,at = seq(0.5,length(threats.to.use)-0.5,1),
       labels = threats.to.use,las=1,cex.axis=0.7,tick = F)
  axis(1,padj=-2, cex.axis=0.7,
       at = pretty(x.values))
}


PlotSingleMapInMultiPanel<-function(d,cell.i,n.breaks,lower.cut,range,axis.range,color.scheme,sqrt = F,i,
                                    axis=T,axis.only = F,horiz=T)
{
  t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  t.labels = ''
  plot.norm = normalise_set_b(d[cell.i], t.breaks,t.labels)
  
  features.to.plot = which(plot.norm[[1]] > 0)

  
  if(!axis.only)
  {
    # par(oma=c(0,0,0,0))
     par(mar=c(0,0,0,0))
    # par(xpd=NA)
    # par(cex=1)
    
    plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
    #plot(bound,lwd=0.5,col = NA,border = "black", add=T)
    
    cols <- if(n.breaks == 2) colorRampPalette(color.scheme)(3) else colorRampPalette(color.scheme)(n.breaks) #if(color.scheme == "Browns") colorRampPalette(browns)(n.breaks) else (colorRampPalette((brewer.pal(8,color.scheme)))(n.breaks))
    palette(cols)
    plot(grid[cell.i[features.to.plot],],col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
    #plot(c.map,border ="black",add=T)
    
    mtext(text = letters[i],cex=1.2,side = 3,adj = 0.1,line = -2)
  }
  
  if(axis || axis.only)
  {
    if(horiz)
    {  
      par(mar=c(2,5,0,5))
      
      x.lim = axis.range
      y.lim = c(0,1)
      
      width = 1/n.breaks
      x_step = seq(x.lim[1], x.lim[2], width)
      
      x.values = seq(x.lim[1],x.lim[2],length.out = 6)
      
      x_steps<-t.breaks
      if(sqrt) x_steps<-x_steps^2
      
      plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
      sapply(seq(1,n.breaks),
             function(i) {
               rect((x_steps[i]), 0, (x_steps[i+1]), 1, col=i, border=NA)
               #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
             })
      axis(1,padj=-2, cex.axis=0.8,
           at = pretty(x.values))
      
    }
    else
    {
      par(mar=c(3,0,3,4.5))
      x.lim = c(0,1)
      y.lim = axis.range
      
      y.values = seq(y.lim[1],y.lim[2],length.out = 6)
      y_steps <- t.breaks
      if(sqrt) y_steps<-y_steps^2
      
      plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
      sapply(seq(1,n.breaks),
             function(i) {
               rect(0,(y_steps[i]), 1, (y_steps[i+1]), col=i, border=NA)
               #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
             })
      axis(4, cex.axis=1.0,
           at = pretty(y.values),las = 1)
    }
    
    
  }
  
}
  

PlotSingleMap<-function(d,cell.i,n.breaks,lower.cut,range,axis.range,fname,color.scheme,sqrt = F)
{
  t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  t.labels = ''
  plot.norm = normalise_set_b(d[cell.i], t.breaks,t.labels)
  
  features.to.plot = which(plot.norm[[1]] > 0)
  
  pdf(fname,width = 18/2.54,height = 15/2.54)
  #png(fname,width = 22,height = 11, units = 'cm',res = 600)
  
  par(oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  par(xpd=NA)
  par(cex=1)
  
  layout(matrix(c(1,2),nrow = 2), heights = c(0.85,0.15))
  
  plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
  #plot(bound,lwd=0.5,col = NA,border = "black", add=T)
  
  cols <- if(n.breaks == 2) colorRampPalette(color.scheme)(3) else colorRampPalette(color.scheme)(n.breaks) #if(color.scheme == "Browns") colorRampPalette(browns)(n.breaks) else (colorRampPalette((brewer.pal(8,color.scheme)))(n.breaks))
  palette(cols)
  plot(grid[cell.i[features.to.plot],],col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
  #plot(c.map,border ="black",add=T)
  
  par(mar=c(3,6,0,6))
  
  x.lim = axis.range
  y.lim = c(0,1)
  
  width = 1/n.breaks
  x_step = seq(x.lim[1], x.lim[2], width)
  
  x.values = seq(x.lim[1],x.lim[2],length.out = 6)
  
  plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
  
  #rect(0, 0, 1, 1, col="light grey", border=NA)
  
  
  #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
  
  
  if(sqrt)
  {
    sapply(seq(1,n.breaks),
           function(i) {
             rect((t.breaks[i])^2, 0, (t.breaks[i+1])^2, 1, col=i, border=NA)
             #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
           })
    
  }
   else
   {
     sapply(seq(1,n.breaks),
            function(i) {
              rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
              #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
            })
    
   }
  
  axis(1,padj=0, cex.axis=0.6,
       at = pretty(x.values))
  
  dev.off()
}


PlotCompositeMapInMultiPanel<-function(d,cell.i,n.breaks,lower.cut,range,axis.range,cols, t.breaks = NULL)
{
  if(is.null(t.breaks)) t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  
  # par(oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  # par(xpd=NA)
  # par(cex=1)
  
  plot(WL.moll,lwd=0.05,col = "lightgrey", border = "lightgrey")#add=T
  #plot(bound,lwd=0.5,col = NA,border = "black", add=T)
  
  for(th in 1:length(cols))
  {
    use.i<-which(d[cell.i] == th)
    
    if(is.null(t.breaks)) t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
    t.labels = ''
    th.d<-grid@data[[grep(paste0(t,"_",names(cols)[th]),colnames(grid@data))[1]]]
    plot.norm = normalise_set_b(th.d[cell.i][use.i], t.breaks,t.labels)
    features.to.plot = which(plot.norm[[1]] > 0)
    
    th.cols <- if(n.breaks == 2) colorRampPalette(rev(cols[[th]]))(4) else colorRampPalette(cols[[th]])(n.breaks)
    #th.cols <- colorRampPalette(cols[[th]])(n.breaks)# colorRampPalette(browns)(n.breaks) else (colorRampPalette((brewer.pal(8,cols[[th]])))(n.breaks))
    palette(th.cols)
    plot(grid[cell.i[use.i][features.to.plot],],
         col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
  }
  
}


PlotCompositeMap<-function(d,cell.i,n.breaks,lower.cut,range,axis.range,fname,cols, t.breaks = NULL)
{
  if(is.null(t.breaks)) t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  # t.labels = ''
  # plot.norm = normalise_set_b(d[cell.i], t.breaks,t.labels)
  # 
  # features.to.plot = which(plot.norm[[1]] > 0)
  # 
  png(fname,width = 17,height = 8, units = 'cm',res = 600)
  
  par(oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  par(xpd=NA)
  par(cex=1)
  
  #layout(matrix(c(1,2),nrow = 2), heights = c(0.9,0.1))
  
  plot(WL.moll,lwd=0.05,col = "lightgrey", border = "lightgrey")#add=T
  #plot(bound,lwd=0.5,col = NA,border = "black", add=T)
  
  for(th in 1:length(cols))
  {
    use.i<-which(d[cell.i] == th)
    
    if(is.null(t.breaks)) t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
    t.labels = ''
    th.d<-grid@data[[grep(paste0(t,"_",names(cols)[th]),colnames(grid@data))[1]]]
    plot.norm = normalise_set_b(th.d[cell.i][use.i], t.breaks,t.labels)
    features.to.plot = which(plot.norm[[1]] > 0)
    
    th.cols <- if(n.breaks == 2) colorRampPalette(rev(cols[[th]]))(4) else colorRampPalette(cols[[th]])(n.breaks)
    #th.cols <- colorRampPalette(cols[[th]])(n.breaks)# colorRampPalette(browns)(n.breaks) else (colorRampPalette((brewer.pal(8,cols[[th]])))(n.breaks))
    palette(th.cols)
    plot(grid[cell.i[use.i][features.to.plot],],
         col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
  }  
  # palette(cols)
  # 
  # plot(grid[cell.i,],col = d[cell.i], border = NA,add=T)
  #plot(c.map,border ="black",add=T)
  
  par(mar=c(2,3,0,3))
  
  x.lim = axis.range
  y.lim = c(0,1)
  
  width = 1/n.breaks
  x_step = seq(x.lim[1], x.lim[2], width)
  
  # plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
  # 
  # #rect(0, 0, 1, 1, col="light grey", border=NA)
  # 
  # sapply(seq(1,n.breaks),
  #        function(i) {
  #          rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
  #          #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
  #        })
  # #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
  # 
  # x.values = seq(x.lim[1],x.lim[2],length.out = 6)
  # axis(1,padj=-2, cex.axis=0.5,
  #      at = x.values)
  # 
  
  dev.off()
}

PlotCompositeMap2<-function(d,cell.i,n.breaks,lower.cut,range,axis.range,fname,cols,d2, cell.i.2 = NULL)
{
  t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  # t.labels = ''
  # plot.norm = normalise_set_b(d[cell.i], t.breaks,t.labels)
  # 
  # features.to.plot = which(plot.norm[[1]] > 0)
  # 
  png(fname,width = 17,height = 12, units = 'cm',res = 600)
  
  par(oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  par(xpd=NA)
  par(cex=1)
  
  layout(matrix(c(1,2),nrow = 2), heights = c(0.9,0.1))
  
  plot(WL.moll,lwd=0.05,col = "lightgrey", border = "lightgrey")#add=T
  #plot(bound,lwd=0.5,col = NA,border = "black", add=T)
  
  for(th in 1:length(cols))
  {
    use.i<-which(d[cell.i] == th)
    
    t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
    t.labels = ''
    th.d<-grid@data[[grep(paste0(t,"_",names(cols)[th]),colnames(grid@data))[1]]]
    plot.norm = normalise_set_b(th.d[cell.i][use.i], t.breaks,t.labels)
    features.to.plot = which(plot.norm[[1]] > 0)
    
    th.cols <- (colorRampPalette((brewer.pal(8,cols[[th]])))(n.breaks))
    palette(th.cols)
    plot(grid[cell.i[use.i][features.to.plot],],
         col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
  } 
  
  if(!is.null(cell.i.2))
  {
    for(th in 1:length(cols))
    {
      use.i<-which(d2[cell.i.2] == th)
      
      t.breaks = seq(range[1],range[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
      t.labels = ''
      th.d<-grid@data[[grep(paste0(t,"_",names(cols)[th]),colnames(grid@data))[1]]]
      plot.norm = normalise_set_b(th.d[cell.i.2][use.i], t.breaks,t.labels)
      features.to.plot = which(plot.norm[[1]] > 0)
      
      th.cols <- (colorRampPalette((brewer.pal(8,cols[[th]])))(n.breaks))
      palette(th.cols)
      if(length(use.i) > 0)
      {
        grid.centres<-rgeos::gCentroid(grid[cell.i.2[use.i][features.to.plot],],byid = T)
        
        plot(grid.centres,
             col = plot.norm[[1]][features.to.plot]+1, pch=16,add=T, cex=0.075)
      }
    }
  }
  
  
  # palette(cols)
  # 
  # plot(grid[cell.i,],col = d[cell.i], border = NA,add=T)
  #plot(c.map,border ="black",add=T)
  
  par(mar=c(2,3,0,3))
  
  x.lim = axis.range
  y.lim = c(0,1)
  
  width = 1/n.breaks
  x_step = seq(x.lim[1], x.lim[2], width)
  
  # plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
  # 
  # #rect(0, 0, 1, 1, col="light grey", border=NA)
  # 
  # sapply(seq(1,n.breaks),
  #        function(i) {
  #          rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
  #          #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
  #        })
  # #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
  # 
  # x.values = seq(x.lim[1],x.lim[2],length.out = 6)
  # axis(1,padj=-2, cex.axis=0.5,
  #      at = x.values)
  # 
  
  dev.off()
}

PlotBivariateMap<-function(d1,d2,cell.i,lower.cut,range1,range2,axis.range,fname,n.breaks)
{
  t.breaks.1 = seq(range1[1],range1[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  t.labels = ''
  plot.norm.1 = normalise_set_b(d1, t.breaks.1,t.labels)
  t.breaks.2 = seq(range2[1],range2[2],length.out = n.breaks)#c(0,seq(lower.cut,1,length.out = n.breaks))
  t.labels = ''
  plot.norm.2 = normalise_set_b(d2, t.breaks.2,t.labels)
  
  features.to.plot = which(plot.norm.1[[1]] > 0)
  
  png(fname,width = 17,height = 12, units = 'cm',res = 600)
  
  #par(oma=c(0,0,0,0))
  par(mar=c(0,0,0,0))
  #par(xpd=T)
  par(cex=1)
  
  layout(matrix(c(rep(1,20),rep(2,5)),nrow = 5,ncol=5,byrow = T), heights = 0.2)
  
  par(bg = "lightgrey")
  plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA,bg = "#636363")#add=T
  
  
  max.1 = max(plot.norm.1[[1]][features.to.plot])/(n.breaks-1)
  
  blue = 1/(n.breaks-1)
  cols = rgb(max.1,
             1-plot.norm.2[[1]][features.to.plot]/(n.breaks-1),
             (1-max.1),
             (plot.norm.1[[1]][features.to.plot]/(n.breaks-1)))
  
  
  
  
  plot(grid[cell.i[features.to.plot],],col = cols, border = NA,add=T)
  
  
  unique.1 = sort(unique(plot.norm.1[[1]][features.to.plot]),decreasing = F)
  unique.2 = sort(unique(plot.norm.2[[1]][features.to.plot]),decreasing = F)
  
  unique.1 = unique.1/(n.breaks-1)
  unique.2 = unique.2/(n.breaks-1)
  
  temp.mat<-matrix(seq_len(length(unique.1)*length(unique.2)),nrow=length(unique.1),ncol=length(unique.2))
  col.mat<-temp.mat
  for(i in seq_len(dim(temp.mat)[1]))
  {
    for(j in seq_len(dim(temp.mat)[2]))
    {
      col.mat[i,j] = rgb((max(unique.1)),
                         1-(unique.2[j]),
                         1- max(unique.1),
                         unique.1[i])
    }
  }
  
  
  
  par(mar=c(2,17,1,17))
  #par(mar=rep(0,4))
  #par('bg' = 'lightgrey')
  #plot(NA,xlim = c(0,1),ylim = c(0,1),axes = F,xlab = '',ylab = '')
  i=image(x = seq_len(dim(temp.mat)[1]),y= seq_len(dim(temp.mat)[2]),
          z=temp.mat,col=col.mat,xaxt="n",yaxt="n",xlab="",ylab="")
  axis(side = 1, at = seq_len(dim(temp.mat)[1]+1)-0.5,
       labels = t.breaks.1[seq_len(dim(temp.mat)[1]+1)],cex.axis=0.5, tcl = -0.1, line = 0,mgp=c(0,0,0.5))
  axis(side = 2, at = seq_len(dim(temp.mat)[2]+1)-0.5,
       labels = t.breaks.2[seq_len(dim(temp.mat)[2]+1)],cex.axis=0.5, tcl = -0.1,las = 1)
  #axis(side = 1, at = c(0,0.5,1), labels = labels.mrr, cex.axis=0.4, padj = -5, tcl = -0.1)
  mtext("Richness",side=2,cex=0.5, line = 2)
  mtext("Probability of threat occurence",side=1,cex=0.5, line = 1)
  
  
  dev.off()
}



GetThreatCodes<-function()
{
  threat.codes = list(
    'Hunting' = c(5.1),
    'Logging' = c(5.3),
    'Invasives'= c(8.1),
    'Fire'= c(7.1),
    'Human disturbance'  = c(6),
    'Pollution' = c(9),
    'Residential' = c(1),
    'Agriculture' = c(2),
    'Habitat conversion' = c(1,2,5.3),
    'Gathering plants' = c(5.2),
    'Fishing/Harvesting aquatic resources' = c(5.4), 
    'Climate change' = c(11)
  )
  threat.codes
}


GetThreatCodesPhasesText<-function()
{
  threat.codes.1 = list(
    'Habitat conversion' = c(1,2,5.3),
    'Residential' = c(1),
    'Agriculture' = c(2),
    'Logging' = c(5.3),
    'Hunting' = c(5.1),
    'Invasives'= c(8.1),
    'Pollution' = c(9), 
    'Climate change' = c(11),
    'AnyThreat' = c(1,2,5.3,5.1,8.1,9,11)
  )
  #Already run from Phase 1
  
  threat.codes.2 = list( 
    'Fire'= c(7.1),
    'Human disturbance'  = c(6)
  )
  
  
  
  list(phase1 = threat.codes.1,
       phase2 = threat.codes.2)
}



GetThreatCodesPhases<-function()
{
  threat.codes.1 = list(
    'Habitat conversion' = c(1,2,5.3),
    'Residential' = c(1),
    'Agriculture' = c(2),
    'Logging' = c(5.3),
    'Hunting' = c(5.1),
    'Invasives'= c(8.1),
    'Pollution' = c(9), 
    'Climate change' = c(11),
    'AnyThreat' = c(1,2,5.3,5.1,8.1,9,11)
  )
  #Already run from Phase 1
  
  threat.codes.2 = list( 
    'Fire'= c(7.1),
    'Human disturbance'  = c(6)
  )
  
  
  
  list(phase1 = threat.codes.1,
       phase2 = threat.codes.2)
}

CalculateSampledThreatCorrelation<-function(m.a.threat.lik.log, cells)
{
  #choose points that are on land
  #land.cells = which(grid@data$Realm == 1)
  cell.sample = sample(cells, size = length(cells)/2, replace = F)
  
  mod = list()
  
  for(th1 in seq_len(ncol(m.a.threat.lik.log)))
  {
    mod[[th1]] = list()
    for(th2 in seq(min(ncol(m.a.threat.lik.log),th1+1),ncol(m.a.threat.lik.log)))
    {
      mod[[th1]][[th2]] = lm(m.a.threat.lik.log[cell.sample,th1] ~ m.a.threat.lik.log[cell.sample,th2])
    }
  }
  
  r.sq = array(NA, dim=c(ncol(m.a.threat.lik.log),ncol(m.a.threat.lik.log)))
  for(th1 in seq_len(ncol(m.a.threat.lik.log)))
  {
    for(th2 in seq(min(ncol(m.a.threat.lik.log),th1+1),ncol(m.a.threat.lik.log)))
    {
      r.sq[th1,th2] = summary(mod[[th1]][[th2]])$r.squared
    }
  }
  r.sq
}




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


ByCellThreatenedSpeciesBinomialsDraws<-function(
  threats,
  range.sizes,draw)
{
  r.threats = threats[(threats$binomial %in% range.sizes$binomial),]

  th.sp.binomials <- r.threats$binomial[which(r.threats[[paste0("Draw.",d)]] == 1)]
  th.sp.binomials
}

Threatened.Unthreatened.SpeciesBinomials<-function(threat.species,threats,threat.codes)
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
    taxa.th.sp.ids[[th]] = th.sp.ids[(th.sp.ids %in% threat.species$Species.ID)]
  }
  
  th.sp.binomials = list()
  un.th.sp.binomials = list()
  for(th in names(threat.codes))
  {
    th.sp.binomials[[th]] = threat.species$binomial[threat.species$Species.ID %in% taxa.th.sp.ids[[th]]]
    un.th.sp.binomials[[th]] = threat.species$binomial[!(threat.species$Species.ID %in% taxa.th.sp.ids[[th]])]
  }
  
  
  
  list(threatened = th.sp.binomials,
       unthreathened= un.th.sp.binomials)
}


ByCellThreatenedSpeciesBinomials<-function(threat.species.info, range.sizes,threats,threat.codes)
{
  threat.species = threat.species.info[(threat.species.info$binomial %in% range.sizes$binomial),]
  
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
    taxa.th.sp.ids[[th]] = th.sp.ids[(th.sp.ids %in% threat.species$Species.ID)]
  }
  
  th.sp.binomials = list()
  for(th in names(threat.codes))
  {
    th.sp.binomials[[th]] = threat.species$binomial[threat.species$Species.ID %in% taxa.th.sp.ids[[th]]]
  }
  
  th.sp.binomials
}

ByCellThreatenedBirdBinomials<-function(ss, ranges,threat.codes)
{
  ss.i<-match(as.character(ranges$binomil),ss$Scientific.name)
  threat.sp.binomials = list()
  
  for(th in names(threat.codes))
  {
    th.sp.binomials = list()
    for(code in threat.codes[[th]])
    {
      if(length((grep("\\.", code))) > 0)
      {
        code.rows = grep(paste0("^",code,"{1}"),ss$ThreatCode)
      } else
      {
        code.rows = grep(paste0("^",code,"{1}\\."),ss$ThreatCode)  
      }
      
      if(length(code.rows) > 0) th.sp.binomials = c(th.sp.binomials,ss$Scientific.name[code.rows])
      
    }
    th.sp.binomials = unique(unlist(th.sp.binomials))
    threat.sp.binomials[[th]] = th.sp.binomials
  }
  threat.sp.binomials
}

ThreatLikelihood<-function(a.array, threat.codes,range.sizes,th.sp.binomials,non.zero.sp.i,sp.i)
{
  
  threat.lik = array(NA,dim = c(nrow(a.array),length(threat.codes)))
  
  for(i in seq_len(length(non.zero.sp.i)))
  {
    c = non.zero.sp.i[i] 
    if(i %% floor(length(non.zero.sp.i)/100) == 0) print(paste(i,"out of",length(non.zero.sp.i)))
    
    #For this cell construct a table with the range area and area of intersection for each species intersecting
    intersection.areas = a.array[c,][sp.i[[c]]]
    df.sp.cell = data.frame(range.sizes[sp.i[[c]],],intersection.areas)
    colnames(df.sp.cell)[3] = "Intersection area"
    df.sp.cell$Inters.over.range = df.sp.cell$`Intersection area`/df.sp.cell$RangeArea
    
    # large.ranged = which(df.sp.cell$RangeArea > 250E3)
    # if(length(large.ranged) > 0)  df.sp.cell = df.sp.cell[-large.ranged,]
    # 
    if(nrow(df.sp.cell) > 0)
    {
      for(th in seq_len(length(threat.codes)))
      {
        threatened = which(df.sp.cell$binomial %in% th.sp.binomials[[th]])
        if(length(threatened) > 0)
        {
          numerator = sum((df.sp.cell$Inters.over.range[threatened]),na.rm=T)
          #denominator = numerator + nrow(df.sp.cell)-length(threatened)
          denominator = nrow(df.sp.cell)
          threat.lik[c,th] = (numerator/denominator)
          #mammal.threat.lik[c,th] = numerator/denominator
        }
        else
        {
          threat.lik[c,th] = 0
        }
      }
    }
  }
  
  threat.lik
}



UpperBoundededLogProfile<-function(x)
{
  i.log = 1/log(x)
  if(x < exp(1))
  {
    1
  }
  else
  {
    i.log
  }
  
}


ThreatLikelihoodLog<-function(a.array, threat.codes,range.sizes,th.sp.binomials,non.zero.sp.i,sp.i)
{
  
  threat.lik = array(NA,dim = c(nrow(a.array),length(threat.codes)))
  
  for(i in seq_len(length(non.zero.sp.i)))
  {
    c = non.zero.sp.i[i] 
    if(i %% floor(length(non.zero.sp.i)/100) == 0) print(paste(i,"out of",length(non.zero.sp.i)))
    
    #For this cell construct a table with the range area and area of intersection for each species intersecting
    intersection.areas = a.array[c,][sp.i[[c]]]
    df.sp.cell = data.frame(range.sizes[sp.i[[c]],],intersection.areas)
    colnames(df.sp.cell)[3] = "Intersection area"
    df.sp.cell$Inters.over.range = df.sp.cell$`Intersection area`/df.sp.cell$RangeArea
    df.sp.cell$Range.over.inters = df.sp.cell$RangeArea/df.sp.cell$`Intersection area`
    df.sp.cell$Inv.log.Range.over.inters = sapply(df.sp.cell$Range.over.inters,UpperBoundededLogProfile)
    
    # large.ranged = which(df.sp.cell$RangeArea > 250E3)
    # if(length(large.ranged) > 0)  df.sp.cell = df.sp.cell[-large.ranged,]
    # 
    if(nrow(df.sp.cell) > 0)
    {
      for(th in seq_len(length(threat.codes)))
      {
        threatened = which(df.sp.cell$binomial %in% th.sp.binomials[[th]])
        if(length(threatened) > 0)
        {
          #numerator = sum(1/(log(df.sp.cell$Range.over.inters[threatened])),na.rm=T)
          numerator = sum(df.sp.cell$Inv.log.Range.over.inters[threatened],na.rm = T)
          #denominator = numerator + nrow(df.sp.cell)-length(threatened)
          denominator = nrow(df.sp.cell)
          threat.lik[c,th] = (numerator/denominator)
          #mammal.threat.lik[c,th] = numerator/denominator
        }
        else
        {
          threat.lik[c,th] = 0
        }
      }
    }
  }
  
  threat.lik
}


ThreatLikelihoodMaxRange<-function(a.array, threat.codes,range.sizes,th.sp.binomials,non.zero.sp.i,sp.i,max.range)
{
  
  threat.lik = array(NA,dim = c(nrow(a.array),length(threat.codes)))
  
  for(i in seq_len(length(non.zero.sp.i)))
  {
    c = non.zero.sp.i[i] 
    if(i %% floor(length(non.zero.sp.i)/100) == 0) print(paste(i,"out of",length(non.zero.sp.i)))
    
    #For this cell construct a table with the range area and area of intersection for each species intersecting
    intersection.areas = a.array[c,][sp.i[[c]]]
    df.sp.cell = data.frame(range.sizes[sp.i[[c]],],intersection.areas)
    colnames(df.sp.cell)[3] = "Intersection area"
    df.sp.cell$Inters.over.range = df.sp.cell$`Intersection area`/df.sp.cell$RangeArea
    
    
    #Exclude any species whose ranges are greater than the supplied max.range
    large.ranged = which(df.sp.cell$RangeArea > max.range)
    if(length(large.ranged) > 0)  df.sp.cell = df.sp.cell[-large.ranged,]
     
    if(nrow(df.sp.cell) > 0)
    {
      for(th in seq_len(length(threat.codes)))
      {
        threatened = which(df.sp.cell$binomial %in% th.sp.binomials[[th]])
        if(length(threatened) > 0)
        {
          numerator = sum((1/df.sp.cell$RangeArea[threatened]))
          #denominator = numerator + nrow(df.sp.cell)-length(threatened)
          denominator = nrow(df.sp.cell)
          threat.lik[c,th] = (numerator/denominator)
          #mammal.threat.lik[c,th] = numerator/denominator
        }
        else
        {
          threat.lik[c,th] = 0
        }
      }
    }
  }
  
  threat.lik
}



ThreatLikelihoodOmni<-function(a.array, threat.codes,range.sizes,th.sp.binomials,non.zero.sp.i,sp.i)
{
  
  threat.lik = array(NA,dim = c(nrow(a.array),length(threat.codes)))
  
  for(i in seq_len(length(non.zero.sp.i)))
  {
    c = non.zero.sp.i[i] 
    if(i %% floor(length(non.zero.sp.i)/100) == 0) print(paste(i,"out of",length(non.zero.sp.i)))
    
    #For this cell construct a table with the range area and area of intersection for each species intersecting
    intersection.areas = a.array[c,][sp.i[[c]]]
    df.sp.cell = data.frame(range.sizes[sp.i[[c]],],intersection.areas)
    colnames(df.sp.cell)[3] = "Intersection area"
    df.sp.cell$Inters.over.range = df.sp.cell$`Intersection area`/df.sp.cell$RangeArea
    
    # large.ranged = which(df.sp.cell$RangeArea > 250E3)
    # if(length(large.ranged) > 0)  df.sp.cell = df.sp.cell[-large.ranged,]
    # 
    if(nrow(df.sp.cell) > 0)
    {
      for(th in seq_len(length(threat.codes)))
      {
        threatened = which(df.sp.cell$binomial %in% th.sp.binomials[[th]])
        if(length(threatened) > 0)
        {
          numerator = length(threatened)
          #denominator = numerator + nrow(df.sp.cell)-length(threatened)
          denominator = nrow(df.sp.cell)
          #threat.lik[c,th] = exp(mean(log(df.sp.cell$Inters.over.range[threatened]),na.rm=T))
          threat.lik[c,th] = numerator/denominator
        }
        else
        {
          threat.lik[c,th] = 0
        }
      }
    }
  }
  
  threat.lik
}




PlotThreatLikelihoodMapsLog<-function(threat.lik,out.dir, n.breaks,threat.codes,suffix)
{
  
  #Assign any values greater than 1 to 1 to avoid spurious rounding issues
  threat.lik[threat.lik > 1] = 1
  
  
  
  set.seed(2016)
  
  for(th in c(1))#seq_len(length(threat.codes)))
  {
    
    
    height.cm = 10
    width.cm = 14
    #png(paste(outdir,"WellsPipes_moll_bg_eez_grey.png",sep=""),width=width.cm,height=height.cm,units="cm",res=600)
    pdf(file = paste(out.dir,names(threat.codes)[th],"_g_mean",suffix,".pdf",sep=""),width=width.cm/2.54,height=height.cm/2.54)
    
    t = log10(threat.lik[,th])
    #t.lims = range(t,na.rm = T)
    #t.breaks = (seq(t.lims[1],t.lims[2],length.out = n.breaks))
    #t.breaks.jenks = getJenksBreaks(t, n.breaks+1)
    #t.breaks.kmeans = classIntervals(t,n = n.breaks-1, style = "kmeans",rtimes = 100)
    #t.breaks.quants = classIntervals(t,n = n.breaks-1, style = "quantile")
    #t.breaks.hclust = classIntervals(t,n = n.breaks-1, style = "bclust")
    #t.breaks = t.breaks.kmeans$brks
    lower.cut = -3
    t.breaks = c(-13,seq(lower.cut,0,length.out = n.breaks))
    t.labels = ''
    t.norm = normalise_set_b(t, t.breaks,t.labels)
    
    par(oma=c(0,0,0,0))
    par(mar=c(0,0,0,0))
    par(xpd=NA)
    par(cex=1)
    
    layout(matrix(c(1,2),nrow = 2), heights = c(0.85,0.15))
    
    plot.norm = t.norm
    
    features.to.plot = which(plot.norm[[1]] > 0)
    #features.to.plot = which(sr > 1100)
    #plot(eez,border = "black", add = T)
    plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
    plot(bound,lwd=0.5,col = NA,border = "black", add=T)
    cols = col2rgb(rev(brewer.pal(length(plot.norm[[2]]),"RdYlBu")))
    #Make the colours darker
    cols = round(cols*0.8)
    cols = apply(cols,2,function(x) paste(as.hexmode(x),collapse =""))
    cols = paste("#",cols,sep="")
    
    cols <- (colorRampPalette(rev(brewer.pal(8,"RdYlBu")))(n.breaks))
    
    palette(cols)
    plot(grid[features.to.plot,],col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
    #plot(grid[plot.norm[[1]] == 0,],col = "#B04573", border = NA,add=T)
    #plot(WL.moll, add=T, col = NA, border= "black",lwd=0.05)
    
    par(mar=c(2,3,0,3))
    
    # x.lim = range(t.breaks)
    
    
    # width <- (x.lim[2]-x.lim[1])/(n.breaks-1)
    # x_step <- seq(x.lim[1], x.lim[2], width)
    
    x.lim = c(t.breaks[2]-(2*(t.breaks[3]-t.breaks[2])),0)
    y.lim = c(0,1)
    
    width = 1/n.breaks
    x_step = seq(0, 1, width)
    
    plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
    
    sapply(seq(2,n.breaks),
           function(i) {
             rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
             #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
           })
    #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
    
    x.values = seq(lower.cut,0)
    # x.axis.pos= c(width,
    #               (1:length(x.values))*(1 - (2*width))/length(x.values),
    #               1-width)#seq(x.lim[1],x.lim[2],length.out = 11)
    #x.labs = c(pretty10exp(10^(head(test,length(test)-1)+1),drop.1 = T),1)
    #x.labs = c(pretty10exp(10^x.values,digits = 1,drop.1 = T))
    x.labs = format(10^seq(lower.cut,0),scientific = F,drop0trailing = T)
    x.labs[1] = paste0('<',x.labs[1])
    #pretty10exp(c(0,t.breaks[seq(2,n.breaks-1,length.out = 9)],1),digits = 1,sub10 = 2)
    #x.labs = format(c(0,t.breaks[seq(2,n.breaks-1,length.out = 9)],1),digits = 1)
    
    #axis(1,at=x.axis.pos,labels = exp(x.axis.pos), tick = F, line=-1,cex.axis=0.7)
    #Complete the x axis
    #eaxis(1,padj=-1, cex.axis=0.8)
    axis(1,padj=-2, cex.axis=0.5,
         at = x.values,
         labels = x.labs)
    
    mtext(side = 1, names(threat.codes[th]), line = 1, cex = 0.8)
    dev.off()
  }
  
}


PlotThreatLikelihoodMapsLn<-function(threat.lik,out.dir, n.breaks,threat.codes,suffix)
{
  
  #Assign any values greater than 1 to 1 to avoid spurious rounding issues
  threat.lik[threat.lik > 1] = 1
  
  
  
  set.seed(2016)
  
  for(th in c(1))#seq_len(length(threat.codes)))
  {
    
    
    height.cm = 10
    width.cm = 14
    #png(paste(outdir,"WellsPipes_moll_bg_eez_grey.png",sep=""),width=width.cm,height=height.cm,units="cm",res=600)
    pdf(file = paste(out.dir,names(threat.codes)[th],"_g_mean",suffix,".pdf",sep=""),width=width.cm/2.54,height=height.cm/2.54)
    
    t = log(threat.lik[,th])
    #t.lims = range(t,na.rm = T)
    #t.breaks = (seq(t.lims[1],t.lims[2],length.out = n.breaks))
    #t.breaks.jenks = getJenksBreaks(t, n.breaks+1)
    #t.breaks.kmeans = classIntervals(t,n = n.breaks-1, style = "kmeans",rtimes = 100)
    #t.breaks.quants = classIntervals(t,n = n.breaks-1, style = "quantile")
    #t.breaks.hclust = classIntervals(t,n = n.breaks-1, style = "bclust")
    #t.breaks = t.breaks.kmeans$brks
    lower.cut = -4
    t.breaks = c(-13,seq(lower.cut,0,length.out = n.breaks))
    t.labels = ''
    t.norm = normalise_set_b(t, t.breaks,t.labels)
    
    par(oma=c(0,0,0,0))
    par(mar=c(0,0,0,0))
    par(xpd=NA)
    par(cex=1)
    
    layout(matrix(c(1,2),nrow = 2), heights = c(0.85,0.15))
    
    plot.norm = t.norm
    
    features.to.plot = which(plot.norm[[1]] > 0)
    #features.to.plot = which(sr > 1100)
    #plot(eez,border = "black", add = T)
    plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
    plot(bound,lwd=0.5,col = NA,border = "black", add=T)
    cols = col2rgb(rev(brewer.pal(length(plot.norm[[2]]),"RdYlBu")))
    #Make the colours darker
    cols = round(cols*0.8)
    cols = apply(cols,2,function(x) paste(as.hexmode(x),collapse =""))
    cols = paste("#",cols,sep="")
    
    cols <- (colorRampPalette(rev(brewer.pal(8,"RdYlBu")))(n.breaks))
    
    palette(cols)
    plot(grid[features.to.plot,],col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
    #plot(grid[plot.norm[[1]] == 0,],col = "#B04573", border = NA,add=T)
    #plot(WL.moll, add=T, col = NA, border= "black",lwd=0.05)
    
    par(mar=c(2,3,0,3))
    
    # x.lim = range(t.breaks)
    
    
    # width <- (x.lim[2]-x.lim[1])/(n.breaks-1)
    # x_step <- seq(x.lim[1], x.lim[2], width)
    
    x.lim = c(0,1)
    y.lim = c(0,1)
    
    width = 1/n.breaks
    x_step = seq(0, 1, width)
    
    plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
    
    rect(0, 0, 1, 1, col="light grey", border=NA)
    
    sapply(seq(2,n.breaks),
           function(i) {
             rect(exp(t.breaks[i]), 0, exp(t.breaks[i+1]), 1, col=i, border=NA)
             #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
           })
    #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
    
    x.values = seq(x.lim[1],x.lim[2],0.2)
    # x.axis.pos= c(width,
    #               (1:length(x.values))*(1 - (2*width))/length(x.values),
    #               1-width)#seq(x.lim[1],x.lim[2],length.out = 11)
    #x.labs = c(pretty10exp(10^(head(test,length(test)-1)+1),drop.1 = T),1)
    #x.labs = c(pretty10exp(10^x.values,digits = 1,drop.1 = T))
    x.labs = format(10^seq(lower.cut,0),scientific = F,drop0trailing = T)
    x.labs[1] = paste0('<',x.labs[1])
    #pretty10exp(c(0,t.breaks[seq(2,n.breaks-1,length.out = 9)],1),digits = 1,sub10 = 2)
    #x.labs = format(c(0,t.breaks[seq(2,n.breaks-1,length.out = 9)],1),digits = 1)
    
    #axis(1,at=x.axis.pos,labels = exp(x.axis.pos), tick = F, line=-1,cex.axis=0.7)
    #Complete the x axis
    #eaxis(1,padj=-1, cex.axis=0.8)
    axis(1,padj=-2, cex.axis=0.5,
         at = x.values,
         labels = x.values)
    
    mtext(side = 1, names(threat.codes[th]), line = 1, cex = 0.8)
    dev.off()
  }
  
}



PlotThreatLikelihoodMapsLinear<-function(threat.lik,out.dir, n.breaks,threat.codes,suffix)
{
  
  #Assign any values greater than 1 to 1 to avoid spurious rounding issues
  threat.lik[threat.lik > 1] = 1
  
  
  
  set.seed(2016)
  
  for(th in c(1))#seq_len(length(threat.codes)))
  {
    
    
    height.cm = 10
    width.cm = 14
    #png(paste(outdir,"WellsPipes_moll_bg_eez_grey.png",sep=""),width=width.cm,height=height.cm,units="cm",res=600)
    pdf(file = paste(out.dir,names(threat.codes)[th],"_g_mean",suffix,".pdf",sep=""),width=width.cm/2.54,height=height.cm/2.54)
    
    t = (threat.lik[,th])
    t.lims = range(t,na.rm = T)
    t.breaks = (seq(t.lims[1],t.lims[2],length.out = n.breaks))
    #t.breaks.jenks = getJenksBreaks(t, n.breaks+1)
    #t.breaks.kmeans = classIntervals(t,n = n.breaks-1, style = "kmeans",rtimes = 100)
    #t.breaks.quants = classIntervals(t,n = n.breaks-1, style = "quantile")
    #t.breaks.hclust = classIntervals(t,n = n.breaks-1, style = "bclust")
    #t.breaks = t.breaks.kmeans$brks
    # lower.cut = -3
    # t.breaks = c(-13,seq(lower.cut,0,length.out = n.breaks))
    t.labels = ''
    t.norm = normalise_set_b(t, t.breaks,t.labels)
    
    par(oma=c(0,0,0,0))
    par(mar=c(0,0,0,0))
    par(xpd=NA)
    par(cex=1)
    
    layout(matrix(c(1,2),nrow = 2), heights = c(0.85,0.15))
    
    plot.norm = t.norm
    
    features.to.plot = which(plot.norm[[1]] > 0)
    #features.to.plot = which(sr > 1100)
    #plot(eez,border = "black", add = T)
    plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
    plot(bound,lwd=0.5,col = NA,border = "black", add=T)
    cols = col2rgb(rev(brewer.pal(length(plot.norm[[2]]),"RdYlBu")))
    #Make the colours darker
    cols = round(cols*0.8)
    cols = apply(cols,2,function(x) paste(as.hexmode(x),collapse =""))
    cols = paste("#",cols,sep="")
    
    cols <- (colorRampPalette(rev(brewer.pal(8,"RdYlBu")))(n.breaks))
    
    palette(cols)
    plot(grid[features.to.plot,],col = plot.norm[[1]][features.to.plot]+1, border = NA,add=T)
    #plot(grid[plot.norm[[1]] == 0,],col = "#B04573", border = NA,add=T)
    #plot(WL.moll, add=T, col = NA, border= "black",lwd=0.05)
    
    par(mar=c(2,3,0,3))
    
    # x.lim = range(t.breaks)
    
    
    # width <- (x.lim[2]-x.lim[1])/(n.breaks-1)
    # x_step <- seq(x.lim[1], x.lim[2], width)
    
    x.lim = c(0,1)
    y.lim = c(0,1)
    
    width = 1/n.breaks
    x_step = seq(0, 1, width)
    
    plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
    
    rect(0, 0, 1, 1, col="light grey", border=NA)
    
    sapply(seq(2,n.breaks),
           function(i) {
             rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
             #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
           })
    #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
    
    x.values = seq(x.lim[1],x.lim[2],0.2)
    # x.axis.pos= c(width,
    #               (1:length(x.values))*(1 - (2*width))/length(x.values),
    #               1-width)#seq(x.lim[1],x.lim[2],length.out = 11)
    #x.labs = c(pretty10exp(10^(head(test,length(test)-1)+1),drop.1 = T),1)
    #x.labs = c(pretty10exp(10^x.values,digits = 1,drop.1 = T))
    # x.labs = format(10^seq(lower.cut,0),scientific = F,drop0trailing = T)
    # x.labs[1] = paste0('<',x.labs[1])
    #pretty10exp(c(0,t.breaks[seq(2,n.breaks-1,length.out = 9)],1),digits = 1,sub10 = 2)
    #x.labs = format(c(0,t.breaks[seq(2,n.breaks-1,length.out = 9)],1),digits = 1)
    
    #axis(1,at=x.axis.pos,labels = exp(x.axis.pos), tick = F, line=-1,cex.axis=0.7)
    #Complete the x axis
    #eaxis(1,padj=-1, cex.axis=0.8)
    axis(1,padj=-2, cex.axis=0.5,
         at = x.values)
    
    mtext(side = 1, names(threat.codes[th]), line = 1, cex = 0.8)
    dev.off()
  }
  
}

PlotThreatLikelihoodMapsJenks<-function(threat.lik,out.dir, n.breaks,threat.codes,suffix)
{
  
  #Assign any values greater than 1 to 1 to avoid spurious rounding issues
  threat.lik[threat.lik > 1] = 1
  
  
  
  set.seed(2016)
  
  for(th in seq_len(length(threat.codes)))
  {
    
    
    height.cm = 10
    width.cm = 14
    pdf(file = paste(out.dir,names(threat.codes)[th],"_g_mean",suffix,".pdf",sep=""),width=width.cm/2.54,height=height.cm/2.54)
    
    t = c(threat.lik[,th])
    #t.lims = range(t,na.rm = T)
    #t.breaks = (seq(t.lims[1],t.lims[2],length.out = n.breaks))
    t.breaks.jenks = getJenksBreaks(t, n.breaks+1)
    #t.breaks.kmeans = classIntervals(t,n = n.breaks-1, style = "kmeans",rtimes = 100)
    #t.breaks.quants = classIntervals(t,n = n.breaks-1, style = "quantile")
    #t.breaks.hclust = classIntervals(t,n = n.breaks-1, style = "bclust")
    #t.breaks = t.breaks.kmeans$brks
    lower.cut = -4
    t.breaks = unique(t.breaks.jenks)
    t.breaks[length(t.breaks)] = 1.0
    #t.breaks = c(-13,seq(lower.cut,0,length.out = n.breaks))
    t.labels = ''
    t.norm = normalise_set_b(t, t.breaks,t.labels)
    
    par(oma=c(0,0,0,0))
    par(mar=c(0,0,0,0))
    par(xpd=NA)
    par(cex=1)
    
    layout(matrix(c(1,2),nrow = 2), heights = c(0.85,0.15))
    
    plot.norm = t.norm
    
    features.to.plot = which(plot.norm[[1]] >= 0)
    #features.to.plot = which(sr > 1100)
    #plot(eez,border = "black", add = T)
    plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
    plot(bound,lwd=0.5,col = NA,border = "black", add=T)
    # cols = col2rgb(rev(brewer.pal(length(plot.norm[[2]]),"RdYlBu")))
    # #Make the colours darker
    # cols = round(cols*0.8)
    # cols = apply(cols,2,function(x) paste(as.hexmode(x),collapse =""))
    # cols = paste("#",cols,sep="")
    # 
    cols <- (colorRampPalette(rev(brewer.pal(8,"RdYlBu")))(length(t.breaks)-1))
    
    palette(cols)
    plot(grid[features.to.plot,],col = plot.norm[[1]][features.to.plot], border = NA,add=T)
    #plot(grid[plot.norm[[1]] == 0,],col = "#B04573", border = NA,add=T)
    #plot(WL.moll, add=T, col = NA, border= "black",lwd=0.05)
    
    par(mar=c(2,3,0,3))
    
    # x.lim = range(t.breaks)
    
    
    # width <- (x.lim[2]-x.lim[1])/(n.breaks-1)
    # x_step <- seq(x.lim[1], x.lim[2], width)
    
    x.lim = c(0,1)
    y.lim = c(0,1)
    
    width = 1/n.breaks
    x_step = seq(0, 1, width)
    
    plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
    
    rect(0, 0, 1, 1, col="light grey", border=NA)
    
    sapply(seq(1,length(t.breaks)-1),
           function(i) {
             rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
             #rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
           })
    #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
    
    x.values = seq(x.lim[1],x.lim[2],0.2)
    # x.axis.pos= c(width,
    #               (1:length(x.values))*(1 - (2*width))/length(x.values),
    #               1-width)#seq(x.lim[1],x.lim[2],length.out = 11)
    #x.labs = c(pretty10exp(10^(head(test,length(test)-1)+1),drop.1 = T),1)
    #x.labs = c(pretty10exp(10^x.values,digits = 1,drop.1 = T))
    x.labs = format(10^seq(lower.cut,0),scientific = F,drop0trailing = T)
    x.labs[1] = paste0('<',x.labs[1])
    #pretty10exp(c(0,t.breaks[seq(2,n.breaks-1,length.out = 9)],1),digits = 1,sub10 = 2)
    #x.labs = format(c(0,t.breaks[seq(2,n.breaks-1,length.out = 9)],1),digits = 1)
    
    #axis(1,at=x.axis.pos,labels = exp(x.axis.pos), tick = F, line=-1,cex.axis=0.7)
    #Complete the x axis
    #eaxis(1,padj=-1, cex.axis=0.8)
    axis(1,padj=-2, cex.axis=0.5,
         at = x.values,
         labels = x.values)
    
    mtext(side = 1, names(threat.codes[th]), line = 1, cex = 0.8)
    dev.off()
  }
  
}


FindCountryPolygons<-function(c,c.map,type)
{
  c.poly = which(c.map@data$ISO3 %in% c & c.map@data$type == type)
  c.poly
}

SummariseThreatsCountries<-function(c,c.map,type,g,probs,threat.lik)
{
  
  c.poly = c.map[FindCountryPolygons(c,c.map,type),]
  
  
  #Perform logical intersection to see which cells fall within this polygon
  c.intersects = gIntersects(c.poly,g,byid = T)
  
  grid.sub2 <- apply(c.intersects, 1, function(x) {sum(x)}) # test across all polygons in the SpatialPolygon whether it intersects or not
  grid.intersects = which(grid.sub2 > 0)
  
  
  
  c.th.q = apply(threat.lik[grid.intersects,],2,function(x) quantile(x,probs = probs,na.rm=T))
  
  list(threat.summary = c.th.q,
       grid.intersects = grid.intersects)
}


PlotThreatLikelihoodMapsTop<-function(threat.lik,out.dir, q,threat.codes,suffix)
{
  
  #Assign any values greater than 1 to 1 to avoid spurious rounding issues
  threat.lik[threat.lik > 1] = 1
  
  
  
  set.seed(2016)
  
  for(th in seq_len(length(threat.codes)))
  {
    
    
    height.cm = 10
    width.cm = 14
    pdf(file = paste(out.dir,names(threat.codes)[th],"_g_mean_hotspots",suffix,".pdf",sep=""),width=width.cm/2.54,height=height.cm/2.54)
    
    t = c(threat.lik[,th])
    
    t.ind = rep(0,length(t))
    t.ind[which(t >= q)] = 1
    
    par(oma=c(0,0,0,0))
    par(mar=c(0,0,0,0))
    par(xpd=NA)
    par(cex=1)
    
    layout(matrix(c(1,2),nrow = 2), heights = c(0.85,0.15))

    features.to.plot = which(t.ind > 0)
    
    plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
    plot(bound,lwd=0.5,col = NA,border = "black", add=T)
    
    plot(grid[features.to.plot,],col = "red", border = NA,add=T)


    dev.off()
  }
  
}

PlotThreatLikelihoodMapsJenksEasy<-function(threat.lik,out.dir, t.breaks,threat.codes,suffix)
{
  
  #Assign any values greater than 1 to 1 to avoid spurious rounding issues
  threat.lik[threat.lik > 1] = 1
  
  for(th in seq_len(length(threat.codes)))
  {
    
    
    height.cm = 10
    width.cm = 14
    pdf(file = paste(out.dir,names(threat.codes)[th],"_g_mean_easy",suffix,".pdf",sep=""),width=width.cm/2.54,height=height.cm/2.54)
    
    t = c(threat.lik[,th])
    #t.breaks = c(-13,seq(lower.cut,0,length.out = n.breaks))
    t.labels = ''
    t.norm = normalise_set_b(t, t.breaks,t.labels)
    
    par(oma=c(0,0,0,0))
    par(mar=c(0,0,0,0))
    par(xpd=NA)
    par(cex=1)
    
    layout(matrix(c(1,2),nrow = 2), heights = c(0.85,0.15))
    
    plot.norm = t.norm
    
    features.to.plot = which(plot.norm[[1]] >= 0)
    #features.to.plot = which(sr > 1100)
    #plot(eez,border = "black", add = T)
    plot(WL.moll,lwd=0.05,col = "lightgrey", border = NA)#add=T
    plot(bound,lwd=0.5,col = NA,border = "black", add=T)
    # cols = col2rgb(rev(brewer.pal(length(plot.norm[[2]]),"RdYlBu")))
    # #Make the colours darker
    # cols = round(cols*0.8)
    # cols = apply(cols,2,function(x) paste(as.hexmode(x),collapse =""))
    # cols = paste("#",cols,sep="")
    # 
    cols <- (colorRampPalette(rev(brewer.pal(8,"RdYlBu")))(length(t.breaks)-1))
    
    palette(cols)
    plot(grid[features.to.plot,],col = plot.norm[[1]][features.to.plot], border = NA,add=T)
    #plot(grid[plot.norm[[1]] == 0,],col = "#B04573", border = NA,add=T)
    #plot(WL.moll, add=T, col = NA, border= "black",lwd=0.05)
    
    par(mar=c(2.5,3,0,3))
    
    # x.lim = range(t.breaks)
    
    
    # width <- (x.lim[2]-x.lim[1])/(n.breaks-1)
    # x_step <- seq(x.lim[1], x.lim[2], width)
    
    x.lim = c(0,1)
    y.lim = c(0,1)
    
    width = 1/n.breaks
    x_step = seq(0, 1, width)
    
    plot(NA, axes=F, xlab="",ylab="",xlim = x.lim, ylim = y.lim)
    
    rect(0, 0, 1, 1, col="light grey", border=NA)
    
    sapply(seq(1,length(t.breaks)-1),
           function(i) {
             #rect((t.breaks[i]), 0, (t.breaks[i+1]), 1, col=i, border=NA)
             rect(x_step[i], 0, x_step[i+1], 1, col=i, border=NA)
           })
    #rect(x.lim[1], 0, t.breaks[2], 1, col="#B04573", border=NA)
    
    x.values = seq(0.1,0.9,0.2)
    
    axis(1,padj=-2, cex.axis=0.5,
         at = x.values,
         labels = c("Very Low","Low","Medium","High","Very High"),tick = F)
    
    #mtext(side = 1, names(threat.codes[th]), line = 1, cex = 0.8)
    dev.off()
  }
  
}

ConstructOccupancyMatrices<-function(r.array,r.ranges,th.sp.binomials)
{
  threat.code.matrix = array(NA,dim = dim(r.array))
  total.range.matrix = array(NA,dim = dim(r.array))
  cell.intersect.matrix = array(NA,dim = dim(r.array))
  
  for(s in seq_len(ncol(r.array)))
  {
    if(s %% floor(ncol(r.array)/10) == 0) print(s)
    s.binomial= colnames(r.array)[s]
    present.cell.inds = which(!is.na(r.array[,s]))
    
    #Fill in the range area matrix where that species is present in eah grid cell
    total.range.matrix[present.cell.inds,s] = r.ranges$RangeArea[s]
    
    cell.intersect.matrix[present.cell.inds,s] = r.array[present.cell.inds,s]
    
    if(s.binomial %in% th.sp.binomials)
    {
      threat.code.matrix[present.cell.inds,s] = 1
    }
    else
    {
      threat.code.matrix[present.cell.inds,s] = 0
    }
  }
  
  
  
  threat.code.matrix = matrix(threat.code.matrix,ncol = ncol(threat.code.matrix))
  colnames(threat.code.matrix) = colnames(r.array)
  
  total.range.matrix = matrix(total.range.matrix, ncol= ncol(total.range.matrix))
  colnames(total.range.matrix) = colnames(r.array)
  
  cell.intersect.matrix = matrix(cell.intersect.matrix, ncol= ncol(cell.intersect.matrix))
  colnames(cell.intersect.matrix) = colnames(r.array)
  
  
  list(threat.code.matrix = threat.code.matrix,
       total.range.matrix = total.range.matrix,
       cell.intersect.matrix = cell.intersect.matrix)
  
}


ConstructOccupancyMatricesSmall<-function(r.array,r.ranges,th.sp.binomials)
{
  threat.code.matrix = array(NA,dim = dim(r.array))
  
  for(s in seq_len(ncol(r.array)))
  {
    #if(s %% floor(ncol(r.array)/10) == 0) print(s)
    s.binomial= colnames(r.array)[s]
    present.cell.inds = which((r.array[,s]))
    
    
    if(s.binomial %in% th.sp.binomials)
    {
      threat.code.matrix[present.cell.inds,s] = 1
    }
    else
    {
      threat.code.matrix[present.cell.inds,s] = 0
    }
  }
  
  
  
  threat.code.matrix = matrix(threat.code.matrix,ncol = ncol(threat.code.matrix))
  colnames(threat.code.matrix) = colnames(r.array)
  

  list(threat.code.matrix = threat.code.matrix)
  
}






ConstructThreatSeverityScopeMatrices<-function(r.array,ss,r.ranges,th.sp.binomials,
                                               severity.lookup,scope.mid.points,
                                               threat)
{
  #threat.code.matrix = array(NA,dim = dim(r.array))
  severity.score.matrix = array(NA,dim = dim(r.array))
  scope.score.matrix = array(NA, dim = dim(r.array))  
    
  for(s in seq_len(ncol(r.array)))
  {
    if(s %% floor(ncol(r.array)/10) == 0) print(s)
    s.binomial= colnames(r.array)[s]
    
    present.cell.inds = which(r.array[,s] == 1)
    
    
    if(s.binomial %in% th.sp.binomials)
    {
      code.rows = NULL
      for(code in unlist(th))
      {
        
        if(length((grep("\\.", code))) > 0)
        {
          code.rows = c(code.rows,grep(paste0("^",code,"{1}"),ss$ThreatCode))
        } else
        {
          code.rows = c(code.rows,grep(paste0("^",code,"{1}\\."),ss$ThreatCode)) 
        }
        
        
      }
      
      code.rows = intersect(unique(code.rows),
                              which(ss$Scientific.name == s.binomial))
      
      max.ind<-code.rows[which.max(unlist(severity.lookup[ss$Severity[code.rows]]))]
      severity<-severity.lookup[[ss$Severity[max.ind]]]
      scope<-scope.mid.points[[ss$Scope[max.ind]]]
      
      #If threatened then look for the severity and scope value 
      severity.score.matrix[present.cell.inds,s] = severity
      
      scope.score.matrix[present.cell.inds,s] = scope
    }
    else
    {
      severity.score.matrix[present.cell.inds,s] = -1
      
      #c.i<-grep(paste0(s.binomial,"$"),cells.intersected$binomil)
      scope.score.matrix[present.cell.inds,s] = scope.mid.points$`Whole (>90%)`#1/(cells.intersected$cells[c.i]^(1/3))
    }
  }
  
  
  
  severity.score.matrix = matrix(severity.score.matrix,ncol = ncol(severity.score.matrix))
  colnames(severity.score.matrix) = colnames(r.array)
  
  scope.score.matrix = matrix(scope.score.matrix,ncol = ncol(scope.score.matrix))
  colnames(scope.score.matrix) = colnames(r.array)
    
  list(severity.matrix = severity.score.matrix,
       scope.matrix = scope.score.matrix)
  
}




ConstructOccupancyMatricesPresence<-function(r.array,r.ranges,th.sp.binomials)
{
  threat.code.matrix = array(0,dim = dim(r.array))
  total.range.matrix = array(NA,dim = dim(r.array))
  cell.intersect.matrix = array(NA,dim = dim(r.array))
  presence.matrix = array(NA,dim = dim(r.array))
  
  for(s in seq_len(ncol(r.array)))
  {
    if(s %% floor(ncol(r.array)/10) == 0) print(s)
    s.binomial= colnames(r.array)[s]
    present.cell.inds = which(!is.na(r.array[,s]))
    
    #Fill in the range area matrix where that species is present in eah grid cell
    total.range.matrix[present.cell.inds,s] = r.ranges$RangeArea[s]
    
    cell.intersect.matrix[present.cell.inds,s] = r.array[present.cell.inds,s]
    presence.matrix[present.cell.inds,s] = 1
    
    
    if(s.binomial %in% th.sp.binomials)
    {
      threat.code.matrix[present.cell.inds,s] = 1.0
    }
    else
    {
      threat.code.matrix[present.cell.inds,s] = 0.0
    }
  }
  
  
  
  threat.code.matrix = matrix(threat.code.matrix,ncol = ncol(threat.code.matrix))
  colnames(threat.code.matrix) = colnames(r.array)
  
  total.range.matrix = matrix(total.range.matrix, ncol= ncol(total.range.matrix))
  colnames(total.range.matrix) = colnames(r.array)
  
  cell.intersect.matrix = matrix(cell.intersect.matrix, ncol= ncol(cell.intersect.matrix))
  colnames(cell.intersect.matrix) = colnames(r.array)
  
  presence.matrix = matrix(presence.matrix,ncol = ncol(presence.matrix))
  colnames(presence.matrix) = colnames(r.array)
  
  
  list(threat.code.matrix = threat.code.matrix,
       total.range.matrix = total.range.matrix,
       cell.intersect.matrix = cell.intersect.matrix,
       presence.matrix = presence.matrix)
  
}



AggregateToGenus<-function(d)
{
  all.genera = unlist(lapply(strsplit(colnames(d)," ",fixed = T),function(x) x[1]))
  genera = unique(all.genera)
  
  g.d = matrix(NA,nrow = nrow(d),ncol = length(genera))
  g.d.cnt = matrix(NA,nrow = nrow(d),ncol = length(genera))
  g.d.p = matrix(NA,nrow = nrow(d),ncol = length(genera))
  
  for(g in 1:length(genera))
  {
    g.i = which(all.genera == genera[g])
    
    if(length(g.i) > 1)
    {
      g.1s = (apply(d[,g.i],1,function(x) sum(x == 1,na.rm=T)))
      g.0s = (apply(d[,g.i],1,function(x) sum(x == 0, na.rm=T)))
    }
    else
    {
      g.1s = d[,g.i]
      g.0s = d[,g.i]
    }
    
    g.cnt = g.0s + g.1s
    
    # Presence of the group in cells
    g.d[which(g.1s > 0),g] = 1
    
    #Presence of the threat in cells
    g.d[which(g.0s > 0),g] = 0
    
    g.d.cnt[,g] = g.cnt
    g.d.p[,g] = g.0s/g.cnt
  }
  
  list(
    d = g.d,
    cnt = g.d.cnt,
    p = g.d.p
  )
  
}
