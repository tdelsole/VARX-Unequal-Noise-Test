rm(list=ls())

lplotfile    = FALSE

nspace       = 1
num.forcings = 3
p.order      = 2
nharm        = 6
mperiod      = 1
nens.max     = 3
nlon         = 90
nlat         = 89
niter        = 2
first.say    = 'ar'
first.say    = 'forcing'


dir.atlas    = '/Users/delsole/Documents/data/CMIP6_AtlasRegions/'
dir.Rlib     = '/Users/delsole/documents/R/delsole_tools/'
dir.Rlib     = NULL
dir.atlas    = NULL

source(paste(dir.Rlib,'pdf.eps.R',sep=''))
source(paste(dir.Rlib,'deviance.AnyAR.FixNoise.NoiseNotFirst.R',sep=''))
source(paste(dir.Rlib,'simulate.AnyAR.FixNoise.1forcing.R',sep=''))
source(paste(dir.Rlib,'gev.R',sep=''))
source(paste(dir.Rlib,'mic.master.R',sep=''))
source(paste(dir.Rlib,'LjungBox.R',sep=''))
source(paste(dir.Rlib,'timeseries2Varx.cyclo.R',sep=''))
source(paste(dir.Rlib,'n.to.monthly.R',sep=''))
source(paste(dir.Rlib,'plot_latlon_v4.R',sep=''))
source(paste(dir.Rlib,'AtlasMaskProduce.R',sep=''))
source(paste(dir.Rlib,'anncyc.remove.monthly.R',sep=''))
source(paste(dir.Rlib,'diff.regression.nested.mult.NoiseNotFirst.R',sep=''))
source(paste(dir.Rlib,'sync.periods.R',sep=''))

library(lubridate)
library(fields)
library(maps)
library(rworldmap)

if (nspace == 1) {
  title.meta = paste('ARX(',p.order,'); ',nspace,' region',sep='')
} else {
  title.meta = paste('ARX(',p.order,'); ',nspace,' regions',sep='')
}




###########################################
######## READ ATLAS REGIONS
###########################################
fname.atlas = paste(dir.atlas,'IPCC-WGI-reference-regions-v4_R.rda',sep='')
load(fname.atlas)
region.name.short = IPCC_WGI_reference_regions_v4[["Acronym"]]
region.name.long  = IPCC_WGI_reference_regions_v4[["Name"]]
nregions          = nrow(IPCC_WGI_reference_regions_v4@data)
region.type       = IPCC_WGI_reference_regions_v4$Type

lon               = seq(from=0,to=360,length.out=nlon+1)[-(nlon+1)]
lat               = seq(from=-90,to=90,length.out = nlat+2)[-c(1,nlat+2)]

mask.list = AtlasMaskProduce(IPCC_WGI_reference_regions_v4, lon, lat)

ipic = 58
plot_latlon_v4(lon,lat,mask.list$lgood[,ipic])

lon.stack  = rep(lon,nlat)
lat.stack  = rep(lat,each=nlon)

index.n.l = NULL
index.s.l = NULL
index.n.o = NULL
index.s.o = NULL
index.n   = NULL
index.s   = NULL
index.l   = NULL
index.o   = NULL
index.glo = 1:nregions
for (n in 1:nregions) if (all(lat.stack[mask.list$lgood[,n]] > 0) & all(lat.stack[mask.list$lgood[,n]] <=  90) & region.type[n] == 'Land' ) index.n.l = c(index.n.l,n)
for (n in 1:nregions) if (all(lat.stack[mask.list$lgood[,n]] < 0) & all(lat.stack[mask.list$lgood[,n]] >= -90) & region.type[n] == 'Land' ) index.s.l = c(index.s.l,n)
for (n in 1:nregions) if (all(lat.stack[mask.list$lgood[,n]] > 0) & all(lat.stack[mask.list$lgood[,n]] <=  90) & region.type[n] == 'Ocean') index.n.o = c(index.n.o,n)
for (n in 1:nregions) if (all(lat.stack[mask.list$lgood[,n]] < 0) & all(lat.stack[mask.list$lgood[,n]] >= -90) & region.type[n] == 'Ocean') index.s.o = c(index.s.o,n)

for (n in 1:nregions) if (all(lat.stack[mask.list$lgood[,n]] > 0) & all(lat.stack[mask.list$lgood[,n]] <=  90) ) index.n = c(index.n,n)
for (n in 1:nregions) if (all(lat.stack[mask.list$lgood[,n]] < 0) & all(lat.stack[mask.list$lgood[,n]] >= -90) ) index.s = c(index.s,n)

for (n in 1:nregions) if ( region.type[n] == 'Land'  ) index.l = c(index.l,n)
for (n in 1:nregions) if ( region.type[n] == 'Ocean' ) index.o = c(index.o,n)

index.e.o = c(49,52,56)

### REMOVE TROPICS
index.o   = setdiff(index.o,index.e.o)

if (nspace == 5) {
  space.list = list('S Land'  = index.s.l, 
                    'N Land'  = index.n.l, 
                    'S Ocean' = index.s.o, 
                    'N Ocean' = index.n.o, 
                    'E Ocean' = index.e.o)
} else if (nspace == 1) {
  space.list = list('GLO' = index.glo)
} else if (nspace == 2) {
  space.list = list('N Hem' = index.n, 'S Hem' = index.s)
  space.list = list('Land'  = index.l, 'Ocean' = index.o)
} else if (nspace == 3) {
  space.list = list('N Hem' = index.n, 'S Hem' = index.s, 'E Ocean' = index.e.o)
  space.list = list('Land'  = index.l, 'Ocean' = index.o, 'E Ocean' = index.e.o)
} else stop()


if (nspace != length(space.list)) stop('space.list has inconsistent length')


suffix.meta = paste('S',nspace,'.H',nharm,'.M',mperiod,
                    '.F',num.forcings,'.P',p.order,'.1st',first.say,sep='')


space.map  = as.numeric(rep(NA,nlon*nlat))
for (ns in 1:nspace) {
  lgood = rowSums(mask.list$lgood[,space.list[[ns]]])
  lgood = lgood > 0
  space.map[lgood] = ns - 0.5 
}

fout = paste('./figures.twins/Atlas.S',nspace,'.map',sep='')
if (lplotfile) pdf.eps(fout,'pdf',width=8,height=5)
plot.breaks = seq(from=0,by=1,length.out=nspace+1)
plot.cols   = 1:(nspace)
cbar.list = list(plot.breaks = plot.breaks, plot.cols = plot.cols)
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
par(mfcol=c(1,1),mar=c(9,5,2,1))
plot_latlon_v4(lon,lat,space.map,cbar.list = cbar.list)
title(main='Regions Included in Deviance Calculation',line=0.5)
if (lplotfile) dev.off()

compress.regions = function(pc) {
  for (ns in 1:nspace) {
    iget    = space.list[[ns]]
    pc.dum  = colSums( t(pc)[iget,] * mask.list$weight.per.region[iget]) / sum(mask.list$weight.per.region[iget])
    if (ns == 1) pc.ts = array(NA,dim=c(length(pc.dum),nspace))
    pc.ts[,ns] = pc.dum
  }
  pc.ts
}



#######################################
##### READ OBSERVATIONS: ERA5
#######################################
fread = paste(dir.atlas,'ipcc.global.2t_Amon_CMIP6_ERA5.consolidate2.RData',sep='')

load(fread)
obs.list  = master.list
obs.ts    = compress.regions(obs.list$model.list$ERA5$pc.val)
obs.targ  = obs.list$model.list$ERA5$pc.time
ntime.obs = length(obs.targ)

obs.targ      = as.Date(obs.targ)
day(obs.targ) = 15


#######################################
##### PLOT OBSERVATIONS: ERA5
#######################################
fout = paste('./figures.twins/',suffix.meta,'.obs',sep='')
if (lplotfile) pdf.eps(fout,'pdf',width=8,height=8)
nrow.plot  = ceiling(nspace/2)
par(mfrow = c(nrow.plot,2),mar=c(3,3,2,0.1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
for (ns in 1:nspace) {
  plot(obs.targ,anncyc.remove.monthly(obs.ts[,ns]),type='l',xlab='',ylab='')
  title(main=names(space.list)[ns],line=0.5)
}
if (lplotfile) dev.off()

#######################################
##### READ HISTORICAL ATLAS REGIONS
#######################################
fread = paste(dir.atlas,'ipcc.global.ts_Amon_CMIP6_historical.consolidate.RData',sep='')
load(fread)
hist.list = master.list


names.hist = names(hist.list$model.list)
nmod.hist  = length(names.hist)

### REPAIR THE DATES
for (nm in 1:nmod.hist) {
  hist.list$model.list[[nm]]$pc.time = as.Date(hist.list$model.list[[nm]]$pc.time)
  day(hist.list$model.list[[nm]]$pc.time)  = 15
}

#### DEFINE OVERLAPPING TIME PERIOD
sync.list = sync.periods(obs.targ,hist.list$model.list[[1]]$pc.time)

obs.ts    = obs.ts  [sync.list$nst1:sync.list$nnd1,]
obs.targ  = obs.targ[sync.list$nst1:sync.list$nnd1]
ntime.obs = length(obs.targ)

ntot.hist = sync.list$ndates.common






nens.hist = NULL
pc.hist   = NULL
for (nm in 1:nmod.hist) {
  iget = hist.list$model.list[[nm]]$pc.time >= sync.list$dates.common[1] &
         hist.list$model.list[[nm]]$pc.time <= sync.list$dates.common[sync.list$ndates.common]
  if (sum(iget) != ntot.hist) stop('incomplete historical run')
  nens = length(hist.list$model.list[[nm]]$ens.unique)

  ndim = length(dim(hist.list$model.list[[nm]]$pc.val))
  pc   = array(NA,dim=c(ntot.hist,nspace,nens))
  if (ndim == 3) {
    for (ne in 1:nens) pc[,,ne] = compress.regions(hist.list$model.list[[nm]]$pc.val[iget,,ne])
  } else if (ndim == 2) {
    pc[,,1] = compress.regions(hist.list$model.list[[nm]]$pc.val[iget, ])
  } else stop('dimension of pc.val is not 2 or 3')

  hist.list$model.list[[nm]]$pc.space = pc
  hist.list$model.list[[nm]]$nens     = nens
  for (ne in 1:min(nens.max,nens)) pc.hist = cbind(pc.hist,pc[,,ne])
  
  nens.hist = c(nens.hist,min(nens.max,nens))

}

hist.targ = hist.list$model.list[[nm]]$pc.time[iget]

nmod.hist.all = sum(nens.hist)
if (length(pc.hist) != ntot.hist * nspace * nmod.hist.all) stop('pc.hist incorrectly dimensioned')
dim(pc.hist) = c(ntot.hist,nspace,nmod.hist.all)

names.all.hist    = NULL
for (nm in 1:nmod.hist) names.all.hist = c(names.all.hist,rep(names.hist[nm],nens.hist[nm]))









#######################################################
######### READ FORCING: forcing.list
#######################################################
fread = paste('/Users/delsole/Documents/data/CMIP6_forcing/CMIP6.forcing.Annex.AllScenarios.F',num.forcings,'.RData',sep='')
load(fread)

forcing.say   = forcing.list[[1]]$forcing.say
forcing.pc    = forcing.list[[1]]$forcing.mon
forcing.targ  = forcing.list[[1]]$forcing.targ
nforcing      = length(forcing.say)
ntime.forcing = length(forcing.targ) 

#######################################################
######### REMOVE TIME MEAN FROM THE FORCINGS
#######################################################
forcing.pc = as.matrix(forcing.pc)
forcing.pc = t( t(forcing.pc) - colMeans(forcing.pc) )


#######################################################
######### COMPUTE DEVIANCES
#######################################################
print(Sys.time())

if (first.say == 'noise') {
  test.order.call = c('noise','forcing','ar','cycle','intercept')
} else if (first.say == 'forcing') {
  test.order.call = c('forcing','ar','cycle','intercept')
} else if (first.say == 'ar') {
  test.order.call = c('ar','forcing','cycle','intercept')
} else stop('do not understand first.say')

noise.var.diag  = array(NA,dim=c(nspace,nmod.hist.all+1,nmod.hist.all+1))
noise.var.gev   = array(NA,dim=c(nspace,nmod.hist.all+1,nmod.hist.all+1))
dev.hist        = array(NA,dim=c(length(test.order.call   )+1,nmod.hist.all+1,nmod.hist.all+1))
dev.eq.noise    = array(NA,dim=c(niter+1,length(test.order.call),nmod.hist.all+1,nmod.hist.all+1))
for (nm1 in 0:nmod.hist.all) for (nm2 in nm1:nmod.hist.all) if (nm1 != nm2) {
  if (nm2 == nm1 + 1) print(paste(nm1,'out of',nmod.hist.all))
  if (nm1 == 0) {
    pc1    = obs.ts 
    ntime1 = ntime.obs
    targ1  = obs.targ
  } else {
    pc1    = pc.hist[,,nm1]
    ntime1 = ntot.hist
    targ1  = hist.targ
  }
  if (nm2 == 0) {
    pc2    = obs.ts 
    ntime2 = ntime.obs
    targ2  = obs.targ
  } else {
    pc2    = pc.hist[,,nm2]
    ntime2 = ntot.hist
    targ2  = hist.targ
  }
  
  if (first.say == 'noise') {
    dev.list = deviance.AnyAR.FixNoise.1Forcing(
      pc1,pc2,targ1,targ2,ntime1,ntime2,1,1,nspace,
      forcing.pc,forcing.targ,ntime.forcing,nforcing,forcing.say,
      p.order,nharm,test.order.call=test.order.call)
  } else {
    dev.list = deviance.AnyAR.FixNoise.NoiseNotFirst(
      pc1,pc2,targ1,targ2,ntime1,ntime2,1,1,nspace,
      forcing.pc,forcing.targ,ntime.forcing,nforcing,forcing.say,
      p.order,nharm, niter = niter, test.order.call = test.order.call)
    dev.eq.noise[,,nm1+1,nm2+1] = dev.list$dev.eq.noise
  }
  dev.hist[ ,nm1+1,nm2+1] = dev.list$dev.table[,'deviance']
  noise.var.diag[,nm1+1,nm2+1] = diag(dev.list$qmat1)/diag(dev.list$qmat2)
  noise.var.gev [,nm1+1,nm2+1] = gev(dev.list$qmat1,dev.list$qmat2)$lambda
  
}
print(Sys.time())



### SYMMETRIZE THE DEVIANCE MATRIX
for (nm1 in 0:nmod.hist.all) for (nm2 in nm1:nmod.hist.all) if (nm1 != nm2) dev.hist[,nm2+1,nm1+1] = dev.hist[,nm1+1,nm2+1]
for (nm1 in 0:nmod.hist.all) dev.hist[,nm1+1,nm1+1] = 0

for (nm1 in 0:nmod.hist.all) for (nm2 in nm1:nmod.hist.all) if (nm1 != nm2) noise.var.diag[,nm2+1,nm1+1] = noise.var.diag[,nm1+1,nm2+1]
for (nm1 in 0:nmod.hist.all) for (nm2 in nm1:nmod.hist.all) if (nm1 != nm2) noise.var.gev [,nm2+1,nm1+1] = noise.var.gev [,nm1+1,nm2+1]

dev.crit     = dev.list$dev.table[,'crit(chisq)']
dev.rownames = rownames(dev.list$dev.table)
dev.dof      = dev.list$dev.table[,'dof']

dev.names    = 'ERA5'
for (nm in 1:nmod.hist) dev.names = c(dev.names,rep(names.hist[nm],nens.hist[nm]))
nmod.dev     = length(dev.names)
matrix.break = cumsum(c(1,nens.hist))

if (first.say != 'noise') dimnames(dev.eq.noise) = list(dimnames(dev.list$dev.eq.noise)[[1]],
                              dimnames(dev.list$dev.eq.noise)[[2]],
                              dev.names,dev.names)

##### STRIP OUT FIRST 7 CHARACTERS IN DEV.ROWNAMES
dev.rownames.new = NULL
for (n in 1:length(dev.rownames)) {
  name.strip       = substr(dev.rownames[n],8,nchar(dev.rownames[n]))
  dev.rownames.new = c(dev.rownames.new,name.strip)
}
dev.rownames = dev.rownames.new


#################################################
####### PLOT A FEW CASES OF ITERATIVE SOLUTIONS
#################################################
suffix  = paste('.S',nspace,'.H',nharm,'.M',mperiod,'.F',num.forcings,'.P',p.order,sep='')
num.hyp = 2
num.ref = 1




##########################################################
###### PLOT HEAT MAP FOR JUST NOISE + AR + FORCING
###### INCLUDE ENSEMBLE MEMBERS IN DEVIANCE ARRAY

###### ENSURE PROJECTION TIME SERIES INCLUDES THE BASE MODEL!
###### PLOT PROJECTIONS FOR TWINS
###### COMPUTE DEVIANCE FOR AR + FORCING ONLY (NO NOISE)
###### COMPUTE F-RATIOS
###### REPRODUCE HEAT MAP FOR S = 1 (GLOBAL MEAN)
###### TRY FEWER REGIONS: N+S; N+S+TROPICS
###### INCLUDE SPEAR DATA
###### INCLUDE ATLAS FOR BERKELEY DATA.
###### NOISE-FREE SIMULATIONS FOR EACH MODEL
###### MODIFY DEVIANCE TO NOT START WITH NOISE COMPARISON
###### COMPARE ATLAS REGIONS WITH GITHUB DATA
#########################################################





########################################
######### PLOT PAIR-WISE DEVIANCES (HEAT MAPS, DENDOGRAM)
########################################
alpha.subdev = 1.e-5
alpha.subdev = 0.01
alpha.subdev = 0.05

suffix.alpha = paste('S',nspace,'.H',nharm,'.M',mperiod,
                     '.F',num.forcings,'.P',p.order,
                     '.1st',first.say,
                     paste('.a',alpha.subdev,sep=''),sep='')

for (dtype in c('ebm','arpar','forcing','cycle')) {
  if (dtype == 'total') {
    rowpic = which(dev.rownames == 'Total')
    dof    = dev.dof[rowpic]
    dtype.say = 'Total'
  } else if (dtype == 'arpar') {
    rowpic = which(dev.rownames == 'AR')
    dof    = dev.dof[rowpic]
    dtype.say = 'AR'
  } else if (dtype == 'noise') {
    rowpic = which(dev.rownames == 'Noise')
    dof    = dev.dof[rowpic]
    dtype.say = 'Noise'
  } else if (dtype == 'cycle') {
    rowpic = which(dev.rownames == 'Cycle')
    dof    = dev.dof[rowpic]
    dtype.say = 'Annual Cycle'
  } else if (dtype == 'forcing') {
    rowpic = which(dev.rownames == 'Forcing')
    dof    = dev.dof[rowpic]
    dtype.say = 'Forcing'
  } else if (dtype == 'intercept') {
    rowpic = which(dev.rownames == 'Intercept')
    dof    = dev.dof[rowpic]
    dtype.say = 'Intercept'
  } else if (dtype == 'ebm') {
    row.noise = which(dev.rownames == 'Noise')
    row.ar    = which(dev.rownames == 'AR')
    row.for   = which(dev.rownames == 'Forcing')
    dtype.say = 'EBM'

    dmat      = dev.hist[row.ar,,] + dev.hist[row.for,,]
    dof       = sum(dev.dof[c(row.ar,row.for)])
    dcrit     = qchisq(alpha.subdev,dof,lower.tail=FALSE)
    
  } else stop('do not understand dtype')
  
  if (dtype != 'ebm') {
    if (length(rowpic) != 1) stop('could not identify proper row in deviance')
    dmat  = dev.hist[rowpic,,]
    dcrit = dev.crit[rowpic]   
    dcrit = qchisq(alpha.subdev,dof,lower.tail=FALSE)
  } 
  
  
  dmat.raw = dmat
  dmat.raw [dmat.raw < 1.e-10] = 0
  
  if (max(abs(dmat - t(dmat)),na.rm=TRUE) > 1.e-6) stop('dmat not symmetric!')
  rownames(dmat.raw) = dev.names
  colnames(dmat.raw) = dev.names
  

  
  
  ### (PRETTY 2-PANEL) PLOT DEVIANCE AND CRITICAL VALUES RELATIVE TO OBSERVATIONS
  mref   = which(dev.names == 'ERA5')
  y      = dmat [ ,mref]
  y.crit = dcrit
  yrange = range(y,y.crit)
  xrange = c(1:nmod.dev)
  norder = order(y)
  nhalf  = floor(nmod.dev/2)
  n1st   = norder[1:nhalf]
  n2nd   = norder[(nhalf+1):nmod.dev]
  yrange = range(y)
  matrix.layout = matrix(c(1,1,2,3),byrow=TRUE,ncol=2)
  if (dtype == 'forcing') {
    ftitle.top = paste('Test Equality of Transfer Coefficients for Radiative Forcing')
  } else {
    ftitle.top = paste('Test Equality of',dtype.say,'Coefficients Relative to Observations')
  }
  fout = paste('./figures.twins/',suffix.meta,'.dev.line.',dtype,'2panel',sep='')
  if (lplotfile) pdf.eps(fout,'pdf',width=8.5,height=8)
  layout(matrix.layout,height=c(1,6))
  par(mar=rep(0.1,4))
  plot(0,0,type='n',xlab='',ylab='',axes=FALSE)
  ftitle.bot = paste(title.meta,sep='')
  text.say   = paste(ftitle.top,'\n',ftitle.bot)
  text(0,0,text.say,cex=1.5)
  par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
  par(mar=c(5,8,0.1,0.2))
  plot(y[n1st],1:nhalf,cex=0.5,pch=19,xlab='deviance wrt observations',ylab='',yaxt='n',xlim=yrange)
  par(cex.axis=0.7)
  axis(2,at=1:nhalf,dev.names[n1st],las=2)
  abline(v=y.crit[1],col='red')
  abline(v=y.crit[2],col='red')
  text(y.crit[1],0.1,'significant difference',col='red',pos=4)
  par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
  par(mar=c(5,8,0.1,0.2))
  plot(y[n2nd],1:length(n2nd),cex=0.5,pch=19,xlab='deviance wrt observations',ylab='',yaxt='n',xlim=yrange)
  par(cex.axis=0.7)
  axis(2,at=1:length(n2nd),dev.names[n2nd],las=2)
  abline(v=y.crit[1],col='red')
  abline(v=y.crit[2],col='red')
  text(y.crit[1],0.1,'significant difference',col='red',pos=4)
  if (lplotfile) dev.off()
  
  
  ### PLOT HEAT MAPS
  fout = paste('./figures.twins/',suffix.alpha,'.matrix.all.',dtype,sep='')
  if (lplotfile) pdf.eps(fout,'pdf',width=8,height=8)
  par(mfcol=c(1,1),mar=c(5,5,5,5))
  col.image = c('white','grey70')
  brk.image = c(0,1,100000)
  dmat.plot = dmat.raw[,(nmod.dev):1]/dcrit
  image(1:(nmod.dev),1:(nmod.dev),dmat.plot,axes=FALSE,xlab='',ylab='',col=col.image,breaks=brk.image)
  axis(1,1:(nmod.dev),rownames(dmat.raw),cex.axis=0.2,las=3)
  axis(2,1:(nmod.dev),rev(colnames(dmat.raw)),cex.axis=0.2,las=1)
  abline(v=matrix.break+0.5,col='red')
  abline(h=nmod.dev-matrix.break+0.5,col='red')
  ftitle.top = paste('Matrix of ',dtype.say,' Deviances',sep='')
  ftitle.bot = paste(title.meta,'; alpha= ',alpha.subdev,sep='')
  title(ftitle.top,line=2.0)
  title(ftitle.bot,line=0.5)
  if (lplotfile) dev.off()

  

}


