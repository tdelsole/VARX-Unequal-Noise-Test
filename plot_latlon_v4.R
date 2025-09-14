plot_latlon_v4 = function(lon,lat,field,nbreaks=10,plot.legend=TRUE,title.say=NULL,
    lmaskzero=TRUE,shrinkdomain=FALSE,plot.continents=TRUE,force.pm=FALSE,
    cbar.list=NULL,suppress.axislab=FALSE,rotate.land=FALSE,labelFontSize=1.5,include.zero=TRUE) {
	
  ## GENERATES SPATIAL PLOTS ON A GLOBE WITH SUPERIMPOSED CONTINENTAL MAP
  ## lon[nlon]: longitude coordinates
  ## lat[nlat]: latitude coordinates
  ## field[nlon,nlat]: the spatial field to plot
  ## nbreaks: number of break points for color bar (default=10)
  ## plot.legend: =true to plot the legend, =false to suppress legend
  ## title.say: character string for the title
  ## lmaskzero: =true to set zero to 'white', otherwise set to 'grey'
  ## shrinkdomain: =true to re-adjust lat/lon to show only defined points
  ## plot.continents: = true to plot continents, otherwise continents are not plotted
  ## force.pm: force color bar to range from negative to positive values
  ## cbar.list: a list that specifies color bar and break points (optional)
  ## suppress.axislab: if true, lat and lon coordinates will not be plotted
  ## rotate.land: shift the longitude so that land is not cut off
  ## labelFontSize: font size of the axis legend
  ## include.zero: force color bar to include 0 (TRUE) or adapt to the data (FALSE)
  

nlon = length(lon)
nlat = length(lat)
if ( length(field) %% (nlon*nlat) != 0) stop('field not dimensioned correctly')

all.na        = all(is.na(field))
# if (!all.na) all.constant  = (var(as.numeric(field),na.rm=TRUE)/max(abs(field),na.rm=TRUE) <= 1.e-15 ) else all.constant = TRUE
if (!all.na) {
	if (all(field[!is.na(field)] == 0)) {
		all.constant = TRUE
	} else {
		all.constant  = (var(as.numeric(field),na.rm=TRUE)/max(abs(field),na.rm=TRUE) <= 1.e-15 ) 
	}
} else all.constant = TRUE

nbreaks.half  = floor(nbreaks/2)

#####################################################################
## IMPORT COLOR BAR PARAMETERS
#####################################################################
if (!is.null(cbar.list)) {
	plot.cols   = cbar.list$plot.cols
	plot.breaks = cbar.list$plot.breaks
}

#####################################################################
## IF FIELD IS ALL NAS, MAKE UP DATA SO THAT A BLANK PLOT IS FORMED
#####################################################################
if (all.na) {
	dim(field) = c(length(lon),length(lat))
	field[1,1] = 0
	plot.breaks = c(-1,0,1)
	plot.cols   = c("white","white")	

} else if (is.null(cbar.list)) {
	#####################################################################
	## IF FIELD IS A CONSTANT, THEN MAKE UP A RANGE OF DATA
	#####################################################################
	if ( all.constant) {
		constant.value = mean(as.numeric(field),na.rm=TRUE)
		if (constant.value == 0 ) {
			plot.breaks = seq(from=0,to=1,length.out=nbreaks.half) 
		} else { 
			plot.breaks = seq(from = 0, to = constant.value,length.out = nbreaks.half)
			plot.breaks[nbreaks.half] = plot.breaks[nbreaks.half]*1.001
		}
	} else {
		plot.breaks = pretty(abs(field),n=floor(nbreaks.half))
	}
	
	#####################################################################
	## GENERATE SYMMETRIC BREAKS ABOUT 0, OR INSERT ZERO IN THE BREAKS
	#####################################################################
	l.pm   = min(field,na.rm=TRUE)*max(field,na.rm=TRUE) < 0
	if (l.pm | force.pm) {
		plot.breaks = c(-rev(plot.breaks),plot.breaks)
		plot.breaks = plot.breaks[ plot.breaks != 0 ] ; # get rid of zero break points
	} else {
		if (include.zero) plot.breaks = c(0,plot.breaks)
	}
	plot.breaks = sort(union(plot.breaks,plot.breaks)); #get rid of redundant breaks
	
	#####################################################################
	## GENERATE COLORS
	#####################################################################
	if (!all.na) {
		nb          = length(plot.breaks)
		if (l.pm) {
			plot.cols   = colorRampPalette(c('blue','deepskyblue','white','orange','red'),space='rgb')(nb-1)
		} else {
			plot.cols   = colorRampPalette(c('orange','red','darkred'),space='rgb')(nb-1)			
		}
	}
	
	#####################################################################
	## SET ZERO INTERVAL TO WHITE IF DESIRED
	#####################################################################
	nb                 = length(plot.breaks)
	break.include.zero = plot.breaks[1:(nb-1)] <=0 & plot.breaks[2:nb] >= 0
	if (lmaskzero) {
		plot.cols[break.include.zero] = 'white'
	} else {
		plot.cols[break.include.zero] = 'grey'
	}
	
}


#####################################################################
## DEFINE DUMMY ARRAY FOR PLOTTING
#####################################################################
plot.var = array(field,dim=c(length(lon),length(lat)))

#####################################################################
## REVERSE LATITUDES IF THEY GO FROM NORTH TO SOUTH
#####################################################################
if ( lat[2] - lat[1] < 0 ) {
	plot.var = plot.var[,length(lat):1]
	lat.say  = lat[length(lat):1]
	print('latitudes reversed')
} else {
	lat.say  = lat
}

#####################################################################
## ROTATE TO AVOID CUTTING ACROSS LAND MASSES
#####################################################################
if (rotate.land & lon[1] >= 0 ) {
	nlon     = length(lon)
	nst      = max(which(lon <= 180))
	npic     = c((nst+1):nlon,1:nst)
	plot.var = plot.var[npic,]
	lon      = c(lon[(nst+1):nlon] - 360,lon[1:nst])
}


#####################################################################
## SHRINK DOMAIN IF DATA IS MISSING 
#####################################################################
if (suppress.axislab) {xaxt='n';yaxt='n'} else {xaxt='s';yaxt='s'}
if (shrinkdomain) {
	good  = which(!is.na(plot.var),arr.ind=TRUE)
	xpic  = range(good[,'row'])
	ypic  = range(good[,'col'])
	xlim = lon[xpic]; ylim = lat.say[ypic]
	image(lon,lat.say,plot.var,col=plot.cols,breaks=plot.breaks,xlab="",ylab="",xlim=xlim,ylim=ylim,xaxt=xaxt,yaxt=yaxt)
} else {
	image(lon,lat.say,plot.var,col=plot.cols,breaks=plot.breaks,xlab="",ylab="",xaxt=xaxt,yaxt=yaxt)
}


#####################################################################
## ADD COLOR BAR 
#####################################################################
if (plot.legend) addMapLegend(cutVector=plot.breaks,colourVector=plot.cols,legendLabels='all',
   legendMar=5,horizontal=TRUE,labelFontSize=labelFontSize,catMethod='pretty')

#####################################################################
## ADD TITLE
#####################################################################
if (!is.null(title.say)) title(main=as.character(title.say))

#####################################################################
## ADD CONTINENTS
#####################################################################
if (plot.continents) if ( lon[1] < 0 ) map('world',interior=FALSE,add=TRUE) else map('world2',interior=FALSE,add=TRUE)
# if (plot.continents) if ( lon[1] < 0 ) map('world',interior=FALSE,add=TRUE,col='gray',fill=TRUE) else map('world2',interior=FALSE,add=TRUE)


#####################################################################
## OUTPUT PLOTTING PARAMETERS
#####################################################################
list(plot.breaks=plot.breaks,plot.cols=plot.cols)


}