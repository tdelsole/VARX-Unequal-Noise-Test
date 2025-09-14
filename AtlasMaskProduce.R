AtlasMaskProduce = function( PolyDataFrame, grid.lon, grid.lat, field = NA ) {
## THIS FUNCTION GENERATES MASK FILES CORRESPONDING TO EACH CMIP6 ATLAS REGION
## INPUT:
##  PolyDataFrame: A SpatialPolygonsDataFrame corresponding to the CMIP Atlas Regions
##  grid.lon[nlon]: longitudes of the CMIP model grid
##  grid.lat[nlat]: latitudes  of the CMIP model grid
##  field[nlon,nlat,ntime]: data field to be projected onto the regions (ignored if = NA)
## OUTPUT:
  
  nlon = length(grid.lon)
  nlat = length(grid.lat)
  
  region.name.short = PolyDataFrame[["Acronym"]]
  region.name.long  = PolyDataFrame[["Name"]]
  nregions          = nrow(PolyDataFrame@data)
  
  lon.stack  = rep(grid.lon,nlat)
  lat.stack  = rep(grid.lat,each=nlon)
  
  ###########################################
  ######## CONSTRUCT MASK FOR EACH REGION
  ###########################################
  range.lon      = range(lon)
  grid.0.360     = range.lon[1] >=    0 & range.lon[2] <= 360
  grid.n180.p180 = range.lon[1] >= -180 & range.lon[2] <= 180
  if (!grid.0.360 & !grid.n180.p180) stop("longidtude range is in neither (0,360) nor (-180,180)")
  
  lgood = array(NA,dim=c(nlon*nlat,nregions))
  for (nr in 1:nregions) {
    npoly      = length(PolyDataFrame@polygons[[nr]]@Polygons)
    mask.temp  = rep(0,nlon*nlat)
    for (np in 1:npoly) {
      poly.coord = coordinates(PolyDataFrame@polygons[[nr]]@Polygons[[np]])
      # print(paste(nr,region.name.long[nr],npoly,paste(range(poly.coord[,1]),collapse=', '),paste(range(poly.coord[,2]),collapse=', '),sep='; '))
      shift.grid.none = (grid.0.360 & (max(poly.coord[,1]) > 0)) | grid.n180.p180
      shift.grid.left =  grid.0.360 & (min(poly.coord[,1]) < 0)
      if (shift.grid.none) mask.temp = mask.temp + point.in.polygon(lon.stack    ,lat.stack,poly.coord[,1],poly.coord[,2],mode.checked=FALSE)
      if (shift.grid.left) mask.temp = mask.temp + point.in.polygon(lon.stack-360,lat.stack,poly.coord[,1],poly.coord[,2],mode.checked=FALSE)    
    }
    lgood[,nr]  = as.logical(mask.temp > 0.5) ## CONVERT BINARY TO LOGICAL
    if (all(!lgood[,nr])) stop('region has no valid points')
  }
  
  ### ADD AN EXTRA MASK TO ACCOUNT THE GAPS, IF NECESSARY
  gap.points           = rowSums(lgood[,1:nregions]) < 0.5
  if (any(gap.points)) {
    lgood              = cbind(lgood,gap.points)
    region.name.short  = c(region.name.short,'gaps')
    region.name.long   = c(region.name.short,'gaps')    
    nregions           = nregions + 1
  } 

  
  #### DEFINE WEIGHTS FOR GLOBAL MEAN (WEIGHT.PER.CELL) AND REGIONAL MEAN (WEIGHT.PER.REGION)
  weight.per.cell   = 1/rowSums(lgood) * cos(lat.stack * pi / 180)
  weight.per.region = as.numeric(rep(NA,nregions))
  for (nr in 1:nregions) weight.per.region[nr] = sum(weight.per.cell[lgood[,nr]])
  
 
  #### COMPUTE TIME SERIES
  if (!all(is.na(field))) {
    if (length(field) %% (nlon * nlat) != 0 ) stop("field grid inconsistent with lon/lat provided")
    nlead      = length(field) / (nlon * nlat)
    pc         = array(NA,dim=c(nlead,nregions))
    dim(field) = c(nlon*nlat,nlead)
    for (nr in 1:nregions) pc[,nr] = colSums(field[lgood[,nr],,drop=FALSE] * weight.per.cell[lgood[,nr]]) / weight.per.region[nr]   
  } else pc = 'field was not provided'
  


  return(list(weight.per.cell = weight.per.cell, lgood = lgood, pc   = pc,   neof = nregions,
              lon  = lon,  lat   = lat,   nlon = nlon, nlat = nlat,
              weight.per.region = weight.per.region, 
              region.name.short = region.name.short, 
              region.name.long = region.name.long)) 
  
}