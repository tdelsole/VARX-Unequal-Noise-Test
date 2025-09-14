anncyc.remove.monthly = function(z.ts,leave.mean=FALSE) {
### REMOVE ANNUAL CYCLE FROM MONTHLY DATA
### INPUT:  Z.TS     [NTIME,NSPACE]
### OUTPUT: Z.TS.ANOM[NTIME,NSPACE]
  
  z.ts     = as.matrix(z.ts)
  z.mean   = colMeans(z.ts)
  ntot     = dim(z.ts)[1]
  nspace   = length(z.ts) / ntot
  if (ntot %% 12 != 0) {
    ncomplete = ceiling(ntot/12) * 12
    ngap      = ncomplete - ntot
    z.ts      = rbind(z.ts,array(NA,dim=c(ngap,nspace)))
  } else {
    ncomplete = ntot
  }
  nyrs = ncomplete / 12
  dim(z.ts) = c(12,nyrs,nspace)
  for (nv in 1:nspace) z.ts[,,nv] = z.ts[,,nv] - rowMeans(z.ts[,,nv],na.rm=TRUE)
  dim(z.ts) = c(ncomplete,nspace)
  z.ts = z.ts[1:ntot,]
  
  if (leave.mean) z.ts = t (t(z.ts) + z.mean)
  
  return(z.ts)
}