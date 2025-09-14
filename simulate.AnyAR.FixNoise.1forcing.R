simulate.AnyAR.FixNoise.1forcing = function(initial.ts,initial.date,nlead,
        f.ts,f.date,fnames,nfor,dev.list,which.model,add.noise=TRUE,suppress.anncyc=FALSE) {
### SIMULATE TIME SERIES USING FITTED CYCLOSTATIONARY MODEL
### INPUT:
###		INITIAL.TS[P.ORDER,NSPACE]: INITIAL VALUES OF THE AR PROCESS
###		INITIAL.DATE  [1]: DATE OF INITIAL.TS[1]. SUBSEQUENT DATES ONE MONTH APART 
###   NLEAD: NUMBER OF STEPS TO INTEGRATE ('MAXIMUM LEAD')
###		F.TS   [FTIME,NFOR]: FORCING TIME SERIES FOR 'NFOR' FORCINGS
###   F.DATE [FTIME]: DATES OF FORCING TIME SERIES
###   FNAMES[NFOR]: CHARACTER STRING GIVING NAME OF EACH FORCING
###		NFOR:	NUMBER OF FORCINGS.  SET TO ZERO TO IGNORE FORCINGS
###   DEV.LIST: OUTPUT FROM dev.AnyAR.FixNoise.1forcing
###   WHICH.MODEL: EITHER 1 OR 2, INDICATING IF MODEL TO BE SIMULATED WAS 1ST OR 2ND IN DEVIANCE CALL
###   ADD.NOISE: TRUE MEANS ADD RANDOM NOISE; FALSE MEANS DO NOT ADD NOISE
### COMMENTS
###  1) REQUIRES LUBRIDATE PACKAGE
  
p.order = dev.list$p.order
nspace  = dev.list$nspace
nharm   = dev.list$nharm
mperiod = dev.list$mperiod

if (suppress.anncyc) nharm = 0


if (length(initial.date) != 1) stop('initial.date should be length 1') 
if (length(initial.ts) != p.order * nspace) stop('initial incorrectly dimensioned')
if (nfor > 0) {
  ftime = length(f.date)
  if (length(f.ts  ) != ftime * nfor ) stop('f.ts incorrectly dimensioned')
  if (length(fnames) != nfor         ) stop('fnames incorrectly dimensioned')
}

y.date   = initial.date %m+% months(1:nlead-1)

if (which.model == 1) {
  beta.hat  = dev.list$lm1$coefficients
  cov.noise = dev.list$qmat1 / dev.list$lm1$df.residual	  
} else {
  beta.hat  = dev.list$lm2$coefficients
  cov.noise = dev.list$qmat2 / dev.list$lm2$df.residual	   
}


beta.hat    = as.matrix(beta.hat)
beta.names  = rownames(beta.hat)
dim(y.date) = c(nlead)

if (substr(beta.names[1],1,9) == 'xmat.star') prefix = 'xmat.star' else prefix = 'xmat'

if (add.noise) {
	cov.eigen   = eigen(cov.noise)
	noise.trans = cov.eigen$vectors %*% diag(sqrt(cov.eigen$values),nrow=nspace,ncol=nspace)
} else noise.trans = NA

if (nfor > 0) {
  ist   = which(y.date[1]     == f.date)
  ind   = which(y.date[nlead] == f.date)
  if (length(ist) != 1 | length(ind) != 1) stop('cannot align forcing and simulation dates')
  f.ts  = as.matrix(f.ts)
  f.sim = as.matrix(f.ts[ist:ind,])
}


y.sim             = array(NA,dim=c(nlead,nspace))
y.sim[1:p.order,] = initial.ts
phase             = (month(y.date) - 1) %% mperiod + 1
omega             = 2 * pi / 12

for (n in (1+p.order):nlead) {
	y.sim[n,] = beta.hat[paste(prefix,'intercept',sep=''),]
	if (p.order > 0) for (lag in 1:p.order) {
		npic = which( beta.names == paste(prefix,'S1lag',lag,'M',phase[n],sep=''))
		y.sim[n,] = y.sim[n,] + y.sim[n-lag,] %*% beta.hat[npic + 1:nspace-1,]
	}
	if (nharm > 0) {
		nh.max = min(5,nharm)
		for (nh in 1:nh.max) y.sim[n,] = y.sim[n,] + beta.hat[paste(prefix,'cos',nh,sep=''),] * cos(omega * nh * month(y.date[n]))
		for (nh in 1:nh.max) y.sim[n,] = y.sim[n,] + beta.hat[paste(prefix,'sin',nh,sep=''),] * sin(omega * nh * month(y.date[n]))
		if (nharm == 6)      y.sim[n,] = y.sim[n,] + beta.hat[paste(prefix,'cos6'  ,sep=''),] * cos(omega * 6  * month(y.date[n]))
	} 
	if (nfor > 0) for (nf in 1:nfor) y.sim[n,] = y.sim[n,] + f.sim[n,nf] * beta.hat[paste(prefix,fnames[nf],sep=''),]
	if (add.noise)                   y.sim[n,] = y.sim[n,] + noise.trans %*% rnorm(nspace)
}

y.sim

}