timeseries2VARX.cyclo = function(y.ts,y.date,f.ts,nlead,nens,nspace,nfor,forcing.names,order.pic,nharm,mperiod=1) {
#### GIVEN TIMESERIES Y, CREATE MATRICES FOR ARX(P) WITH ANNUAL CYCLE AND FORCING
#### INPUT:
####	Y.TS   [NLEAD,NENS,NSPACE]: TIME SERIES OF LENGTH NLEAD; NENS: NUMBER OF ENSEMBLE MEMBERS
####    Y.DATE [NLEAD,NENS]: DATES ASSOCIATED WITH Y (SEE COMMENTS BELOW FOR FURTHER DETAILS)
####	F.TS   [NLEAD,NENS,NFOR]: FORCING; IGNORED IF F.TS = NA; NFOR = NUMBER OF FORCINGS
####    NLEAD:		NUMBER OF LEADS
#####	NENS: ENSEMBLE SIZE
#####	NSPACE: DIMENSION OF THE STATE SPACE (E.G., NUMBER OF SPACIAL ELEMENTS)
##### 	NFOR: NUMBER OF FORCINGS; FORCINGS ARE IGNORED IF NFOR = 0
#####	FORCING.NAMES[NFOR]: NAME OF THE FORCINGS
#####   ORDER.PIC:	ORDER OF THE AUTOREGRESSIVE MODEL
#####   NHARM:		NUMBER OF ANNUAL HARMONICS
#####	MPERIOD:	PERIOD OF CYCLOSTATIONARY CYCLE (1 FOR STATIONARY [DEFAULT]; 12 FOR MONTHLY CYCLOSTATIONARY)
##### OUTPUT:LIST$
#####	Y.LHS[NTOT]: TIME SERIES AFTER OMITTING FIRST ORDER.PIC STEPS
#####	Y.LAG[NTOT,ORDER.PIC*PERIOD]: LAGGED PREDICTOR MATRIX FOR Y; EQUALS NULL IF ORDER.PIC = 0
#####	Y.CYC[NTOT,2*NHARM]:    COSINE/SINE PREDICTOR MATRIX; EQUALS NULL IF NHARM = 0
#####	Y.FOR[NTOT,NFOR]: FORCING TIME SERIES
#####	X.INT[NTOT]: VECTOR OF ONES (FOR THE INTERCEPT)
##### COMMENTS
##### 1) LIBRARY(LUBRIDATE) IS REQUIRED
##### 2) THE ANNUAL HARMONIC IS ALWAYS PERIOD 12.  IT DOES NOT NECESSARILY EQUAL MPERIOD (the cyclostationary period)
##### 3) FOR MONTHLY DATA, y.date = start.date %m+% months(1:nlead-1); e.g., start.date = as.Date('15-01-1900',format='%d-%m-%Y')

if (length(y.ts) != nlead * nens * nspace) stop('y.ts incorrectly dimensioned in timeseries2Varx.cyclo')

x.names     = NULL
dim(y.ts)   = c(nlead*nens,nspace)
dim(y.date) = c(nlead*nens)

if (nfor != 0) dim(f.ts) = c(nlead*nens,nfor)

### THIS PHASE REFERS TO CYCLOSTATIONARY PHASE, NOT PHASE OF ANNUAL CYCLE
phase          = (month(y.date) - 1) %% mperiod + 1  ## GENERATES INTEGERS 1,2, ..., MPERIOD
lead.all       = rep(1:nlead,nens)
num.per.phase  = as.numeric(rep(NA,mperiod))
num.predictors = as.numeric(rep(0 ,mperiod))

y.lhs         = NULL
for (m in 1:mperiod) {
	npic    = which( phase == m & lead.all > order.pic)
	y.lhs   = rbind(y.lhs,as.matrix(y.ts[npic,]))
	num.per.phase[m] = length(npic)
}

x.ar        = NULL
names.ar    = NULL
space.say   = paste('S',1:nspace,sep='')
if (order.pic > 0) for (m in 1:mperiod) {
	npic    = which( phase == m & lead.all > order.pic)
	x.all.lag = NULL
	for (lag in 1:order.pic) x.all.lag = cbind(x.all.lag,y.ts[npic-lag,])
	
	if (m == 1      ) zero.pre = NULL else zero.pre = array(0,dim=c(length(npic),nspace*order.pic*(m-1)))
	if (m == mperiod) zero.suf = NULL else zero.suf = array(0,dim=c(length(npic),nspace*order.pic*(mperiod-m)))
	x.ar = rbind(x.ar,cbind(zero.pre,x.all.lag,zero.suf))

	for (lag in 1:order.pic) names.ar = c(names.ar,paste(space.say,'lag',lag,'M',m,sep=''))
	num.predictors[m] = num.predictors[m] + dim(as.matrix(x.all.lag))[2]
}

x.cyc     = NULL
names.cyc = NULL
if (nharm > 0) for (m in 1:mperiod) {
	npic    = which(phase == m & lead.all > order.pic)
	x.temp  = NULL
	omega.t = 2 * pi * month(y.date[npic]) / 12
	for (nh in 1:nharm) x.temp = cbind(x.temp,cos(omega.t * nh),sin(omega.t * nh))
	if (nharm == 6) x.temp = x.temp[,-12]
	x.cyc = rbind(x.cyc,x.temp)
	if (m == 1) {
		for (nh in 1:nharm) names.cyc = c(names.cyc,paste('cos',nh,sep=''),paste('sin',nh,sep=''))
		if (nharm == 6) names.cyc = names.cyc[-12]
	}
	num.predictors[m] = num.predictors[m] + dim(as.matrix(x.temp))[2]
}

x.for     = NULL
if (nfor > 0) for (m in 1:mperiod) {
	npic  = which(phase == m & lead.all > order.pic)
	x.for = rbind(x.for,f.ts[npic,,drop=FALSE])
	num.predictors[m] = num.predictors[m] + dim(as.matrix(x.for))[2]
}

x.int = as.matrix(rep(1,dim(y.lhs)[1]))
num.predictors = num.predictors + 1

if (order.pic > 0) colnames(x.ar)  = names.ar
if (nharm     > 0) colnames(x.cyc) = names.cyc
if (nfor      > 0) colnames(x.for) = forcing.names
colnames(x.int)  = 'intercept'


list(y.lhs=y.lhs,x.ar=x.ar,x.cyc=x.cyc,x.for=x.for,x.int=x.int,num.per.phase=num.per.phase,num.predictors = num.predictors)

}