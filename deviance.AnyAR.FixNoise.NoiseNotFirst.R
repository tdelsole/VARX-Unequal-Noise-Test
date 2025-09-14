deviance.AnyAR.FixNoise.NoiseNotFirst = function(
		ts1,ts2,date1,date2,ntime1,ntime2,nens1,nens2,nspace,
		forcing.ts,forcing.date,ftime,nfor,forcing.names,
		p.order,nharm,niter,mperiod=1,perpetual.pic=NA,perpetual.year=NA,
		alpha=0.05,monte.carlo=FALSE,test.order.call=NA,lag.lbtest=10) {
### COMPUTE EQUALITY OF VARX(P) CYCLOSTATIONARY MODELS
### INPUT:
###		TS1   [NTIME1,NENS1,NSPACE]: TIME SERIES ARRAY 1
###     TS2   [NTIME2,NENS2,NSPACE]: TIME SERIES ARRAY 2 (MUST HAVE SAME SPATIAL DIMENSION AS TS1)
###     DATE1 [NTIME1,NENS1]: DATES FOR TIME SERIES 1
###		DATE2 [NTIME2,NENS2]: DATES FOR TIME SERIES 2
###		NTIME1: LENGTH OF TIME SERIES 1 PER ENSEMBLE MEMBER
###		NTIME2: LENGTH OF TIME SERIES 2 PER ENSEMBLE MEMBER
###		NENS1:  NUMBER OF ENSEMBLE MEMBERS PER LEAD FOR TIME SERIES 1
###		NENS2:  NUMBER OF ENSEMBLE MEMBERS PER LEAD FOR TIME SERIES 2
###		NSPACE: SPATIAL DIMENSION
###		FORCING.TS  [FTIME,NFOR]: FORCING FOR TIME SERIES
### 	FORCING.DATE[FTIME]
###		FTIME: LENGTH OF FORCING TIME SERIES
###		NFOR: NUMBER OF FORCINGS; SET TO 0 TO IGNORE FORCINGS
###     FORCING.NAMES[NFOR]
###     P.ORDER: ORDER OF THE VAR MODEL (can be 0)
###     NHARM: NUMBER OF HARMONICS (can be 0)
###		MPERIOD: PERIOD OF CYCLOSTATIONARY CYCLE (1 FOR STATIONARY [DEFAULT]; 12 FOR MONTHLY CYCLOSTATIONARY)
###		PERPETUAL.PIC[2]: INDICATES TIME SERIES 1 OR 2 TO USE PERPETUAL FORCING.  EG.(1), (2), (1,2). IGNORED IF NA
###     PERPETUAL.YEAR: SINGLE YEAR OF PERPETUAL FORCING.  EG (1850).  
###		ALPHA: SIGNIFICANCE LEVEL
###		MONTE.CARLO: =TRUE TO COMPUTE SIGNIFICANCE THRESHOLDS BY MONTE CARLO, =FALSE TO SKIP
###		TEST.ORDER.CALL: ORDER OF THE TESTS FOR EQUALITY.  DEFAULT: NOISE, AR, CYCLE, FORCING, INTERCEPT
### OUTPUT:
### COMMENTS:
###		1) MONTE CARLO: ONLY ONE CALL NECESSARY IF TIME-SPACE DIMENSION NEVER CHANGE BETWEEN CALLS
###		2) THE FORCING AND DATA TIME SERIES ARE ALIGNED IN TIME AND THUS NEED ONLY ONE TIME ARRAY
###		3) TS1 AND TS2 SHOULD BE TEMPORALLY SEQUENTIAL IN FIRST INDEX; 2ND INDEX IS UNCONSTRAINED
###		4) TO CREATE MONTHLY TIME SERIES OF LENGTH NTOT: date1 = start.date %m+% months(1:ntot-1)
###		5) THE NOISE VARIANCE IS ASSUMED TO BE STATIONARY; ONLY THE AR-COEFFICIENTS DEPEND ON PHASE

alpha.extra   = alpha * c(1,1/5)

######################################################
########## CHECK DIMENSIONS ARE CONSISTENT WITH EXPECTATIONS
######################################################
if (length(ts1)   != ntime1 * nens1 * nspace) stop('ts1 incorrectly dimensioned')
if (length(ts2)   != ntime2 * nens2 * nspace) stop('ts2 incorrectly dimensioned')
if (length(date1) != ntime1 * nens1         ) stop('date1 incorrectly dimensioned')
if (length(date2) != ntime2 * nens2         ) stop('date2 incorrectly dimensioned')
if (nfor > 0 && length(forcing.ts) != ftime * nfor) stop('forcing.ts incorrectly dimensioned')
if (nfor > 0 && nfor               != length(forcing.names)) stop('nfor inconsistent with forcing.names')

##############################################################################
########## CHECK TEST.ORDER INCLUDES ONLY PERMISSIBLE VALUES
#############################################################################
default.order   = c('forcing','ar','cycle','intercept') 
if (any(is.na(test.order.call))) test.order.call = default.order
if (length(setdiff(test.order.call, default.order)) > 0) stop(paste('test.order.call contains unrecognized items; items should be',paste(default.order,collapse=', ')))
if (any(test.order.call == 'noise')) stop('this function does not test equality of noise')
test.order.all = c(test.order.call,setdiff(default.order,test.order.call))



##############################################################################
########## CHECK TEST.ORDER EXCLUDES ITEM WHEN (ORDER, NHARM, NFOR) EQUAL 0
#############################################################################
for (n in 1:length(test.order.call)) {
	if (test.order.call[n] == 'ar'      & p.order == 0) stop('remove ar from test.order.call if p.order = 0')
	if (test.order.call[n] == 'cycle'   & nharm == 0  ) stop('remove cycle from test.order.call if nharm = 0')
	if (test.order.call[n] == 'forcing' & nfor == 0   ) stop('remove forcing from test.order.call if nfor = 0')
}


######################################################
########## SYNCHRONIZE FORCINGS WITH TIME SERIES
### IF PERPETUAL.PIC = 1, THEN F1.TS IS NOT USED SO IT IS SET TO NA
### AND F2.TS = F2.TS - F2.TS[PERPETUAL.YEAR].  ANALOGOUS FOR PERPETUAL.PIC = 2
### CHECK THAT DEVIANCE ON 'FORCING' IS NOT REQUESTED IF PERPETUAL FORCING IS CALLED
######################################################
if (!all(is.na(perpetual.pic)) & any(test.order.call == 'forcing')) stop('cannot test deviance in forcing if perpetual forcing is called')

if (nfor != 0) {
	if (!all(is.na(perpetual.pic))) {
		iyst.persist  = which(as.Date(paste('15-01-',perpetual.year,sep=''),format='%d-%m-%Y') == forcing.date)
		if (length(iyst.persist) != 1) stop('cannot find perpetual year in forcing dates')
	}

	if (!any(perpetual.pic == 1,na.rm=TRUE)) {
		n.forcing.ts1.format = n.to.monthly(date1,forcing.date)
		if (any(date1 != forcing.date[n.forcing.ts1.format])) stop('forcing and ts1 dates misaligned')
		f1.ts = as.matrix(forcing.ts[n.forcing.ts1.format,])
		if (any(perpetual.pic == 2,na.rm=TRUE)) for (nf in 1:nfor) f1.ts[,nf] = f1.ts[,nf] - forcing.ts[iyst.persist + 1:12 - 1, nf]
		nfor1 = nfor
	} else {
		f1.ts = NA
		nfor1 = 0
	}
	
	if (!any(perpetual.pic == 2,na.rm=TRUE)) {
		n.forcing.ts2.format = n.to.monthly(date2,forcing.date)
		if (any(date2 != forcing.date[n.forcing.ts2.format])) stop('forcing and ts2 dates misaligned')
		f2.ts = as.matrix(forcing.ts[n.forcing.ts2.format,])		
		if (any(perpetual.pic == 1,na.rm=TRUE)) for (nf in 1:nfor) f2.ts[,nf] = f2.ts[,nf] - forcing.ts[iyst.persist + 1:12 -1, nf]
		nfor2 = nfor
	} else {
		f2.ts = NA
		nfor2 = 0
	} 
} else {
  nfor1 = 0
  nfor2 = 0
  f1.ts = NA
  f2.ts = NA
}



######################################################
########## SET UP TIME SERIES MATRICES
######################################################
xy1  = timeseries2VARX.cyclo(ts1,date1,f1.ts,ntime1,nens1,nspace,nfor1,forcing.names,p.order,nharm,mperiod)
xy2  = timeseries2VARX.cyclo(ts2,date2,f2.ts,ntime2,nens2,nspace,nfor2,forcing.names,p.order,nharm,mperiod)

######################################################
########## SET UP BREAK POINTS AND DEVIANCE LABELS
######################################################
xmat1         = NULL
xmat2         = NULL
x.break       = NULL

noise.first = FALSE
if (noise.first) {
  dev.row.names = 'D[0:1] Noise'
  step          = 0
} else {
  dev.row.names = NULL
  step          = -1  
}

########################################################
for (n in 1:length(test.order.all)) {
	if (test.order.all[n] == 'ar' & p.order > 0) {
		step          = step+1
		dev.row.names = c(dev.row.names,paste('D[',step,':',step+1,'] AR',sep=''))
		xmat1         = cbind(xmat1,xy1$x.ar)
		xmat2         = cbind(xmat2,xy2$x.ar)
		x.break       = c(x.break,dim(as.matrix(xy1$x.ar))[2])
	}
	if (test.order.all[n] == 'cycle' & nharm > 0) {
		step          = step+1
		dev.row.names = c(dev.row.names,paste('D[',step,':',step+1,'] Cycle',sep=''))	
		xmat1         = cbind(xmat1,xy1$x.cyc)
		xmat2         = cbind(xmat2,xy2$x.cyc)
		x.break       = c(x.break,dim(as.matrix(xy1$x.cyc))[2])
	}
	if (test.order.all[n] == 'forcing' & nfor > 0) {
		step          = step+1
		dev.row.names = c(dev.row.names,paste('D[',step,':',step+1,'] Forcing',sep=''))
		xmat1         = cbind(xmat1,xy1$x.for)
		xmat2         = cbind(xmat2,xy2$x.for)
		x.break       = c(x.break,nfor)
	}
	if (test.order.all[n] == 'intercept') {
		step          = step+1
		dev.row.names = c(dev.row.names,paste('D[',step,':',step+1,'] Intercept',sep=''))		
		xmat1         = cbind(xmat1,xy1$x.int)
		xmat2         = cbind(xmat2,xy2$x.int)
		x.break       = c(x.break,1)
	}
}

#### ELIMINATE TERMS THAT WERE NOT REQUESTED TO BE COMPARED
ntest = length(test.order.all)
ncall = length(test.order.call)
if (ncall < ntest) {
	ncut          = ntest - ncall
	x.break       = x.break[-c(length(x.break) - 1:ncut + 1)]
	dev.row.names = dev.row.names[-c(length(dev.row.names) - 1:ncut + 1)]
} else if (ncall > ntest) stop('test.order contains unrecognized items 2')

dev.row.names = c(dev.row.names,paste('D[0:',length(dev.row.names),'] Total',sep=''))

names(x.break) = test.order.all[1:length(x.break)]
######################################################
########## PERFORM DIFFERENCE-IN-REGRESSION PARAMETER TEST
######################################################
nest.list = diff.regression.nested.mult.NoiseNotFirst(
  xy1$y.lhs,xy2$y.lhs,xmat1,xmat2,x.break = x.break, test.order.call = test.order.call,
  alpha = alpha, alpha.extra = alpha.extra,monte.carlo = monte.carlo, niter = niter)

rownames(nest.list$dev.table)       = dev.row.names
rownames(nest.list$dev.table.extra) = dev.row.names

diagnose.list = 'no diagnosis '
nest.list$diagnose.list = diagnose.list

###### PORTMANTAU TEST
lbtest1 = LjungBox(residuals(nest.list$lm1),lags=lag.lbtest,order=p.order)
lbtest2 = LjungBox(residuals(nest.list$lm2),lags=lag.lbtest,order=p.order)
nest.list$lbtest1 = lbtest1
nest.list$lbtest2 = lbtest2

###### MIC (WHICH REQUIRES REMOVING THE INTERCEPT)
n.int  = which(colnames(xmat1) == 'intercept')
if (length(n.int) != 1) stop('cannot find intercept in xmat1')
if (colnames(xmat1)[n.int] != colnames(xmat2)[n.int]) stop('xmat2 has interecept in unexpected column')
mic1   = mic.master(xmat1[,-n.int],xy1$y.lhs,dim(xmat1)[1],alpha)
mic2   = mic.master(xmat2[,-n.int],xy2$y.lhs,dim(xmat2)[1],alpha)
nest.list$mic1 = mic1
nest.list$mic2 = mic2

######################################################
########## APPEND LIST WITH AUXILIARY DATA
######################################################
nest.list$f1.ts = f1.ts
nest.list$f2.ts = f2.ts

nest.list$nspace  = nspace
nest.list$p.order = p.order
nest.list$nharm   = nharm
nest.list$mperiod = mperiod
#######
            
nest.list


}