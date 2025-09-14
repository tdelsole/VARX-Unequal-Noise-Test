diff.regression.nested.mult.NoiseNotFirst = 
  function(yvec,yvec.star,xmat,xmat.star,x.break,test.order.call,
           alpha=0.05,monte.carlo=FALSE,ntrials=1000,
           alpha.extra=c(0.05,0.01),niter=10) {
###### TEST EQUALITY OF REGRESSION MODELS THROUGH NESTED HYPOTHESES:
###### YVEC = X2  B2  + ... + XL  BL  + E
###### YVEC = X2* B2* + ... + XL* BL* + E*
###### IMPORTANT: USER SHOULD INCLUDE INTERCEPT IN XMAT AND XMAT.STAR
###### HYPOTHESES: 
######	OMEGA.0: UNRESTRICTED
######  OMEGA.1: EQUALITY OF NOISE VARIANCE (VAR[E] = VAR[E*]):  X.COMM[1] = 0
######  OMEGA.2: OMEGA.1 AND B2 = B2*; LENGTH(B2) = X.BREAK[1] = X.COMM[2]
######  OMEGA.3: OMEGA.2 AND B3 = B3*; LENGTH(B3) = X.BREAK[2] = X.COMM[3]
######  ...
###### 	OMEGA.L: OMEGA.L-1 AND BL = BL*; TESTED ONLY IF X.BREAK[L-1] IS SPECIFIED
###### INPUT:
###### 	YVEC     [NTOT,NSPACE]: 	Y      TIMESERIES
###### 	YVEC.STAR[NTOT.STAR,NSPACE]:Y.STAR TIMESERIES
######	XMAT     [NTOT,KTOT]:		MATRIX [X2 ,X3 ,...,XL ]
######	XMAT.STAR[NTOT,KTOT.STAR]:	MATRIX [X2*,X3*,...,XL*]
######  X.BREAK  [LTOT]:			NUMBER OF COLUMNS IN X2, X3, ... XL; SUM(X.BREAK) <= KTOT
######	ALPHA:						SIGNIFICANCE LEVEL; DEFAULT IS 5%
######  MONTE.CARLO[LOGICAL]:		SET TO 'TRUE' TO DO MONTE CARLO ESTIMATES OF CRITICAL VALUES (EXPENSIVE!)
######	NTRIALS:					NUMBER OF MONTE CARLO TRIALS FOR ESTIMATING DEVIANCE [OMEGA.0 - OMEGA.1]
######  ALPHA.EXTRA[NALPHA]:		ADDITIONAL SIGNIFICANCE LEVELS TO ASSESS; DEFAULT IS 5% AND 1%
######  NITER: NUMBER OF ITERATIONS FOR SOLUTION
###### OUTPUT:
######	DEV.TABLE[LTOT+1,7]: ROWS: LTOT SUB-DEVIANCES, PLUS TOTAL DEVIANCE.  COLUMNS: DEV, CRITICAL (ASYM), PVAL (ASYM), ALPHA, CRITICAL (MC)
######	LM.OMEGA[[LTOT]]: OUTPUT CAPTURED FROM 'LM' FOR HYPOTHESES OMEGA_1,...,OMEGA_LTOT
###### 	LM1: OUTPUT CAPTURED FROM 'LM' FOR FITTING MODEL 1 ALONE
###### 	LM2: OUTPUT CAPTURED FROM 'LM' FOR FITTING MODEL 2 ALONE
###### 	ALPHA.STEP: THE INDIVIDUAL ALPHA LEVEL TO CONTROL FAMILYWISE ERROR RATE
###### 	DEV.TABLE.EXTRA[LTOT+1,4]: SIGNIFICANCE LEVELS FOR LTOT SUB-DEVIANCES, PLUS TOTAL DEVIANCE, FOR ALPHA AND ALPHA/5, ASYMPTOTIC AND MC.   
###### 	QMAT1, QMAT2, QMAT.OMEGA, XMAT.OMEGA1: ARRAYS USED IN DIAGNOSIS
###### 	CDA.NOISE.RATIO.MAX: ALPHA/2 SIGNIFICANCE THRESHOLD FOR MAXIMUM DISCRIMINANT NOISE VARIANCE RATIO 
###### 	CDA.NOISE.RATIO.MIN: ALPHA/2 SIGNIFICANCE THRESHOLD FOR MAXIMUM DISCRIMINANT NOISE VARIANCE RATIO 
###### 	DEV.CRIT.PARAMDIFF[LTOT]: SIGNIFICANCE THRESHOLD FOR MAXIMUM DISCRIMINANT FOR D[I:I+1].  I=1 IS NA BECAUSE ONLY RATIO IS REPORTED
###### 	CDA.ALL[[LTOT]]: LIST OF (EIGENVALUES, DEVIANCES, Q, P) FROM CDA OF D[I:I+1]
###### 	KTOT: TOTAL NUMBER OF PREDICTORS IN FIRST MODEL
###### 	KTOT.STAR: TOTAL NUMBER OF PREDICTORS IN 2ND MODEL
###### 	NCOMM: NUMBER OF PREDICTOR BLOCKS TO BE TESTED FOR EQUALITY ("KAPPA")
###### 	DOF1: DEGREES OF FREEDOM FOR FITTING 1ST MODEL
###### 	DOF2: DEGREES OF FREEDOM FOR FITTING 2ND MODEL
###### 	DOF.SUM: DOF1 + DOF2


##############################
###### METADATA
##############################
# x.comm    = c(0,x.break) # append 0 at the start, corresponding to no common predictors
x.comm    = x.break
xmat      = as.matrix(xmat)
xmat.star = as.matrix(xmat.star)
yvec      = as.matrix(yvec)
yvec.star = as.matrix(yvec.star)

ncomm     = length(x.comm)
ntot      = dim(yvec)[1]
ntot.star = dim(yvec.star)[1]
nspace    = dim(yvec)[2]
if (sum(x.comm) > dim(xmat)     [2]) stop('x.comm inconsistent with xmat')
if (sum(x.comm) > dim(xmat.star)[2]) stop('x.comm inconsistent with xmat.star')

alpha.step  = 1 - (1-alpha)^(1/ncomm)
param.num0  = NA

nalpha           = length(alpha.extra)
alpha.extra.step = 1 - (1-alpha.extra)^(1/ncomm)

nstep      = 1:ntot
nstep.star = 1:ntot.star + ntot

########################################
####### EQUALITY OF VARIANCES (OMEGA_1)
########################################
lm1     = lm(yvec      ~ xmat       - 1)
lm2     = lm(yvec.star ~ xmat.star  - 1)
qmat1   = t(residuals(lm1)) %*% residuals(lm1) 
qmat2   = t(residuals(lm2)) %*% residuals(lm2) 
dof1    = lm1$df.residual
dof2    = lm2$df.residual
dof.sum = dof1 + dof2

param.num0 = nspace * (ncol(xmat) + ncol(xmat.star)) + nspace*(nspace+1)

loglik0    = 
  dof1 * determinant(as.matrix(qmat1/dof1))$modulus + 
  dof2 * determinant(as.matrix(qmat2/dof2))$modulus

########################################
####### COMPUTE COVARIANCE MATRICES FOR NESTED HYPOTHESES
########################################
y.all          = rbind(yvec,yvec.star)
qmat.omega     = array(NA,dim=c(nspace,nspace,ncomm))
dof.omega      = as.numeric(rep(NA,ncomm))
dev.omega      = as.numeric(rep(NA,ncomm))
param.num      = as.numeric(rep(NA,ncomm))
phi            = as.numeric(rep(NA,ncomm))
pval.asym      = as.numeric(rep(NA,ncomm))
loglik.all     = array(NA,dim=c(niter+1,ncomm))
gamma.all      = array(NA,dim=c(nspace,nspace,ncomm))
gamma.all.star = array(NA,dim=c(nspace,nspace,ncomm))
param.diff     = as.numeric(rep(NA,ncomm))
dev.crit.asym  = as.numeric(rep(NA,ncomm))
dev.crit.extra.asym = array(NA,dim=c(nalpha,ncomm))
dev.omega.check     = as.numeric(rep(NA,ncomm))

for (nb in 1:ncomm) {
	ncol.comm      = sum(x.comm[1:nb])
	ncol.diff      = ncol(xmat)      - ncol.comm
	ncol.diff.star = ncol(xmat.star) - ncol.comm
	
	if (ncol.comm == 0) {
		x.all = NULL
	} else {
		x.all  = rbind(xmat[,1:ncol.comm,drop=FALSE],xmat.star[,1:ncol.comm,drop=FALSE])
	}
	
	if (ncol.diff > 0) {
		x.diff      = cbind(xmat     [,1:ncol.diff      + ncol.comm])
		x.all       = cbind(x.all,rbind(x.diff,array(0,dim=c(ntot.star,ncol.diff))))
	}
	if (ncol.diff.star > 0) {
		x.diff.star = cbind(xmat.star[,1:ncol.diff.star + ncol.comm])		
		x.all       = cbind(x.all,rbind(array(0,dim=c(ntot,ncol.diff.star)),x.diff.star))
	}
	
	lm.nest          = lm(y.all ~ x.all - 1)
	qmat.omega[,,nb] = t(residuals(lm.nest)) %*% residuals(lm.nest)
	dof.omega[nb]    = lm.nest$df.residual
	if (nb == 1) lm.omega = list(lm.nest) else lm.omega[nb] = list(lm.nest)
	param.num[nb]    = nspace * ncol(x.all) + nspace*(nspace+1)
	phi[nb]          = ncol(x.all)
	weight           = ntot      / (ntot + ntot.star)
	weight.star      = ntot.star / (ntot + ntot.star)
	loglik.all[1,nb] = (dof1 + dof2) * determinant(as.matrix(qmat.omega[,,nb]/(dof1+dof2)))$modulus
	
	##################################################
	########## BEGIN ITERATIVE SOLUTION
	##################################################
	
	big.x      = x.all[nstep     ,]
	big.x.star = x.all[nstep.star,]
	beta.all   = lm.nest$coefficients
	beta.rownames = rownames(beta.all)
	dim.beta   = dim(beta.all)
	for (ni in 1:niter) {
	  err        = y.all - x.all %*% beta.all
	  gamma      = t(err[nstep     ,]) %*% err[nstep     ,] / dof1
	  gamma.star = t(err[nstep.star,]) %*% err[nstep.star,] / dof2
	  gamma.inv  = chol2inv(chol(gamma))
	  gamma.inv.star = chol2inv(chol(gamma.star))
	  lhs        = kronecker( gamma.inv, t(big.x) %*% big.x) + kronecker( gamma.inv.star , t(big.x.star) %*% big.x.star)
	  rhs        = t(big.x) %*% y.all[nstep,] %*% gamma.inv + t(big.x.star) %*% y.all[nstep.star,] %*% gamma.inv.star
	  beta.all   = chol2inv(chol(lhs)) %*% as.numeric(rhs)
	  dim(beta.all) = dim.beta
	  loglik.all[ni+1,nb] = 
	    dof1 * determinant(as.matrix(gamma     ))$modulus + 
	    dof2 * determinant(as.matrix(gamma.star))$modulus
	  ### CHECK
 	  # check = (t(big.x) %*% big.x) %*% beta.all %*% gamma.inv + 
 	  #   (t(big.x.star) %*% big.x.star) %*% beta.all %*% gamma.inv.star
#     print(max(abs(check-rhs)/abs(rhs)))
	}
	gamma.all     [,,nb] = gamma
	gamma.all.star[,,nb] = gamma.star
	rownames(beta.all)   = beta.rownames
	model.list           = list(beta = beta.all,gamma = gamma,gamma.star = gamma.star)
	
	if (nb == 1) {
	  gamma.pre = qmat1/dof1; gamma.pre.star = qmat2/dof2; param.pre = param.num0
	  model.omega = list(model.list)
	} else {
	  gamma.pre = gamma.all[,,nb-1]; gamma.pre.star = gamma.all.star[,,nb-1]; param.pre = param.num[nb-1]
	  model.omega[nb] = list(model.list)
	}
	gev1 = gev(gamma.all     [,,nb],gamma.pre     )
	gev2 = gev(gamma.all.star[,,nb],gamma.pre.star)
	dev.omega[nb] = dof1 * sum(log(gev1$lambda)) + dof2 * sum(log(gev2$lambda))
	
	#### CHECK DEVIANCE
	dev.omega.check[nb] = 
	  dof1 * (determinant(as.matrix(gamma.all     [,,nb]))$modulus - determinant(as.matrix(gamma.pre     ))$modulus) +
	  dof2 * (determinant(as.matrix(gamma.all.star[,,nb]))$modulus - determinant(as.matrix(gamma.pre.star))$modulus)
	
	dev.crit.asym      [ nb] = qchisq(alpha.step      ,param.pre - param.num[nb],lower.tail=FALSE)
	dev.crit.extra.asym[,nb] = qchisq(alpha.extra.step,param.pre - param.num[nb],lower.tail=FALSE)
	pval.asym          [ nb] = pchisq(dev.omega[nb]   ,param.pre - param.num[nb],lower.tail=FALSE)
	param.diff         [ nb] = param.pre - param.num[nb]
}



param.diff.total = param.num0 - param.num[ncomm]

#### TOTAL DEVIANCE
dev.total                 = sum(dev.omega)
dev.total.pval.asym       = pchisq(dev.total  ,param.num0 - param.num[ncomm],lower.tail=FALSE)
dev.total.crit.asym       = qchisq(alpha      ,param.num0 - param.num[ncomm],lower.tail=FALSE)
dev.total.crit.extra.asym = qchisq(alpha.extra,param.num0 - param.num[ncomm],lower.tail=FALSE)


########################################
#### SUMMARY TABLE
########################################
dev.table = cbind(dev.omega,dev.crit.asym,pval.asym,alpha.step,param.diff)
dev.table = rbind(dev.table,c(dev.total,dev.total.crit.asym,dev.total.pval.asym,alpha,param.diff.total))
rownames(dev.table) = c(paste('D[',1:ncomm-1,':',1:ncomm,']',sep=''),paste('D[0:',ncomm,']',sep=''))
colnames(dev.table) = c('deviance','crit(chisq)','pval(chisq)','alpha','dof')

if (length(names(x.break)) > 0) {
	r.names = c(paste('D[',1:ncomm-1,':',1:ncomm,']',sep=''),paste('D[0:',ncomm,']',sep=''))
	r.names = paste(r.names,c(names(x.break),'total'))
	rownames(dev.table) = r.names
}

dev.table.extra = cbind(t(cbind(dev.crit.extra.asym,dev.total.crit.extra.asym)))
rownames(dev.table.extra) = rownames(dev.table)
colnames(dev.table.extra) = c(paste('asym',(1-alpha.extra)*100,'%',sep=''))

########################################
#### CLEAN UP LOG-LIKELIHOODS
########################################
dev.eq.noise           = loglik.all - loglik0
rownames(loglik.all)   = paste('iter',0:niter)
colnames(loglik.all)   = test.order.call
rownames(dev.eq.noise) = paste('iter',0:niter)
colnames(dev.eq.noise) = test.order.call

########################################
#### OUTPUT THE RESULTS
########################################
list(dev.table = dev.table, lm.omega = lm.omega, lm1=lm1, lm2=lm2, alpha=alpha,
  alpha.step=alpha.step, dev.table.extra=dev.table.extra,
  qmat1 = qmat1, qmat2 = qmat2, qmat.omega = qmat.omega, 
  ncomm = ncomm, dof1 = dof1, dof2 = dof2, dof.sum = dof.sum,
  param.num0 = param.num0, param.num = param.num,
  dev.eq.noise = dev.eq.noise, 
  model.omega = model.omega, x.break = x.break)

}

### ELIMINATE EXTRANEOUS OUTPUT LIST
### COMPARE WITH OLD DEVIANCE