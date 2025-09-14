mic.master = function(x,y,nsamp,alpha=0.05) {
########################################
### COMPUTE MIC FOR Y = XB + MU + E
### INPUT
###   X[NSAMP,XDIM]: X-ARRAY OF PREDICTORS
###   Y[NSAMP,YDIM]: Y-ARRAY OF PREDICTANDS
###   APLPHA: SIGNIFICANCE LEVEL FOR TESTING B = 0 
### OUTPUT:
###   MIC
###   MIC.CRIT: THE ALPHA * 100% SIGNIFICANCE LEVEL FOR B = 0
###   PVAL: P-VALUE OF MIC FOR TESTING B = 0
########################################
if (length(y) %% nsamp != 0) stop('y dimensioned incorrectly')
if (length(x) %% nsamp != 0) stop('x dimensioned incorrectly')
ydim = length(y) / nsamp
xdim = length(x) / nsamp

if (nsamp-xdim-ydim-2 <= 1) return(list(mic=NA,mic.crit=NA,pval=NA))

dim(x) = c(nsamp,xdim)
dim(y) = c(nsamp,ydim)

x     = t(t(x)-colMeans(x))
y     = t(t(y)-colMeans(y))

x.svd = svd(x)
y.svd = svd(y)

#################################################
## CHECK FOR ZERO SINGULAR VALUE AND REMOVE IT
## TO ACCOUNT FOR INTERCEPT IN X-MATRIX
#################################################
if (x.svd$d[xdim]/x.svd$d[1] < 1.e-10) {
	x.svd$d = x.svd$d[ -xdim]
	x.svd$u = x.svd$u[,-xdim]
	x.svd$v = x.svd$v[,-xdim]
	xdim    = xdim - 1
}

if (x.svd$d[xdim]/x.svd$d[1] < 1.e-10) stop('X is not full rank')
if (y.svd$d[ydim]/y.svd$d[1] < 1.e-10) stop('Y is not full rank')

xy.svd = svd(t(x.svd$u) %*% y.svd$u)

penalty  = (nsamp+1) * ( (xdim+ydim)/(nsamp-(xdim+ydim)-2) - xdim/(nsamp-xdim-2) - ydim/(nsamp-ydim-2) )

if (any(xy.svd$d >= 1)) {
	mic      = NA
	mic.crit = NA
	pval     = NA
} else {
	lambda   = sum(log(1-xy.svd$d^2))
	mic      = lambda + penalty
	
	nu.e     = nsamp - xdim - 1
	nu.h     = xdim
	fctr     = nu.e-(ydim-nu.h+1)/2
	mic.crit = -qchisq(alpha,df=ydim*xdim,lower.tail=FALSE)/fctr + penalty
	pval     = pchisq(-(mic-penalty)*fctr,df=ydim*xdim,lower.tail=FALSE)
	
	mic      = mic * nsamp
	mic.crit = mic.crit * nsamp	
}

list(mic = mic, mic.crit = mic.crit, pval = pval, penalty = penalty)

}