# "C:/Users/John/Desktop/Code for thesis/SIS_data_simulation.R"
# source( , echo=T)
# note to self: Georgia font 14 point

## url for neuroblastoma data
#	ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24080/
#	matrix/GSE24080_series_matrix.txt.gz

cat("\nDate/time  ..............", date(), "\n")

# simulating the data for SIS paper <<<<<<<<<<
#		Fan and Lv (2008)

# simulationI: independent features

# ######## define lqnorm function
	# finds the Lq norm of a numeric vector, column matrix, or
	# row matrix, q defaults to 2
lqnorm = function(x, q=2) {
	if (q<0) {cat("inproper q"); return(numeric(0))}
	x = abs(as.numeric(x))
	if (length(x) < 1) {cat("improper x"); return(numeric(0))}
	if (q==0) return(sum(x != 0))
	if (q==Inf) return(max(x))
	return(  ((sum((x**q)))^(1/q))  )
}  # ## end of function


### reconstruct.long.beta(logical.vector, short.beta) is a function
### for turning short
### sparse beta vectors into sparse vectors the size of the
### true beta for purposes of comparing the two vectors
reconstruct.long.beta = function(logical.vector, short.beta) {
	if (is.null(logical.vector) || is.null(short.beta) || length(logical.vector)
		== 0 || length(short.beta) == 0 || !is.logical(logical.vector) || 
		!is.numeric(short.beta)) return(NULL)
	short.beta = as.numeric(short.beta)
 	if (length(short.beta) != sum(logical.vector)) {cat("\nincompatable 
		inputs to reconstruct.long.beta\n"); return(NULL) }
	output.beta = numeric(0)
	while (length(logical.vector > 0) )  {
		if (!logical.vector[1]) {output.beta = c(output.beta, 0)}
		else {output.beta = c(output.beta, short.beta[1]); 
			short.beta = short.beta[-1]}
		logical.vector = logical.vector[-1]
	} # end while loop
	return(output.beta)	
} # end reconstruct.long.beta.function   ##################

### SIS as a function that outputs a sparse logical vector
### whose "TRUE" entries correspond to important features
## x is the n by p data matrix, y is the response column matrix
## d is the number of important variables to be selected
SIS = function(x, y, d) {
	x.reg = scale(x)
	y.reg = scale(y)
	omeg = t(x.reg)%*%y.reg
	abs.omeg = abs(omeg)
	if (d >= length(omeg)) return(rep(TRUE, length(omeg)))
	threshold = sort(abs.omeg, decreasing=TRUE)[d]
	return((abs.omeg >= threshold))
} ### end of SIS function

#### CODE FOR BIC-SCAD FOLLOWS: #############
## based on "Tunig parameter selectors for smoothly clipped 
## absolute deviation method";  Wang, Hanshen; Li, Runze; and
## Tsai, Chih-Ling. (Date of publication unknown)

library(ncvreg);

p.prime.lambda = function(theta, lambda, a=3.7){
	return( lambda*((lambda >= theta) + ((a*lambda-
	theta)*(a*lambda >= theta)*(theta > lambda)/(lambda*(a-1)))))
} ## end of p.prime.lambda function

reduce.matrix = function(x, beta, min.step=.Machine$double.eps) {
	beta = abs(beta)
	beta = beta > min.step
	x = x[,beta]
	return(x)
} ## end of reduce.matrix function

reduce.vector = function(beta1, beta2=beta1, min.step=.Machine$double.eps) {
	beta2 = abs(beta2)
	beta2 = beta2 > min.step
	beta1 = beta1[beta2]
	return(beta1)
} ## end of reduce.vector function

sigma.lambda = function(beta, lambda) {
	beta=abs(beta)
	step1 = p.prime.lambda(theta=beta, lambda=lambda)
	step2 = step1/beta
	step3 = length(step2)
	if (step3 < 2) return(matrix(data=1, nrow=1, ncol=1))
	else {  return( diag(step2 ) )} # end else
	
} # end of sigma.lambda function

df.lambda = function(x, beta, n , lambda) {
	x = as.matrix(x)
	if (length(x) < 1) return(0)
	if (length(beta) < 1) return(0)
	if (length(beta) != dim(x)[2]) return (0)
	x.red = reduce.matrix(x=x, beta=beta)
	beta.red = reduce.vector(beta1=beta)
	red.sigma.lambda = reduce.matrix(x=sigma.lambda(beta=
		beta, lambda=lambda), beta=beta)
	bob =
	sum(
		diag(
			x.red %*% solve(
				t(x.red) %*% x.red + (n*red.sigma.lambda)
			) %*% t(x.red)
		)
	)
	return(bob)
} ## end of df.lambda function

sig.hat.sqr.lamb = function(x,  y, beta, n, lambda) { 
	if(length(x) < 1) return(lqnorm(x=y)/n)
	else return(
		(1/n)*(	
			lqnorm(
				(y - x%*%beta)
			)
		)**2
	)
} ## end of sig.hat.sqr.lamb function

my.bic.scad = function(x,  y, lambda=3.7*(999:1)/1000, truebeta) {
	x.des = cbind(rep(1,dim(x)[1]), x)
	n = nrow(x)
	fit.beta <- ncvreg(X=x, y=y, penalty="SCAD", max.iter =
		max( ncol(x), 1000), lambda=lambda, gamma=3.7)$beta
	BIC.vals = ncol(fit.beta)
	for (i in 1:ncol(fit.beta)) {
		BIC.vals[i] = 
		log(sig.hat.sqr.lamb
			(x=x.des, y=y, beta = 
				fit.beta[,i], n=n, lambda=as.numeric(colnames(
				fit.beta)[i]) 
			)
		) + 
		df.lambda(
			x=reduce.matrix(x=x.des, beta=fit.beta[,i]), 
			beta=reduce.vector(beta1 = fit.beta[,i]),
			n=n, 
			lambda=as.numeric(colnames(fit.beta))[i] 
		) * 
			log(n)/n
	} # end of for loop
	BIC.vals
	opt.bic = min(BIC.vals)
	opt.bic.position = 1
	while (opt.bic != BIC.vals[opt.bic.position]) {opt.bic.position = 	opt.bic.position + 1}

	return(fit.beta[,opt.bic.position])
#	model.size =  lqnorm(fit.beta[,opt.bic.position] , q=0)
#	estimation.error = lqnorm(fit.beta[,opt.bic.position] - c(0,
#		truebeta))
#	return(list(model.size=model.size, estimation.error =
#		estimation.error))
} ## end of my.bic.scad function
##############################################

#### PARAMETERS FOR SIMULATION ###
#simIa independent features, p=1000
n = 200 # <<<< 200
p = 1000 # <<<<< 1000
s = 8
num.sims = 500 # <<<<< 500

a = 4*log(n)/sqrt(n)

model.size.ds=numeric(num.sims)
estimation.error.ds = model.size.ds
model.size.lasso = model.size.ds
estimation.error.lasso = model.size.ds
betastar.L2norm = model.size.ds
model.size.sis.scad = model.size.ds
estimation.error.sis.scad = model.size.ds
model.size.sis.ds = model.size.ds
estimation.error.sis.ds = model.size.ds
model.size.sis.ds.scad =model.size.ds
estimation.error.sis.ds.scad = model.size.ds




set.seed(72)
for (i in 1:num.sims) {
## GENERATE THE NON-ZERO ELEMENTS OF betastar, THE
## TRUE BETA FOR THE MODEL
betastar = ((-1)**rbinom(n=s, size=1, prob=0.4)) * 
	(a + abs(rnorm(n=s, mean=0, sd=1)))
betastar.L2norm[i] = lqnorm(betastar)
### GENERATE THE FIRST s  COLUMNS OF x
x = matrix(data = rnorm(n=n*s, mean=0, sd=1), nrow = n, ncol =   
  	s, byrow = FALSE, dimnames = NULL)
### GENERATE y = x%*%beta + noise
y = matrix(ncol=1, data=rnorm(n=n, mean=0, sd=1.5)) + (x %*%
	 betastar)
###  CENTER y AT ZERO
y = y - mean(y)
### FILL IN REMAINING ELEMENTS OF betastar
betastar = c(betastar, rep(0, (p-s)) )
### FILL IN REMAINING ELEMENTS OF x
x = cbind(x, matrix(data= rnorm(n=n*(p-s), mean=0, sd=1),
					 nrow=n, ncol=(p-s)))
dimnames(x) = list(NULL, 1:p)

### DANTZIG SELECTOR SECTION #####
#library(fastclime) ## for dantzig selector function
#ds.obj = dantzig(X=x, y=y, lambda = sqrt(2*log(p)), nlambda = 2*p)
#ds.out = ds.obj$BETA0[,ds.obj$validn]
#model.size.ds[i] = sum(ds.out != 0)
#estimation.error.ds[i] = lqnorm(ds.out - betastar)

### LARS LASSO SECTION #############
library(lars)		## for LARS LASSO function
lars.obj = lars(x=x, y=y, type="lasso", normalize=TRUE 
	, use.Gram=FALSE, max.steps = floor(sqrt(4*p))
	)
lars.beta = lars.obj$beta[length(lars.obj$df),]
#cat("lars beta", lars.beta, "<<<<\n")
model.size.lasso[i] = sum(lars.beta != 0)
estimation.error.lasso[i] = lqnorm(lars.beta-betastar)

### SIS called #################################
sis.vec = SIS(x=x, y=y, d=floor(n-1))
sparse.x = x[,sis.vec]

### MY.BIC.SCAD FUNCTION CALLED ###############
out.bic.scad = my.bic.scad(x=sparse.x,  y = y, truebeta =
	betastar[sis.vec])
long.beta = reconstruct.long.beta(logical.vector=c(TRUE, sis.vec),
	 short.beta=out.bic.scad)
estimation.error.sis.scad[i] = lqnorm(long.beta - c(0, betastar))
model.size.sis.scad[i] = sum(abs(out.bic.scad) > sqrt(
	.Machine$double.eps) )


#out.bic.scad[[1]]

#### SIS-DS SECTION ########################

library(fastclime) ## for dantzig selector function
modified.p = 0.1*ncol(sparse.x)
ds.obj = dantzig(X=sparse.x, y=y, lambda = sqrt(2*log(modified.p)),
	nlambda = 2*modified.p)
ds.out = ds.obj$BETA0[,ds.obj$validn]
long.beta = reconstruct.long.beta(logical.vector=sis.vec,
	 short.beta=ds.out)

model.size.sis.ds[i] = sum(abs(ds.out) > sqrt(.Machine$double.eps))
estimation.error.sis.ds[i] = lqnorm(long.beta - betastar)

### SIS-DS-SCAD Section #######################

intermediate.vec = (long.beta != 0)
if (sum(intermediate.vec) == 0) { model.size.sis.ds.scad[i] = 0
	estimation.error.sis.ds.scad[i] = lqnorm(betastar)
	}
else {
	intermediate.x = as.matrix(x[ , intermediate.vec])

	

	out.bic.scad2 = my.bic.scad(x=intermediate.x,  y = y, truebeta
		 =(betastar[intermediate.vec]))
	sis.ds.scad.beta = reconstruct.long.beta(logical.vector=c(TRUE,
		 intermediate.vec), short.beta=out.bic.scad2)

	model.size.sis.ds.scad[i] = sum(abs(sis.ds.scad.beta) > sqrt(
		.Machine$double.eps))
	estimation.error.sis.ds.scad[i] = lqnorm(sis.ds.scad.beta - c(
		0, betastar))
} # end of else


} # end for i= loop 
#### END OF RUN SUMMARIES ###############
summary(betastar.L2norm)
cat("n=", n, "p=", p, "number of simulations =", num.sims)
summary(model.size.ds)
summary(estimation.error.ds)
summary(model.size.lasso)

summary(estimation.error.lasso)
summary(model.size.sis.scad)
summary(estimation.error.sis.scad)
summary(model.size.sis.ds)
summary(estimation.error.sis.ds)
summary(model.size.sis.ds.scad)
summary(estimation.error.sis.ds.scad)

date()
cat("\n")
# ###