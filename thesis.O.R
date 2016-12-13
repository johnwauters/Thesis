start.time = date()
set.seed(2)

#source( "C:/Users/John/Desktop/Code for thesis/thesis.O.R", echo=T)

#O		Using simulated data, fit sets of data with "vanilla" 
#			iterative SIS to a cox model
#P		fit same simulated data using Variant 1 iterative SIS, 
#			compare in terms of mean model size, L2 norm of 
#			difference of true and estimated betas, and frequency 
#			of including true model.

# x = read.table(file=file.choose(), comment.char="!", header=
#	TRUE, row.names=1)

library(SIS)
library(survival)
#library(glmnet)
library(MASS)

### construction of frauda

n = 400		# <<<<<<<< 400
p = 1000		# <<<<<<<< 1000
reps = 100 		# <<<<<<<< 100

make.frauda = function(n = 6, p = 8,
		betastar =  c(
	-1.6328, 1.3988, -1.6497, 1.6353, -1.4209, 1.7022, (rep(0, p-6) )),
	sigma=diag(p), censor.rate=0.1) {

	x = mvrnorm(n=n, mu= rep(0,p), Sigma=sigma)
	const = exp(as.numeric(x %*% matrix(data=betastar,
		 ncol=1)))
	U = runif(n=n)
	times = -log(1-U)/const
	censoring.times = rexp(n=n, rate=censor.rate)
	delta = as.numeric(times > censoring.times)
	y = pmin(times, censoring.times)
	list(x=x, y=y, delta=delta, times=times,
		censoring.times=censoring.times, censoring_rate = 
		mean(delta))
} # end make.frauda function

sigma.case5 = matrix(data=0.5, nrow=p, ncol=p)
	diag(sigma.case5) = rep(1, p)

sigma.case6 = matrix(data=0.5, ncol=p, nrow=p)
	sigma.case6[4,] = rep(sqrt(0.5), p)
	sigma.case6[,4] = rep(sqrt(0.5), p)
	sigma.case6[5,] = rep(0, p)
	sigma.case6[,5] = rep(0, p)
	diag(sigma.case6) = rep(1, p)

beta.case5 = c( -1.5140, 1.2799, -1.5307, 1.5164, -1.3020, 1.5833,
		rep(0, (p-6)) )

beta.case6 = c(4, 4, 4, (-6*sqrt(2)), (4/3), rep(0, (p-5)) )

censor.rate.case5 = 0.045

censor.rate.case6 = 0.03
# case 5 vanilla version records
L1.5 = numeric(reps)
L22.5 = L1.5
P5 = L1.5
MS.5 = P5
#case 6 vanilla version records
L1.6 = P5
L22.6 = P5
P6 = P5
MS.6 = P5
#case 5 variant 1 records
vL1.5 = numeric(reps)
vL22.5 = L1.5
vP5 = L1.5
vMS.5 = P5
#case 6 variant 1 records
vL1.6 = P5
vL22.6 = P5
vP6 = P5
vMS.6 = P5


for (i in 1:reps) {
#case 5
	simulated = make.frauda(n=n, p=p, betastar=beta.case5, 
		censor.rate=censor.rate.case5, sigma=sigma.case5)
	y = Surv(time=simulated$y, event= (1-simulated$delta),
		type='right') 
	x = scale(simulated$x)
# vanilla variant case 5
	van.sis.obj = SIS(x=x, y=y, family='cox', penalty='lasso')
	selected = van.sis.obj$ix
	beta.hat = rep(0, p)
	beta.hat[selected] =  van.sis.obj$coef.est
	beta.diff = beta.hat - beta.case5
	L1.5[i] = sum(abs(beta.diff))
	L22.5[i] = sum(beta.diff^2)
	P5[i] = all(1:6 %in% selected)
	MS.5[i] = length(selected)
# aggressive/variant 1 case 5
	var.sis.obj = SIS(x=x, y=y, family='cox', penalty='lasso', varISIS= 
		'aggr')
	selected = var.sis.obj$ix
	beta.hat = rep(0, p)
	beta.hat[selected] =  var.sis.obj$coef.est
	beta.diff = beta.hat - beta.case5
	vL1.5[i] = sum(abs(beta.diff))
	vL22.5[i] = sum(beta.diff^2)
	vP5[i] = all(1:6 %in% selected)
	vMS.5[i] = length(selected)
# case 6
	simulated = make.frauda(n=n, p=p, betastar=beta.case6, 
		censor.rate=censor.rate.case6, sigma=sigma.case6)
	y = Surv(time=simulated$y, event= (1-simulated$delta),
		type='right') 
	x = scale(simulated$x)
# vanilla variant case 6
	van.sis.obj = SIS(x=x, y=y, family='cox', penalty='lasso')
	selected = van.sis.obj$ix
	beta.hat = rep(0, p)
	beta.hat[selected] =  van.sis.obj$coef.est
	beta.diff = beta.hat - beta.case6
	L1.6[i] = sum(abs(beta.diff))
	L22.6[i] = sum(beta.diff^2)
	P6[i] = all(1:5 %in% selected)
	MS.6[i] = length(selected)
# aggressive/variant 1 case 6
	var.sis.obj = SIS(x=x, y=y, family='cox', penalty='lasso', varISIS= 
		'aggr')
	selected = var.sis.obj$ix
	beta.hat = rep(0, p)
	beta.hat[selected] =  var.sis.obj$coef.est
	beta.diff = beta.hat - beta.case6
	vL1.6[i] = sum(abs(beta.diff))
	vL22.6[i] = sum(beta.diff^2)
	vP6[i] = all(1:5 %in% selected)
	vMS.6[i] = length(selected)
} # end for loop
 # ##############################################  ###

sink("C:/Users/John/Desktop/Code for thesis/thesis.O.txt",
		 append=TRUE, split=TRUE)
cat("\n<><><><><><><><><><>   <><><><><><><><><><>\n")
cat("\n\n\n")

OUT = data.frame(L1.5 = L1.5, L22.5 = L22.5, P5 = P5, MS.5 = MS.5, L1.6 
	= L1.6, L22.6 = L22.6, P6 = P6, MS.6 = MS.6, vL1.5 = vL1.5, vL22.5 
	= vL22.5, vP5 = vP5, vMS.5 = vMS.5, vL1.6 = vL1.6, vL22.6 = 
	vL22.6, vP6 = vP6, vMS.6 = vMS.6)

write.csv(x=OUT, file=
	'C:/Users/John/Desktop/Code for thesis/thesis.O.csv')

median(L1.5)
median(L22.5)
median(P5)
median(MS.5)
median(L1.6)
median(L22.6)
median(P6)
median(MS.6)
median(vL1.5)
median(vL22.5)
median(vP5)
median(vMS.5)
median(vL1.6)
median(vL22.6)
median(vP6)
median(vMS.6)

cat("\nstart time=", start.time, "end time =", date(), "\n")
sink()
#	##############