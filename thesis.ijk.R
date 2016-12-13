##### part i: (an j an k?) comparison sis, isis, and lasso at selecting
# true model, simulated data (4.2.3 exampleIII) p=1000, n=20,50,70,
# rho=.5, 200 replications
# source("C:/Users/John/Desktop/Code for thesis/thesis.i.R", echo=T)
 
#sink(file="C:/Users/John/Desktop/Code for thesis/output.i.txt", 

#split=T, append=  T )

p = 100 # <<<<< 100, 1000
n = 20 # <<<<< 20, 50, 70
rho = 0.5 # <<<<<  0.5
sqr.rho = sqrt(rho) 
replications = 200  # <<<<< 200
set.seed(1984)
library(lars)
library(SIS) # <<<<<<<<
report.lasso = logical(replications)
report.sis = report.lasso
report.i.sis = report.sis
report.alternate = report.sis
check.correlation = numeric(replications)

plain.simple.lasso.coef = function(x, y) {
	A=x
	B=y
	 #First get cross validation score: 
	   test_lasso_cv=cv.lars(A,B, plot.it=FALSE)
	   # Find the best one
	   bestfraction = test_lasso_cv$index[which.min(test_lasso_cv$cv)]
	   #Find Coefficients
	   coef.lasso = predict(lars(A, B), A, s=bestfraction, type=
		 "coefficient", mode="fraction")
	coef.lasso$coefficients
} # end of plain.simple.lasso.coef function


find.i.sis.selected = function(x, y, to.be.selected, already.selected) {
	y.prime = matrix(data=scale(y), ncol=1)
	x.prime = scale(x)
	n = length(y)
	p = ncol(x)
	i.sis.vector = t(x.prime)%*%y.prime
	i.sis.vector = abs(i.sis.vector)
	i.sis.vector = cbind(i.sis.vector, -(1:p))
	i.sis.vector = i.sis.vector[!already.selected,]
	i.sis.vector = i.sis.vector[order(i.sis.vector[,1], i.sis.vector[,2],
		 decreasing = TRUE),]
	i.sis.vector[,2] = -i.sis.vector[,2]
	i.sis.vector[1:(to.be.selected),2]
} # end of find.i.sis.selected function

find.sis.selected = function(x, y, to.be.selected) {
	y.prime = matrix(data=scale(y), ncol=1)
	x.prime = scale(x)
	n = length(y)
	p = ncol(x)
	sis.vector = t(x.prime)%*%y.prime
	sis.vector = abs(sis.vector)
	sis.vector = cbind(sis.vector, -(1:p))
	sis.vector = sis.vector[order(sis.vector[,1], sis.vector[,2], decreasing
		 = TRUE),]
	sis.vector[,2] = -sis.vector[,2]
	sis.vector[1:(to.be.selected),2]
} # end of find.sis.selected function

#sessionInfo()
start.time = date()
sigma1 = matrix(data=0.5, ncol=p, nrow=p)
	sigma1[,4] = rep(sqr.rho, p)
	sigma1[4,] = rep(sqr.rho, p)
	sigma1[,5] = rep(0, p)
	sigma1[5,] = rep(0, p)
	diag(sigma1) = rep(1, p)

my.frauda = function(p=10, n=10, rho=0.5, betastar = c(5, 5, 5,
		-15*sqr.rho, 1, rep(0, p-5)), sigma=sigma1) {
	decomp = svd(sigma)
	sqrt.sigma = decomp$u %*% diag(sqrt(decomp$d)) %*% t(
		decomp$v)
	x = n*p
	x = rnorm(n=x)
	x = matrix(data = x, ncol=p, nrow=n)
	x = x %*% sqrt.sigma
	y = (x%*%betastar)[,1] + rnorm(n=n)
	return(list(x=x, y=y))
} # end my.frauda function


#begin cycle
for (i in 1:replications) {
	#generate frauda
	simulated = my.frauda(p=p, n=n, rho=rho)
	x = simulated$x
	y = simulated$y
check.correlation[i] = cor(y, x[,4])	
	#lasso n-1 features
	#find LASSO, confirm n-1
		lars.obj = lars(x=x, y=y)
		k = nrow(lars.obj$beta)
		numselected = numeric(k)
		for (j in 1:k) {
			numselected[j] = 
			sum(
				abs(lars.obj$beta)[j, ] > 1e-8
			)
		} # end of for loop
		if (any(numselected == n-1)) beta.lasso =
			lars.obj$beta[min(which(numselected==(n-1))),]
		else {
			numselected = abs(numselected - n + 1)
			closest = max(which(numselected ==
				 min(numselected)))
			beta.lasso =lars.obj$beta[closest,]
		} # end else
		rm(numselected)
#report, does lasso select true active features?
		report.lasso[i] = all(abs(beta.lasso[1:5]) > 1e-8)
	#sis n-1 features
		sis.selected = find.sis.selected(x=x, y=y, to.be.selected=(n-1))
		report.sis[i] = all(1:5 %in% sis.selected)

	# i-sis n-1 features
		y.cent = scale(y)
		x.cent = scale(x)
		d = floor(n/log(n))
		# loop to select in batches of size d
		still.needed = n-1
		residual = y.cent
		i.sis.selected = logical(p)
		beta = numeric(p)
			while (still.needed > 0)
			 {
cat("\nstill looping still.needed=", still.needed, "\n")
			new.selected = find.i.sis.selected(x=x.cent, y=residual, 
				to.be.selected = min(d, still.needed),
				 already.selected = i.sis.selected)
			new.selected = sort(new.selected) 
			i.sis.selected[new.selected] = TRUE
			x.prime = x.cent[,i.sis.selected]
			beta[i.sis.selected] = plain.simple.lasso.coef(y = y.cent, x
				 = x.prime)
			#i.sis.selected = (beta != 0)
			still.needed = n - 1 - sum(i.sis.selected)
			residual = as.numeric(y.cent - x.cent[,i.sis.selected] %*%
				matrix(data=beta[i.sis.selected], ncol=1))
		} # end while loop
		report.i.sis[i] = all(1:5 %in% which(i.sis.selected))
sis.obj = SIS(x=x, y=y, family='gaussian', penalty='lasso', tune='cv')
cat("\nalternate =", sis.obj$ix, "\n")
report.alternate[i] = all(1:5 %in% (SIS(x=x, y=y, family='gaussian', penalty='lasso', tune='cv')$ix))

}# end cycle
#report results
cat("\nn = ", n, "p = ", p, "replications=", replications, "\n")
cat("true model selection rate: lasso:", mean(report.lasso), "\n")
cat("true model selection rate: sis:", mean(report.sis), "\n")
cat("true model selection rate: i-sis:", mean(report.i.sis), "\n")
cat("true model selection rate: alternate:", mean(report.alternate), "\n")
cat("\nstart time =", start.time, "end time =", date(), "\n")
########################################