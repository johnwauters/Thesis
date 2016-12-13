
## "C:/Users/John/Desktop/Code for thesis/my.methods.R"

## assumed format of x, y and beta
## x is an n by p matrix
## y is a column matrix n by 1
## beta is a column matrix p by 1
## penalized least squares:
##  pls(beta) = (1/(2*n))sum(as.numeric(y - x%*%beta)**2
##		) + penalty(lambda, beta)
## where penalty(lambda, beta) takes form of
##	sum from j=1 to d of p-sub-lambda-sub-j(|beta-sub-j|)
##
##
##
##

mylq = function(x, q=2) ## my function for computing an Lq norm
{
	if (class(x) != 'numeric' && class(x) != "matrix" && class(x) != 
			'integer')	{cat("incorrect 	type for Lq norm"); 
				return(NULL)}
	if (class(x) == 'matrix' && min(dim(x)) != 1 ) {echo("incorrect 			input type for Lq norm"); return(NULL)}
	if (!is.numeric(q) || (q < 0) ) {cat("improper q"); return(NULL)}
	x = as.numeric(x)
	x = abs(x)
	if (q == 0) return(sum(x > 0))
	if (q == Inf) return(max(x))
	x = x**q
	x = sum(x)
	x = x**(1/q)
	return(x)
} # ## end of function

myls = function(x, y, minstep = ((.Machine$double.eps))) {
 	## my function to find un-penalized least-squares solutions, 
	## assumes y is either a column-matrix or a numeric vector, x is 
	## a matrix with nrows(x) = length(y),  also assumed is that y 
	## and x have both been subject to scale() function
	n = length(y)
	p = dim(x)[2]
	step=1
	loop.control = 0
	beta=numeric(p)
	nextbeta = beta
	x = scale(x)
	y = as.numeric(scale(y))
	while (step >= minstep && log(loop.control) < n*29.0){ 
		for( i in 1:n) {
			betaback = beta
			betaback[i] = betaback[i] - step
			back = sum((y - x%*%betaback)**2) + 0:0
			betafore = beta
			betafore[i] = betafore[i] + step
			fore = sum((y - x%*%betafore)**2) + 0:0
			center = sum((y - x%*%beta)**2) + 0:0
  			if (back == min(back, fore, center)) { nextbeta[i] = 
				beta[i] - step}
  			if (fore == min(back, fore, center)) { nextbeta[i] = 
				beta[i] + step}

  			if (center == min(back, fore, center)) { 
				nextbeta[i]= beta[i]
				 }
		} ## end of foor loop
		if (center == sum((y - x%*%nextbeta)**2)
			) {step = step/2}
		else {step = step*0.95
		 }
		beta = nextbeta
		loop.control= loop.control + 1
cat("beta", beta, "s", step, loop.control, "--\n")
	} ## end of while loop
	if (log(loop.control) > 29*n) {cat("loop overflow error"); 
		return(NULL) }
	for (j in 1:n) {
		abs.beta = abs(beta)
		if (abs.beta[j] <= 2*minstep) {beta[j] = 0}
		} # end of for loop
return(matrix(data=beta, ncol=1))
} ### end of function

# ### END OF FILE ############################