start.time = date()
#source("C:/Users/John/Desktop/Code for thesis/thesis.MN.R", echo=T)

library(SIS)
library(pamr)
library(ncvreg)
library(MASS)

abT = function(x.vec, y){
	x = as.numeric(x.vec)
	y = as.numeric(y)
	x.1 = x[(y==1)]
	x.0 = x[(y==0)]
	abs(mean(x.1)-mean(x.0))/sqrt(var(x.1)/length(x.1) + 
		var(x.0)/length(x.0))
} # end abT function

optnum.fair = function(x, y) { # y is assumed numeric, not factor
					   # y is assumed to have values 0 and 1
					# finds the m.hat.0 from fan & fan pg 8
	x = scale(x)
	n = length(y)
	if (nrow(x) != n) stop('\ninput size mismatch optnum.fair\n')
	p = ncol(x)
	acc = numeric(p) # acc_umulator for values
	alpha.hat.2 = numeric(p)
	n1 = sum(y == 1)
	n2 = n - n1
	if (n1 < n2) {temp = n1; n1 = n2; n2 = temp; rm(temp)}
	alpha.hat.2 = (colMeans(x[which(y == 1),]) - colMeans(x[which(y !=
			 1),]))^2
	alpha.hat.2 = alpha.hat.2[order(alpha.hat.2, decreasing=TRUE)]
	for (m in 1:p) {
		acc[m] = (sum(alpha.hat.2[1:m]) + m*(n2-n1)/(n1*n2))^2
		acc[m] = acc[m] / (n*m/(n1*n2) + sum(alpha.hat.2[1:m]))
	} # end of for loop
	R = cor(x[,order(alpha.hat.2, decreasing=TRUE)])
	m = 2
	while(which.max(acc[1:(m-1)])  > (.5*(m-200)) )  {
			acc[m] = acc[m]/ max(eigen(R[1:m,1:m])$values)
			m = m + 1
	} # end while loop	
	which.max(acc[1:(m-2)])
} # end of optnum.fair

delta.fair = function(train.x, train.y, newx, selected) {
	n = length(train.y)
	m.train = length(selected)
	p = ncol(train.x)
	if (nrow(train.x) != n) {stop('error in delta.fair input training data')}
	if (length(newx) == p) {newx = matrix(data=newx, nrow=1)}
	if (ncol(newx) != ncol(train.x)) stop('inproper newx in delta.fair')
	mu.hat = colMeans(x.train[,selected])
	sig.hat.2 = apply(x.train[,selected], 2, var)
	a.hat = colMeans(x.train[(train.y == 1),selected]) - 
		colMeans(x.train[(train.y==0), selected])
	n.new = nrow(newx)
	pred = numeric(n.new)
	train.x = train.x[,selected]
	newx = newx[,selected]
	for (i in 1:n.new) {
		pred[i] = sum(a.hat*(newx[i,] - mu.hat)/sig.hat.2)
	} # end for loop	
	as.numeric(pred > 0)
} # end delta.fair function

# data from prostate dataset in SIS package
# from paper Sing et al (2002)
# prostate.test 34 observations on 12600 variables, 25 tumor 9 normal
# column 12601 is response, 0 = tumor, 1 = normal
# prostate.train 102 observations on 12600 variables, 52 tumor
# 50 normal

data(prostate.test)
data(prostate.train)
	###
x.test = scale(prostate.test[,1:12600])
y.test = as.numeric(prostate.test[,12601])
x.train = scale(prostate.train[,1:12600])
y.train = as.numeric(prostate.train[,12601])

## section for train error, test error, and num genes retained for
## nearest shrunken centroids

mydata.train = list(x =t(x.train), y = factor(y.train)) 
pamr.obj = pamr.cv(fit = pamr.train(data=mydata.train),
		data=mydata.train)
optim.loc = which.min(pamr.obj$error)
my.threshold = pamr.obj$threshold[optim.loc]
nsc.num.genes = pamr.obj$size[optim.loc]
my.pamr.train = pamr.train( data = mydata.train, threshold =
		my.threshold)
nsc.train.error = sum(pamr.predict(fit=my.pamr.train, newx=t(x.train),
		threshold=my.threshold) != factor(y.train))
nsc.test.error = sum(pamr.predict(fit=my.pamr.train, newx=t(x.test),
		threshold=my.threshold) != factor(y.test))
rm(mydata.train, pamr.obj, optim.loc, my.threshold, my.pamr.train)

## end nearest shrunken centroids part

## FAIR

m.train = optnum.fair(x=x.train, y=y.train)
p = ncol(x.train)
Tvec = numeric(p)
for (i in 1:p) {
	Tvec[i] = abT(x.vec=x.train[,i], y=y.train)
} # end for loop
selected = (order(Tvec, decreasing=TRUE)[1:m.train])


fair.train.error = sum(delta.fair(x.train, y.train, x.train, selected) != y.train)
fair.test.error =  sum(delta.fair(x.train, y.train, x.test, selected) != y.test)
fair.genes.selected = length(selected)

sink("C:/Users/John/Desktop/Code for thesis/thesis.MN.txt",
		 append=TRUE, split=TRUE)
cat("\n<><><><><><><><><><>   <><><><><><><><><><>\n")
cat("\nnsc train error=", nsc.train.error, "test error=", nsc.test.error, 
	"nsc genes selected=", nsc.num.genes, "\n")
cat("\nfair train error=", fair.train.error, "test error=", fair.test.error, 
	"fair genes selected=", fair.genes.selected, "\n")
cat("\n\n\n")
sink()
# ####
x.comb = rbind(x.test, x.train)
y.comb = c(y.test, y.train)
x.1 = x.comb[(y.comb == 1),]
x.0 = x.comb[(y.comb == 0),]
n1 = nrow(x.1)
n0 = nrow(x.0)

# gamma takes values .4, .5, and .6
cycle.length = 100  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 100
# MC cycle 1
MC1nsctesterror = numeric(cycle.length)
MC1fairtesterror = numeric(cycle.length)
MC1difftesterror = numeric(cycle.length)
MC1nscnumgenes = numeric(cycle.length)
MC1fairnumgenes = numeric(cycle.length)
MC1diffnumgenes = numeric(cycle.length)
pull1 = as.integer(.4*n1)
pull0 = as.integer(.4*n0)
for( i in 1:cycle.length) {
	pick1 = sample(1:n1, size=pull1)
	pick0 = sample(1:n0, size=pull0)
	x.train = rbind(x.1[pick1,], x.0[pick0,])
	y.train = c(rep(1, pull1), rep(0, pull0))
	x.test = rbind(x.1[-pick1,], x.0[-pick0,])
	y.test = c(rep(1, n1-pull1), rep(0, n0-pull0))
##nsc
	mydata.train = list(x =t(x.train), y = factor(y.train)) 
	pamr.obj = pamr.cv(fit = pamr.train(data=mydata.train),
			data=mydata.train)
	optim.loc = which.min(pamr.obj$error)
	my.threshold = pamr.obj$threshold[optim.loc]
	MC1nscnumgenes[i] = pamr.obj$size[optim.loc]
	my.pamr.train = pamr.train( data = mydata.train, threshold =
			my.threshold)
	MC1nsctesterror[i] = sum(pamr.predict(fit=my.pamr.train,
		 newx=t(x.test), threshold=my.threshold) != factor(y.test))
##fair
	m.train = optnum.fair(x=x.train, y=y.train)
	p = ncol(x.train)
	Tvec = numeric(p)
	for (J in 1:p) {
		Tvec[J] = abT(x.vec=x.train[,J], y=y.train)
	} # end for loop
	selected = (order(Tvec, decreasing=TRUE)[1:m.train])

	MC1fairtesterror[i] =  sum(delta.fair(x.train, y.train, x.test, 
		selected)  != y.test)
	MC1fairnumgenes[i] = length(selected)

	MC1difftesterror[i] = MC1nsctesterror[i] - MC1fairtesterror[i]
	MC1diffnumgenes[i] = MC1nscnumgenes[i] - MC1fairnumgenes[i]
} # end first MC cycle

######
# MC cycle 2
MC2nsctesterror = numeric(cycle.length)
MC2fairtesterror = numeric(cycle.length)
MC2difftesterror = numeric(cycle.length)
MC2nscnumgenes = numeric(cycle.length)
MC2fairnumgenes = numeric(cycle.length)
MC2diffnumgenes = numeric(cycle.length)
pull1 = as.integer(.5*n1)
pull0 = as.integer(.5*n0)
for( i in 1:cycle.length) {
	pick1 = sample(1:n1, size=pull1)
	pick0 = sample(1:n0, size=pull0)
	x.train = rbind(x.1[pick1,], x.0[pick0,])
	y.train = c(rep(1, pull1), rep(0, pull0))
	x.test = rbind(x.1[-pick1,], x.0[-pick0,])
	y.test = c(rep(1, n1-pull1), rep(0, n0-pull0))
##nsc
	mydata.train = list(x =t(x.train), y = factor(y.train)) 
	pamr.obj = pamr.cv(fit = pamr.train(data=mydata.train),
			data=mydata.train)
	optim.loc = which.min(pamr.obj$error)
	my.threshold = pamr.obj$threshold[optim.loc]
	MC2nscnumgenes[i] = pamr.obj$size[optim.loc]
	my.pamr.train = pamr.train( data = mydata.train, threshold =
			my.threshold)
	MC2nsctesterror[i] = sum(pamr.predict(fit=my.pamr.train,
		 newx=t(x.test), threshold=my.threshold) != factor(y.test))
##fair
	m.train = optnum.fair(x=x.train, y=y.train)
	p = ncol(x.train)
	Tvec = numeric(p)
	for (J in 1:p) {
		Tvec[J] = abT(x.vec=x.train[,J], y=y.train)
	} # end for loop
	selected = (order(Tvec, decreasing=TRUE)[1:m.train])

	MC2fairtesterror[i] =  sum(delta.fair(x.train, y.train, x.test, 
		selected)  != y.test)
	MC2fairnumgenes[i] = length(selected)

	MC2difftesterror[i] = MC2nsctesterror[i] - MC2fairtesterror[i]
	MC2diffnumgenes[i] = MC2nscnumgenes[i] - MC2fairnumgenes[i]
} # end second MC cycle

###### end of MC2
# MC cycle 3
MC3nsctesterror = numeric(cycle.length)
MC3fairtesterror = numeric(cycle.length)
MC3difftesterror = numeric(cycle.length)
MC3nscnumgenes = numeric(cycle.length)
MC3fairnumgenes = numeric(cycle.length)
MC3diffnumgenes = numeric(cycle.length)
pull1 = as.integer(.6*n1)
pull0 = as.integer(.6*n0)
for( i in 1:cycle.length) {
	pick1 = sample(1:n1, size=pull1)
	pick0 = sample(1:n0, size=pull0)
	x.train = rbind(x.1[pick1,], x.0[pick0,])
	y.train = c(rep(1, pull1), rep(0, pull0))
	x.test = rbind(x.1[-pick1,], x.0[-pick0,])
	y.test = c(rep(1, n1-pull1), rep(0, n0-pull0))
##nsc
	mydata.train = list(x =t(x.train), y = factor(y.train)) 
	pamr.obj = pamr.cv(fit = pamr.train(data=mydata.train),
			data=mydata.train)
	optim.loc = which.min(pamr.obj$error)
	my.threshold = pamr.obj$threshold[optim.loc]
	MC3nscnumgenes[i] = pamr.obj$size[optim.loc]
	my.pamr.train = pamr.train( data = mydata.train, threshold =
			my.threshold)
	MC3nsctesterror[i] = sum(pamr.predict(fit=my.pamr.train,
		 newx=t(x.test), threshold=my.threshold) != factor(y.test))
##fair
	m.train = optnum.fair(x=x.train, y=y.train)
	p = ncol(x.train)
	Tvec = numeric(p)
	for (J in 1:p) {
		Tvec[J] = abT(x.vec=x.train[,J], y=y.train)
	} # end for loop
	selected = (order(Tvec, decreasing=TRUE)[1:m.train])

	MC3fairtesterror[i] =  sum(delta.fair(x.train, y.train, x.test, 
		selected)  != y.test)
	MC3fairnumgenes[i] = length(selected)

	MC3difftesterror[i] = MC3nsctesterror[i] - MC3fairtesterror[i]
	MC3diffnumgenes[i] = MC3nscnumgenes[i] - MC3fairnumgenes[i]
} # end third MC cycle

#### end of MC3

## reporting section
sink("C:/Users/John/Desktop/Code for thesis/thesis.MN.txt",
		append=TRUE, split=TRUE)
summary(MC1fairtesterror)
summary(MC1fairnumgenes)
summary(MC1nsctesterror)
summary(MC1nscnumgenes)
summary(MC1difftesterror)
summary(MC1diffnumgenes)

cat("\nstart time =", start.time, "end time=", date(), "\n")

cat("\n#########################################\n")
sink()

OUT = data.frame(MC1fairtesterror = MC1fairtesterror, MC1fairnumgenes = MC1fairnumgenes, MC1nsctesterror = MC1nsctesterror, MC1nscnumgenes = MC1nscnumgenes, MC1difftesterror = MC1difftesterror, MC1diffnumgenes = MC1diffnumgenes,

MC2fairtesterror = MC2fairtesterror, MC2fairnumgenes = MC2fairnumgenes, MC2nsctesterror = MC2nsctesterror, MC2nscnumgenes = MC2nscnumgenes, MC2difftesterror = MC2difftesterror, MC2diffnumgenes = MC2diffnumgenes,

MC3fairtesterror = MC3fairtesterror, MC3fairnumgenes = MC3fairnumgenes, MC3nsctesterror = MC3nsctesterror, MC3nscnumgenes = MC3nscnumgenes, MC3difftesterror = MC3difftesterror, MC3diffnumgenes = MC3diffnumgenes)

write.csv(x=OUT, file=
	'C:/Users/John/Desktop/Code for thesis/thesis.MN.csv')
# ##########################################