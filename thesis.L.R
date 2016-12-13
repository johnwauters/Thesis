start.time = date()
#source("C:/Users/John/Desktop/Code for thesis/thesis.L.R", echo=T)
# Leukemia datasets URLs -- TAB delimited

# may need to get other excell format

#http://www.broadinstitute.org/mpr/publications/projects/Leukemia
#/data_set_ALL_AML_train.txt

#http://www.broadinstitute.org/mpr/publications/projects/Leukemia
#/data_set_ALL_AML_independent.txt

# test.df = 
# read.table(file=file.choose(), sep="\t", header=T, row.names=NULL)


# SIS-SCAD-NB
library(SIS)
# NSC
library(pamr)
# FAIR
#library(fairselect) # does it exist?

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
#browser()
	R = cor(x[,order(alpha.hat.2, decreasing=TRUE)])
	m = 2
cat(acc[1:6], "\n")
	while(which.max(acc[1:(m-1)])  > (.5*(m-200)) )  {
			acc[m] = acc[m]/ max(eigen(R[1:m,1:m])$values)
			m = m + 1
	} # end while loop	
	which.max(acc[1:(m-2)])
} # end of optnum.fair


n = 38
d = floor(2*n/log(n))
 #need leukemia dataset
data(leukemia.test) # first 20 ALL, next 14 AML
data(leukemia.train) # first 27 ALL, next 11 AML
x.test = scale(leukemia.test)
y.test = c( rep(1,20), rep(0, 14) ) #class 1 = ALL, class 0 = AML
x.train = scale(leukemia.train)
y.train = c(rep(1, 27), rep(0, 11))

## section for train error, test error, and num genes retained for
## nearest shrunken centroids

mydata.train = list(x =t(x.train), y = factor(y.train)) 
pamr.obj = pamr.cv(fit = pamr.train(data=mydata.train),
		data=mydata.train)
optim.loc = which.min(pamr.obj$error)
my.threshold = pamr.obj$threshold[optim.loc]
 num.genes = pamr.obj$size[optim.loc]
my.pamr.train = pamr.train( data = mydata.train, threshold =
		my.threshold)
train.error = sum(pamr.predict(fit=my.pamr.train, newx=t(x.train),
		threshold=my.threshold) != factor(y.train))
test.error = sum(pamr.predict(fit=my.pamr.train, newx=t(x.test),
		threshold=my.threshold) != factor(y.test))
rm(mydata.train, pamr.obj, optim.loc, my.threshold, my.pamr.train)

## end nearest shrunken centroids part

## SIS-SCAD PART

#sis.scad.obj = tune.fit(x=x.train, y=y.train, family='gaussian', 
#	penalty='SCAD', tune='bic')

sis.obj = abs(t(x.train)%*%scale(y.train))
sis.obj = as.numeric(sis.obj)


sis.features.chosen = order(sis.obj, decreasing=TRUE)[2:13]


#sis.num.chosen=optnum.fair(x.train, y.train)
 #now do SCAD to prune down to 16
library(ncvreg)
#ncvreg(X=x.train[,sis.features.chosen], y=y.train, penalty='SCAD')
library(MASS)
fit.train = lda(x=x.train[,(sis.features.chosen)], grouping=y.train)
ld.train.error = sum((predict(object=fit.train, newdata= 
		x.train[,sis.features.chosen])$class) != factor(y.train))
ld.test.error = sum((predict(object=fit.train, newdata= 
		x.test[,sis.features.chosen])$class) != factor(y.test))

cat("\n<><><>train", train.error, "test", test.error, "-\n")

## NB section
library(klaR)
nb.obj = NaiveBayes(x = x.train[,sis.features.chosen], grouping= 
			factor(y.train))
nb.train.error = sum(predict(object=nb.obj, newdata= 
		x.train[,sis.features.chosen])$class != factor(y.train))
nb.test.error = sum(predict(object=nb.obj, newdata= 
		x.test[,sis.features.chosen])$class != factor(y.test))

sink("C:/Users/John/Desktop/Code for thesis/thesis.L.txt", append=TRUE, split=TRUE)
cat("\n<><><><><><><><><><>     <><><><><><><><><><>\n")
cat("start time =", start.time, "\n")
cat("\nnsc train error=", train.error, "test error=", test.error, "\n")
cat("\nsis-scad-ld train error=", ld.train.error, "test error=", ld.test.error, 
		"\n")
cat("\nsis-scad-nb train error=", nb.train.error, "test error=",
		 nb.test.error, "\n")
cat("\n\n\n")
sink()
# ####

x.class.1 = rbind(x.train[(y.train == 1),], x.test[(y.test == 1),])
x.class.0 = rbind(x.train[(y.train == 0),], x.test[(y.test == 0),])

set.seed(1984)
## begin improvised loop
Q = 300
improve.genes = numeric(Q) 
improve.train = numeric(Q) 
improve.test = numeric(Q) 
for (i in 1:Q) {

sample1 = sample(x=1:47, size=23)
sample0 = sample(x=1:25, size=12)

x.temp.train = scale(rbind(x.class.1[sample1,], x.class.0[sample0,]))
x.temp.test = scale(rbind(x.class.1[-sample1,], x.class.0[-sample0,]))
y.temp.train = c(rep(1,23), rep(0, 12))
y.temp.test = c(rep(1, 24), rep(0, 13))


mydata.train = list(x =t(x.temp.train), y = factor(y.temp.train)) 
pamr.obj = pamr.cv(fit = pamr.train(data=mydata.train),
		data=mydata.train)
optim.loc = which.min(pamr.obj$error)
my.threshold = pamr.obj$threshold[optim.loc]
nsc.num.genes = pamr.obj$size[optim.loc]
my.pamr.train = pamr.train( data = mydata.train, threshold =
		my.threshold)
nsc.train.error = sum(pamr.predict(fit=my.pamr.train, newx= 
		t(x.temp.train), threshold=my.threshold) != 
		factor(y.temp.train))
nsc.test.error = sum(pamr.predict(fit=my.pamr.train, newx=t(x.temp.test),
		threshold=my.threshold) != factor(y.temp.test))
# ### SIS version

sis.obj = abs(t(x.temp.train)%*%scale(y.temp.train))
sis.obj = as.numeric(sis.obj)

optnum = optnum.fair(x=x.temp.train, y=y.temp.train)

sis.features.chosen = order(sis.obj, decreasing=TRUE)[1:optnum]
x.temp.train = x.temp.train[,sis.features.chosen]
x.temp.test = x.temp.test[,sis.features.chosen]

mydata.train = list(x =t(x.temp.train), y = factor(y.temp.train)) 
pamr.obj = pamr.cv(fit = pamr.train(data=mydata.train),
		data=mydata.train)
optim.loc = which.min(pamr.obj$error)
my.threshold = pamr.obj$threshold[optim.loc]
sisnsc.num.genes = pamr.obj$size[optim.loc]
my.pamr.train = pamr.train( data = mydata.train, threshold =
		my.threshold)
sisnsc.train.error = sum(pamr.predict(fit=my.pamr.train, newx= 
		t(x.temp.train), threshold=my.threshold) != 
		factor(y.temp.train))
sisnsc.test.error = sum(pamr.predict(fit=my.pamr.train, newx=t(x.temp.test),
		threshold=my.threshold) != factor(y.temp.test))
cat("\nnum genes", (nsc.num.genes - sisnsc.num.genes), "train err", 
	(nsc.train.error-sisnsc.train.error), "test err", (nsc.test.error - 
	sisnsc.test.error), "\n") 

	improve.genes[i] = (nsc.num.genes - sisnsc.num.genes)
	improve.train[i] = (nsc.train.error - sisnsc.train.error)
	improve.test[i] = (nsc.test.error - sisnsc.test.error)

} # end improvised loop

sink("C:/Users/John/Desktop/Code for thesis/thesis.L.txt", append=TRUE, split=TRUE)

cat(Q, "\n")
summary(improve.genes)
sd(improve.genes)
quantile(improve.genes, probs=c(.025, .975))
summary(improve.train)
sd(improve.train)
quantile(improve.train, probs=c(.025, .975))
summary(improve.test)
sd(improve.test)
quantile(improve.test, probs=c(.025, .975))



cat("\nstart time =", start.time, "end time=", date(), "\n")

cat("\n#########################################\n")
sink()

# ##########################################