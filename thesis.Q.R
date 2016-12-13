start.time = date()
 set.seed(1984)
#source( "C:/Users/John/Desktop/Code for thesis/thesis.Q.R", echo=T)

# part Q		real data 
#			iterative SIS to a cox model

library(SIS)
library(survival)
library(glmnet)

# data from http://statweb.stanford.edu/~tibs/superpc/staudt.html

file.x = "D:\\staudt.x.txt"
file.cen = "D:\\staudt.status.txt"
file.y = "D:\\staudt.tim.txt"

x = read.table(file=file.x, comment.char="!", header= FALSE)
x = t(as.matrix(x))
x = scale(x)

times = scan(file=file.y)
cens = scan(file=file.cen)
#
times[which(times==0)] = 1e-16

y = Surv(time=times, event= ( cens), type='right') 

sis.obj = SIS(x=x, y=y, family='cox', penalty='lasso')
selected = sis.obj$ix
small.x = x[,selected]
coxph.obj = coxph(y~small.x, data=as.data.frame(small.x))
loglik.model = coxph.obj$loglik[2]

### find loglik for "short" models mm

LOO.models = numeric(length(selected))
for (i in 1:length(selected)) {
	small.x = x[,(selected[-i])]
	LOO.models[i] = coxph(y~small.x, data=
		as.data.frame(small.x))$loglik[2] 
} # end of for loop

### find loglik for extended models
ATE.models  = numeric(20000) #add two extra
potentials = (1:ncol(x))[-selected]
for (i in 1:20000) {
	extras = sample(x = potentials, size = 2)
	small.x = x[,c(selected, extras)]
	ATE.models[i] = coxph(y~small.x, data=
		as.data.frame(small.x))$loglik[2] 
} # end second for loop

quantile(x=(ATE.models - loglik.model), probs=((89:100)/100))

#sink("C:/Users/John/Desktop/Code for thesis/thesis.Q.txt",
#		 append=TRUE, split=TRUE)
#cat("\n<><><><><><><><><><>   <><><><><><><><><><>\n")
#cat("\n\n\n")
#sink()

#OUT = data.frame(found.model=c(loglik.model, rep(NA, 199)), 
#	LOO.models = c(LOO.models, rep(NA, 190)), ATE.models = 
#		ATE.models )

#write.csv(x=OUT, file=
#	'C:/Users/John/Desktop/Code for thesis/thesis.Q.csv')
### 