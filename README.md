# Thesis
Code for my master's thesis "Independence Screening in High-Dimensional Feature Space" 

This repository hold the code I used in my master's thesis, completed December of 2016.
The code is in the R language.
The items coded are replications of selected published methods for independence screening, applied to simulated and real data.

-----------------------------------------------

Items replicated from Fan and Lv (2008) "Sure Independence Screening in Ultra-High Dimensional Data"

Following section 3.3.2	SIS performance compared to conventional feature 
	selection/screening on simulated data: (independent features, p=
	1000, n=200, true model size=8, 500 replications)

A		DS				
B		LASSO				
C		sis-scad		
D		sis-ds			
E		sis-ds-scad	

Following section 4.2.3	SIS compared to conventional on simulated data, dependent 
	features, p=20000, true model size=14, n=800, 200 replications

F		SIS-SCAD	
G		SIS-DS			
H		SIS-DS-SCAD	

Following section 4.2.3 example III, 	comparison sis, isis, and lasso at selecting true model, 
	simulated data p=1000, n=20,50,70, rho=.5, 
	200 replications

I		n=20			
J		n=50			
K		n=70			

Following section 3.3.3,	real data: leukemia data: classification by sis-scad-nb 

L		SIS-SCAD-naive Bayes	done

Items replicated from Fan and Fan (2008) "High Dimensional Classification Using Features Annealed Independence Rules"


M		Following section 5.2.1, Nearest shrunken centroid classification of Leukemia data
			misclassification rate and mean training error,     
			testing error, and number of selected genes
N		Following section 5.2.3, Features Annealed Independence Rules classification of 
			leukemia data, misclassification rate and mean	      done
			training error, testing error, and number of selected
			genes
		
Items replicated from Fan, Feng, and Wu (2010) "High-Dimesnional Variable Selection for Cox's Proportional Hazard Model"

Following Section 5, cases 5 and 6:
O		Using simulated data, fit sets of data with "vanilla" 
			iterative SIS to a cox model
P		fit same simulated data using Variant 1 iterative SIS, 
			compare in terms of mean model size, L2 norm of 
			difference of true and estimated betas, and frequency 
			of including true model.

Following the methods of section 6, but using Lymphoma data from http://statweb.stanford.edu/~tibs/superpc/

Q		fit real Lymphoma data data using "vanilla"
			iterative SIS, compare to published result.

------------------------------

References:

Fan, J., & Fan, Y. (2008). High dimensional classification using features annealed independence rules.
Annals of statistics, 36(6), 2605.

Fan, J., Feng, Y., & Wu, Y. (2010). High-dimensional variable selection for Cox’s proportional hazards model.
In Borrowing Strength: Theory Powering Applications–A Festschrift for Lawrence D. Brown (pp. 70-86). Institute
of Mathematical Statistics.

Fan, J., & Lv, J. (2008). Sure independence screening for ultrahigh dimensional feature space.
Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(5), 849-911.

----------------------------
