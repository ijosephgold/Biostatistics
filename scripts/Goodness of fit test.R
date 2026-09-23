#Goodness-of-fit tests

#Chi-square test for equal  sex ratio in Senorita Fish

#create data
female<-273
male<-318

#set up table of observed and expected value:
observed <- matrix(c(female, male), nrow=1, ncol=2) #created a table of observed values
exp <- matrix(c((female+male)*0.5, (female+male)*0.5), nrow=1, ncol=2)

#create a matrix to view data 
fishdata <- matrix(c(observed, exp), nrow=2, ncol=2)
fishdata
#add row and column names
row.names(fishdata) <- (c("Female", "Male"))
colnames(fishdata) <- (c("Observed", "Expected"))
fishdata

#perform X2 test:
X2test<- sum((observed-exp)^2/exp)
X2test #chi-square test statistic, df = #categories-1=1
1-pchisq(X2test, df=1) #P

#perform G2 test
#calculate G:
G2 <- 2*sum(observed*log(observed/exp))
G2 
#get p-value
1-pchisq(G2, df=1) #P 
