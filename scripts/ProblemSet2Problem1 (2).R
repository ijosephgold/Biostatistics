#Goodness-of-fit tests

#Chi-square test for equal  sex ratio in Hymenoptera Bees

#create data
female<-371
male<-250

#set up table of observed and expected value:
observed <- matrix(c(female, male), nrow=1, ncol=2) #created a table of observed values
exp <- matrix(c((female+male)*0.5, (female+male)*0.5), nrow=1, ncol=2)

#create a matrix to view data 
beedata <- matrix(c(observed, exp), nrow=2, ncol=2)
beedata
#add row and column names
row.names(beedata) <- (c("Female", "Male"))
colnames(beedata) <- (c("Observed", "Expected"))
beedata

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
