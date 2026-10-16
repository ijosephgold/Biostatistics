data <- read.csv("Data/Dihybrid Cross Gtest - Sheet1.csv")
view(data)
#Goodness-of-fit tests

#Chi-square test for equal  sex ratio in Hymenoptera Bees

#create data
yellowround<-96
yellowwrinkled<-18
greenround<-24
greenwrinkled<-7

n<-sum(yellowround, yellowwrinkled, greenround, greenwrinkled)  #created a table of observed values
exp<-matrix(c((n*9/16), (n*3/16), (n*3/16), (n*1/16)), nrow=1, ncol=4)
observed <- matrix(c(yellowround, yellowwrinkled, greenround, greenwrinkled), nrow=1, ncol=4) #created a table of observed values

#create a matrix to view data 
peadata <- matrix(c(observed, exp), nrow=4, ncol=2)
peadata
#add row and column names
row.names(peadata) <- (c("yellowround", "yellowwrinkled", "greenround", "greenwrinkled"))
colnames(peadata) <- (c("Observed", "Expected"))
peadata            
            
#set up table of observed and expected value:

#perform G2 test
#calculate G:
G2 <- 2*sum(observed*log(observed/exp))
G2 
#get p-value
1-pchisq(G2, df=3) #P 
