#H_0: The observed data fits the expected 9:3:3:1 ratio
#H_A: The observed data does not fit the expected 9:3:3:1 ratio
#Critical value: 7.815
data <- read.csv("Data/Dihybrid Cross Gtest - Sheet1.csv")
view(data)
#Goodness-of-fit tests

#create Observed data
yellowround<-96
yellowwrinkled<-18
greenround<-24
greenwrinkled<-7

n<-sum(yellowround, yellowwrinkled, greenround, greenwrinkled)  #created a table of observed values
exp<-matrix(c((n*9/16), (n*3/16), (n*3/16), (n*1/16)), nrow=1, ncol=4)
observed <- matrix(c(yellowround, yellowwrinkled, greenround, greenwrinkled), nrow=1, ncol=4) 
#created a table of observed values

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
G2 #Observe G value is less than critical 7.815
#get p-value
p_value <- pchisq(G2, df = 3, lower.tail = FALSE)
p_value

#Conclusion: Since the p-value is greater than 0.05 (0.077), we fail to reject the null hypothesis 
#and conclude that there's not enough evidence say they differ from expected 9:3:3:1 ratio. 