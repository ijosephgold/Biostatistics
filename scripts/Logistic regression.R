#Logistic regression

#load libraries
library(lme4)
library(car)
library(pscl) #for producing rho-squared
library(ggplot2)

#load data
perchdata <- read.csv("Data/Surfperch data.csv")
perchdata

#create generalized linear model
glm <- glm(mated ~ male.size, family = "binomial"(link="logit"), data=perchdata)
#look at the coefficients (estimate)
summary(glm)
#convert body mass to odds
exp(0.6563) #e^b1
#run log likelihood ratio test
Anova(glm)

#calculate rho-squared
pR2(glm)
#llh is the log likelihood from the fitted model
#llhnull is the log likelihood from the intercept-only model
#try calculating rho-square by hand
rho2 = 1 - (21.7156206/32.0517739)

#G2 method
#This does what Anova() does behind the scenes
G2 <- 2*((-21.7156206)-(-32.0517739)) #same as LR Chisq
1-pchisq(G2, df=1) #gives the P-value

#make simple plot
plot(perchdata$male.size, perchdata$mated, xlab="Male Size", ylab="Probability of Mating")
curve(predict(glm,data.frame(male.size=x),type="resp"),add=TRUE) # draws a curve

#make a nicer plot with ggplot()
predicted.data<-as.data.frame(predict(glm, type="response", se=TRUE))

ggplot(perchdata, aes(x=perchdata$male.size, y=perchdata$mated)) +
  theme_bw() +
  theme(axis.text.x=element_text(color="black", size=12), 
        axis.text.y=element_text(color="black", size=12), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
        geom_point() + geom_smooth(data = predicted.data, aes(y=fit)) +
        labs(x="Male Size (cm)", y="Pr(Mating Success)")
