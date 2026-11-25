#ANCOVA

#load libraies
library(tidyr)
library(dplyr)
library(car)
library(moments)
library(ggplot2)

#let's look at mussel byssal thread force vs. size at two sites
#load data
threads <- read.csv("Data/Byssal threads full.csv")

#create GLM
#y = strength, x1 = length, x2 = site
model1 <- lm(Strength ~ Length*Site, data=threads)

#check residuals
model1res <- resid(model1)
qqPlot(model1res) #not terrible
shapiro.test(model1res) #not good

#test for homogeneity of variance
# call the fitted y values
fitted(model1) 
#plot the fitted values against the residuals
plot(model1res~fitted(model1), ylab = "Residuals", xlab = "Predicted log y") 
#it looks ok

#try log transformation 
threads$logstrength <- log(threads$Strength)
threads$loglength <- log(threads$Length)

#create GLM
model2 <- lm(logstrength ~ loglength*Site, data=threads)

#check residuals
model2res <- resid(model2)
qqPlot(model2res) #not terrible
shapiro.test(model2res) #good

#test for homogeneity of variance
# call the fitted y values
fitted(model2) 
#plot the fitted values against the residuals
plot(model2res~fitted(model2), ylab = "Residuals", xlab = "Predicted log y") 
#it looks a little wedge shaped

#let's analyze model2
anova(model2) #the interaction isn't significant -- slopes are the same
#let's create a model without the interaction
#test is the heights of the lines are the same
model3 <- lm(logstrength~loglength+Site, data=threads)
anova(model3) #line heights are the same
#site doesn't have an effect but larger mussels are stronger

#graph data
predWeight <- predict(model1) #Gets the predicted values from the regression lines in the ANCOVA
graphdata <- cbind(threads, predWeight) #attaches those predictions to the dataset
#make plot
ggplot(data=graphdata, aes(Length, Strength, color=Site)) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=12), 
        axis.text.x=element_text(color="black", size=10), 
        axis.text.y=element_text(color="black", size=10), 
        axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=12),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point() + geom_line(aes(y=predWeight)) +
  labs(x="Length (mm)", y="Strength (N)", fill="Site")

#let's look at chirp pulse rate vs. temperature comparing two different cricket species
#load data
chirps <- read.csv("Data/Cricket data.csv")
head(chirps)

#create GLM
#y = pulse rate, x1 = temperature, x2 = species
model4 <- lm(PulseRate ~ Temperature*Species, data=chirps)
#let's ignore testing assumptions for now since the code for that is above
#run the ANCOVA
anova(model4)
#the interaction is significant
#the slopes of the lines are not the same
#the effect of temperature on chirp rate depends upon species

#run separate linear regressions for each species

#filter the data by species
GrRbdata <- chirps %>%
  filter(Species == "Gryllus rubens")
GrTxdata <- chirps %>%
  filter(Species == "Gryllus texensis")

#create GLMs
GrRbmodl <- lm(PulseRate~Temperature, data=GrRbdata)
anova(GrRbmodl)
summary(GrRbmodl)
GrTxmodl <- lm(PulseRate~Temperature, data=GrTxdata)
anova(GrTxmodl)
summary(GrTxmodl)
#temperature affects chirp rate in both species
#but chirp rate is higher in G. texensis

#graph data
predWeight <- predict(model4) #Gets the predicted values from the regression lines in the ANCOVA
graphdata <- cbind(chirps, predWeight) #attaches those predictions to the dataset
#make plot
ggplot(data=graphdata, aes(Temperature, PulseRate, color=Species)) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=12), 
        axis.text.x=element_text(color="black", size=10), 
        axis.text.y=element_text(color="black", size=10), 
        axis.title.x = element_text(color="black", size=12), 
        axis.title.y = element_text(color="black", size=12),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point() + geom_line(aes(y=predWeight)) +
  labs(x="Temperature (C)", y="Chirp rate", fill="Species")
