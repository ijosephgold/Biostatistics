#Factorial ANOVA

#load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)

#simple two-way ANOVA (Model I)
#we'll use the clipped vs. weeded data

#load data
grassdata <- read.csv("Data/Bromus data.csv")

#make x variables factors
grassdata$ClipTrt <- factor(grassdata$ClipTrt)
grassdata$WeedTrt <- factor(grassdata$WeedTrt)

#create the two-way GLM
model1 <- lm(NoFlrs~ClipTrt*WeedTrt, data = grassdata)
anova(model1)

#check residuals
model1res <- resid(model1)
qqPlot(model1res)
shapiro.test(model1res) #not great

#try sqrt transformation
grassdata$sqrtnoflrs <- sqrt(grassdata$NoFlrs)
model2 <- lm(sqrtnoflrs~ClipTrt*WeedTrt, data = grassdata)
anova(model2)

#check residuals
model2res <- resid(model2)
qqPlot(model2res)
shapiro.test(model2res) #much better

#since the interaction term is significant, we'll perform separate t-tests
#clipped treatment
bartlett.test(sqrtnoflrs~ClipTrt, data=grassdata) #variances are equal
t.test(sqrtnoflrs~ClipTrt, var.equal=TRUE, data=grassdata)
#weed treatment
bartlett.test(sqrtnoflrs~WeedTrt, data=grassdata) #variances are unequal
t.test(sqrtnoflrs~WeedTrt, var.equal=FALSE, data=grassdata)

#make barplot
graphdata <- as.data.frame(emmeans(model2, ~ ClipTrt*WeedTrt))
graphdata 
ggplot(graphdata, aes(x=ClipTrt, y=emmean, fill=factor(WeedTrt), group=factor(WeedTrt))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_bar(stat="identity", position="dodge", size=0.6) + #determines the bar width
  geom_errorbar(aes(ymax=emmean+SE, ymin=emmean-SE), stat="identity", position=position_dodge(width=0.9), width=0.1) + #adds error bars
  labs(x="", y="Sqrt(No Flowers)", fill="") #labels the x and y axes

#mixed model ANOVA (Model III)
#let's look at the white clover genotype study

#load data
clovers <- read.csv("Data/Clover data.csv")

#make x variables factors
clovers$Genotype <- factor(clovers$Genotype)
clovers$Location <- factor(clovers$Location)

#let's not worry about assumptions right now

#create general linear mixed effects models
#random intercept model
model3 <- lmer(Frequency~Genotype + (1|Location), data=clovers) #|1 sets location as random
#interaction model
model4 <- lmer(Frequency~Genotype + (1|Location) + (1|Location:Genotype), data=clovers) #1|Location:Genotype codes the interaction term
#log likelihood ratio test
anova(model3, model4) #genotype frequency depends upon location
anova(model3) #genotype frequencies are significantly different

#make barplot
#first create a data frame that ggplot can call
graphdata <- clovers %>%
  group_by(Location, Genotype) %>%
  summarise(
    mean = mean(Frequency),                  
    SEM = sd(Frequency) / sqrt(n()),          
    n_lawns = n(),                                
    .groups = "drop")
#rearrange order of locations
desired_order<-c("Camarillo", "CSUN", "Antelope Valley")
graphdata$Location<-factor(graphdata$Location, levels=desired_order)

#then use ggplot for plotting
ggplot(graphdata, aes(x=Genotype, y=mean, fill=factor(Location), group=factor(Location),
)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_bar(stat="identity", position="dodge", size=0.6) + #determines the bar width
  geom_errorbar(aes(ymax=mean+SEM, ymin=mean-SEM), stat="identity", position=position_dodge(width=0.9), width=0.1) + #adds error bars
  labs(x="", y="Genotype Frequency", fill="") 
