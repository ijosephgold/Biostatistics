#ANOVA

#load libraries
library(tidyverse)
library(car) #for QQ-plot
library(moments) #for skewness and kurtosis)
library(agricolae) #for Tukey's HSD
library(ggplot2)
library(lme4) #for random factors
library(vegan) #for perMANOVA

#Model I ANOVA
#test California Sagebrush heights in different soil types

#load data
sagedata <- read.csv("Data/Sagebrush data.csv")
sagedata

#fit data to a linear model with soil = x and height = y
sagedata$Soil <- as.factor(sagedata$Soil)
lm <- lm(Height~Soil, data=sagedata) #takes the data as lm(y~x)
lm
anova(lm) #provides the ANOVA table
summary(lm) #provides extra details about the model

#check the assumptions
#recall, we want to look at the residuals to assess normality
lmres <- resid(lm)
qqPlot(lmres, "norm")
hist(lmres, xlab = "Height residuals", main = "")
boxplot(lmres, ylab = "Height residuals")
shapiro.test(lmres) #P > 0.05
skewness(lmres) 
kurtosis(lmres)
#next test equality of variances
bartlett.test(Height~Soil, data=sagedata) #P > 0.05
#or make a boxplot
Sterile <- sagedata %>%
  filter(Soil == "Sterile")
SandDune <- sagedata %>%
  filter(Soil == "Sand Dune")
Mountains <- sagedata %>%
  filter(Soil == "Mountains")
Chaparral <- sagedata %>%
  filter(Soil =="Chaparral")
boxplot(Sterile$Height, SandDune$Height, Mountains$Height, Chaparral$Height,
        names = c("Sterile", "Sand Dune", "Mountains", "Chaparral"), 
        xlab = "Soil type", ylab = "Height (cm)") #spread is evenish among samples

#the levels of the x variable (soil type) are fixed
#perform Tukey's HSD post-hoc comparisons
TukeyHSD <- HSD.test(lm, "Soil", console=TRUE)
#it's also possible to get the p-values for each comparison
model <- aov(Height~Soil, data=sagedata)
Tukeymodel <- TukeyHSD(model)
Tukeymodel
plot(Tukeymodel) #can visualize the confidence intervals

#prepare data for plotting
graphdata <- sagedata %>%
  group_by(Soil) %>%
  summarize(mean=mean(Height, na.rm=TRUE), se=sd(Height, na.rm=TRUE)/sqrt(length(na.omit(Height))))
graphdata
#rearrange order of the levels on the x-axis
graphdata$Soil <- factor(graphdata$Soil, levels = c("Sterile", "Sand Dune", "Mountains", "Chaparral"))

#make a barplot with ggplot()
ggplot(data=graphdata, aes(x=Soil, y=mean, fill=Soil)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(color="black", size=14), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=16, face="bold"), 
        axis.title.y = element_text(color="black", size=16, face="bold"),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_bar(colour="black", fill="grey", width=0.8, stat="identity") + 
  guides(fill=FALSE) + 
  ylab("Height (cm)") +
  xlab("Soil Type") +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), stat="identity", 
        position=position_dodge(width=0.9), width=0.1) +
  geom_text(aes(label=c("a", "bc", "b", "c"), y = mean + 1))
#note that the labels are assigned in the order of levels in "graphdata"


#Model II ANOVA
#test different coral genets for polyp area

#load data
genets <- read.csv("Data/Genet data.csv")
genets

#create a univariate ANOVA with genet as a random effect
genets$Genet <- as.factor(genets$Genet)
#lmer() is like lm() but can take random effects
genetmodel <- lmer(Area~1 + (1|Genet), data = genets) #the code 1|Genet sets genet as a random factor
summary(genetmodel) #gives raw variance components

#check assumptions
hist(resid(genetmodel), breaks = 20, xlab = "Area residuals", main = "")
qqPlot(resid(genetmodel))
shapiro.test(resid(genetmodel)) #all okay
#make boxplots of each genet
genet1 <- genets %>%
  filter(Genet == "Genet 1")
genet2 <- genets %>%
  filter(Genet == "Genet 2")
genet3 <- genets %>%
  filter(Genet == "Genet 3")
#groups should have similar dispersion
boxplot(genet1$Area, genet2$Area, genet3$Area, names = c("1", "2", "3"), 
        xlab = "Genet", ylab = "Polyp area") #not exactly perfect
#can also look at a boxplot of pooled sample residuals
boxplot(resid(genetmodel), ylab = "Area residuals") #not too bad
bartlett.test(Area~Genet, data = genets) #not too good
#try a square root transformation
genets$sqrtarea <-sqrt(genets$Area)
genets

#make a glm with sqrt data
sqrtgenetmodel <- lmer(sqrtarea~1 + (1|Genet), data = genets) #the code 1|Genet sets genet as a random factor
summary(sqrtgenetmodel) #gives raw variance components

#check assumptions
hist(resid(sqrtgenetmodel), breaks = 20, xlab = "Sqrt area residuals", main = "")
qqPlot(resid(sqrtgenetmodel))
shapiro.test(resid(sqrtgenetmodel)) #all okay
#make boxplots of each genet
genet1 <- genets %>%
  filter(Genet == "Genet 1")
genet2 <- genets %>%
  filter(Genet == "Genet 2")
genet3 <- genets %>%
  filter(Genet == "Genet 3")
boxplot(genet1$sqrtarea, genet2$sqrtarea, genet3$sqrtarea, names = c("1", "2", "3"), 
        xlab = "Genet", ylab = "Sqrt (polyp area)")
boxplot(resid(sqrtgenetmodel), ylab = "Sqrt area residuals") #can also look at a boxplot of residuals
bartlett.test(sqrtarea~Genet, data = genets) #better but p is still < 0.05

#create a model without genet
genetmodel0 <- lm(sqrtarea ~ 1, data = genets) 
#perform a log likelihood ratio test
anova(sqrtgenetmodel, genetmodel0) #P < 0.05 -- genet has a significant effect
#can also calculate chi-square using the log likelihood scores
X2 <- -2*(89.862-106.995)
1-pchisq(X2, 1) #generate p-value

#calculate variance components
vc <- as.data.frame(VarCorr(sqrtgenetmodel)) 
ICC <- vc$vcov[1] / sum(vc$vcov)
ICC #Intraclass correlation coefficient (%added variance among groups)
varwithin <- 1-ICC #% variance within groups
varwithin

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Area ~ Genet, data = genets) 
#good but really only shows that genet dispersion differs

#perMANOVA
set.seed(123) #makes the permutation output reproducible
adonis2(dist(genets$Area) ~ genets$Genet, permutations = 999) #performs test
#check homogeneity of dispersion
dist_mat <- dist(genets$Area, method = "euclidean") #calculate distance matrix
bd <- betadisper(dist_mat, genets$Genet) #calculate beta dispersion
anova(bd)
permutest(bd)
#tests of homogeneity of dispersion are significant, but perMANOVA isn't detecting an effect of genet
#these perMANOVA results are driven by unequal within-genet variance rather than differences in means
#mixed model ANOVA with sqrt transformed data is the best choice for this analysis