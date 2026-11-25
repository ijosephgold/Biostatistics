#Rey et al. 2019 post fire bird data

#load libraries
library(tidyverse)
library(emmeans)
library(ggplot2)

#load data
firedata <- read.csv("Data/Reyetal.csv")
firedata

#first I want to perform a two-way ANOVA with site and forest state as x variables and species richness as the response

#make x variables factors
firedata$Site <- factor(firedata$Site)
firedata$foreststate <- factor(firedata$foreststate)

#create the two-way GLM
model1 <- lm(n_spec~Site*foreststate, data = firedata)

#check residuals
model1res <- resid(model1)
qqPlot(model1res)
shapiro.test(model1res) #not great

#try log transformation
firedata$log_nspec <- log(firedata$n_spec)
model2 <- lm(log_nspec~Site*foreststate, data = firedata)

#check residuals
model2res <- resid(model2)
qqPlot(model2res)
shapiro.test(model2res) #not great

#try sqrt transformation
firedata$sqrt_nspec <- sqrt(firedata$n_spec)
model3 <- lm(sqrt_nspec~Site*foreststate, data = firedata)

#check residuals
model3res <- resid(model3)
qqPlot(model3res)
shapiro.test(model3res) #not great

#just run the two way ANOVA with the raw data
anova(model1) #the interaction term isn't significant
#rerun the model without it
model4 <- lm(n_spec~Site + foreststate, data = firedata)
anova(model4)

#make barplot
graphdata <- as.data.frame(emmeans(model1, ~ Site*foreststate))
graphdata 
ggplot(graphdata, aes(x=Site, y=emmean, fill=factor(foreststate), group=factor(foreststate))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_bar(stat="identity", position="dodge", size=0.6) + #determines the bar width
  geom_errorbar(aes(ymax=emmean+SE, ymin=emmean-SE), stat="identity", position=position_dodge(width=0.9), width=0.1) + #adds error bars
  labs(x="", y="Species richness", fill="") #labels the x and y axes


#Next I want to test territory area as the response against site and forest state
model5 <- lm(terr_10ha~Site*foreststate, data = firedata)

#check residuals
model5res <- resid(model5)
qqPlot(model5res)
shapiro.test(model5res) #ok
#test for homogeneity of variance
# call the fitted y values
fitted(model5) 
#plot the fitted values against the residuals
plot(model5res~fitted(model5), ylab = "Residuals", xlab = "Predicted log y") #ew

#try log transformation
firedata$logterr <- log(firedata$terr_10ha)
model6 <- lm(logterr~Site*foreststate, data = firedata)

#check residuals
model6res <- resid(model6)
qqPlot(model6res)
shapiro.test(model6res) #not ok
#test for homogeneity of variance
# call the fitted y values
fitted(model6) 
#plot the fitted values against the residuals
plot(model5res~fitted(model5), ylab = "Residuals", xlab = "Predicted log y") #ew
#this doesn't really help

#I'll just analyze the raw data
anova(model5)
#the interaction is significant
#we could drop it and rerun the analysis or do separate t-tests
#terrain area vs. site
t.test(terr_10ha ~ Site, var.equal = FALSE, data = firedata)
#terrain area vs. forest state (burned or unburned)
t.test(terr_10ha ~ foreststate, var.equal = FALSE, data = firedata)

#make barplot
graphdata <- as.data.frame(emmeans(model5, ~ Site*foreststate))
graphdata 
ggplot(graphdata, aes(x=Site, y=emmean, fill=factor(foreststate), group=factor(foreststate))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_bar(stat="identity", position="dodge", size=0.6) + #determines the bar width
  geom_errorbar(aes(ymax=emmean+SE, ymin=emmean-SE), stat="identity", position=position_dodge(width=0.9), width=0.1) + #adds error bars
  labs(x="", y="Territory Ara (Ha)", fill="") #labels the x and y axes


#next I want to test species richness vs. territory area
#I'm ignoring the assumptions and just using a non-parametric test
cor.test(firedata$n_spec, firedata$terr_10ha, method = "spearman")

#make a scatterplot
ggplot(firedata, aes(x=n_spec, y=terr_10ha)) +
  geom_point(size = 2) +
  xlab("Species richness") +
  ylab("Average territory/10ha") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())


#I also want to know if time since fire has a linear relationship with species richness and territory area
#time vs. richness
cor.test(firedata$years_postfire, firedata$n_spec, method ="spearman")
#time vs. territory area
cor.test(firedata$years_postfire, firedata$terr_10ha, method ="spearman")

#make scatterplots
#time vs. richness
ggplot(firedata, aes(x=years_postfire, y=n_spec)) +
  geom_point(size = 2) +
  xlab("Years since fire") +
  ylab("Species richness") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
#this looks weird but ok

#time vs. territory area
ggplot(firedata, aes(x=years_postfire, y=terr_10ha)) +
  geom_point(size = 2) +
  xlab("Years since fire") +
  ylab("Territory area (Ha)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
