#Repeated measures ANOVA

#load libraries
library(tidyverse)
library(rstatix) #for identifying outliers
library(afex) #for ANOVA
library(ggpubr) #for plotting
library(emmeans)
library(ggplot2)

#load data
biofeedback <- read.csv("Data/Stress data.csv")
biofeedback
inspect

#check model assumptions
qqPlot(biofeedback$AverageBPM, ylab = "Average BPM")
shapiro.test(biofeedback$AverageBPM) #P = 0.1, data has normal distribution
biofeedback %>%
  group_by(Test) %>%
  identify_outliers(AverageBPM) #no outliers

#perform repeated measures ANOVA
#build linear model
rm_anova <- aov_ez(data = biofeedback, dv = "AverageBPM", id = "Subject", within = "Test")
#view the ANOVA & HSD results
summary(rm_anova)
em_results <- emmeans(rm_anova, ~ Test)
pairs(em_results, adjust = "tukey") #Tukey's HSD

#make a boxplot
bxp <- ggboxplot(biofeedback, x = "Test", y = "AverageBPM", add = "point", notch = FALSE)
bxp

#try a bar plot
#calculate means and variance
graphdata <- biofeedback %>%
  group_by(Test) %>%
  summarize(mean=mean(AverageBPM, na.rm=TRUE), se=sd(AverageBPM, na.rm=TRUE)/sqrt(length(na.omit(AverageBPM))))
graphdata

#ggplot
ggplot(data=graphdata, aes(x=Test, y=mean, fill=Test)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold", color="black", size=16), 
        axis.text.y=element_text(face="bold", color="black", size=13), 
        axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_bar(colour="black", fill="DodgerBlue", width=0.7, stat="identity") + 
  guides(fill=FALSE) + 
  ylab("Average Heart Rate (BPM)") +
  xlab("") +
  #scale_x_discrete(labels=c("Baseline" = "Baseline", "1" = " Test 1", "2" = "Test 2", "3" = "Test 3")) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), stat="identity", 
                position=position_dodge(width=0.9), width=0.1) +
  coord_cartesian(ylim = c(75, 100))
