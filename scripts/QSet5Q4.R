#Repeated measures ANOVA
#H_0: The algal symbionts of these Polyps can not acclimate to elevated temperatures
#H_A: At least one can acclimate to elevated temperatures

#load libraries
library(tidyverse)
library(rstatix) #for identifying outliers
library(afex) #for ANOVA
library(ggpubr) #for plotting
library(emmeans)
library(ggplot2)
load_biostats_packages(install_missing = TRUE)

corals <- read.csv("data/Coral data.csv")
inspect_df(corals)
corals

#check model assumptions
qqPlot(corals$FvFm, ylab = "FvFm")
shapiro.test(corals$FvFm) #P = 0.1, data has normal distribution
bartlett.test(FvFm ~ day, data = corals) #P = 0.1, data has equal variance
corals %>%
  group_by(day) %>%
  identify_outliers(FvFm) #no outliers

#perform repeated measures ANOVA
#build linear model
rm_anova <- aov_ez(data = corals, dv = "FvFm", id = "coral", within = "day")
rm_res <- resid(rm_anova, append = TRUE)
hist(rm_res$.residuals, breaks = 20, xlab = "Residuals FvFm", ylab = "Frequency", main = "")
qqPlot(rm_res$.residuals, ylab = "FvFm residuals")
shapiro.test(rm_res$.residuals) #P = 0.1, data has normal distribution
skewness(rm_res$.residuals)
kurtosis(rm_res$.residuals)
bartlett.test(FvFm ~ day, data = corals)
bartlett.test(rm_res$.residuals ~ corals$day, data = corals)
boxplot(rm_res$.residuals, ylab = "FvFm residuals")
var(rm_res$.residuals)
var(corals$FvFm)
inspect_df(rm_res)

#view the ANOVA & HSD results
summary(rm_anova, append = TRUE)
em_results <- emmeans(rm_anova, ~ day)
pairs(em_results, adjust = "tukey") #Tukey's HSD

#make a boxplot
bxp <- ggboxplot(corals, x = "day", y = "FvFm", add = "point", notch = TRUE)
bxp

#try a bar plot
#calculate means and variance
graphdata <- corals %>%
  group_by(day) %>%
  summarize(mean = mean(FvFm, na.rm = TRUE), se = sd(FvFm, na.rm = TRUE) / sqrt(length(na.omit(FvFm))))
graphdata

#ggplot
ggplot(data = graphdata, aes(x = Treatment, y = mean, fill = Treatment)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 16), 
        axis.text.y = element_text(face = "bold", color = "black", size = 13), 
        axis.title.x = element_text(color = "black", size = 14, face = "bold"), 
        axis.title.y = element_text(color = "black", size = 14, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_bar(colour = "black", fill = "DodgerBlue", width = 0.7, stat = "identity") + 
  guides(fill = FALSE) + 
  ylab("Average FvFm") +
  xlab("") +
  scale_x_discrete(labels = c("baseline" = "day 1", "1" = "day 2", "2" = "day 3", "3" = "day 4", "4" = "day 5", "5" = "day 6")) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), stat = "identity", 
                position = position_dodge(width = 0.9), width = 0.1) +
  coord_cartesian(ylim = c(75, 100))

