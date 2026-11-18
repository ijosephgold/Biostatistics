load_biostats_packages()
data <- read.csv("data/Krat data.csv")
view(data)

#first fit a linear model
datalin <- lm(RatDnsty ~ ShrbCvr + SeedPrd + SnkDnsty, data=data, na.action="na.omit")
mfdatalin <- model.frame(datalin)

datalinres <- resid(datalin)
#termplot(datalin, se = TRUE, partial.resid = FALSE)
#plot model
plot(RatDnsty ~ ShrbCvr, data = data,
     xlab = "Shrub cover", ylab = "Rat density")
abline(lm(data$RatDnsty ~ data$ShrbCvr, data = data), col = "red")
plot(RatDnsty ~ SeedPrd, data = data,
     xlab = "Seed production", ylab = "Rat density")
abline(lm(data$RatDnsty ~ data$SeedPrd, data = data), col = "red")
plot(RatDnsty ~ SnkDnsty, data = data,
     xlab = "Snake density", ylab = "Rat density")
abline(lm(data$RatDnsty ~ data$SnkDnsty, data = data), col = "red")

#check the assumptions
#recall, we want to look at the residuals
#first test normality
datalinres <- residuals(datalin)
#qqPlot(datalinres, main = "Q-Q plot of residuals", ylab = "Residuals")

#check the assumptions
#recall, we want to look at the residuals
#first test normality
#model1res <- resid(model1)
#qqp(model1res, "norm")
#hist(model1res)
#shapiro.test(model1res) #P < 0.05. It's a bit non-normal
#test for homogeneity of variance
#call the fitted y values
#fitted(model1) #this is the predicted value of y, for each value of x
#plot the fitted (predicted) values against the residuals
#plot(model1res~fitted(model1), ylab = "Residuals", xlab = "Predicted y") #it looks a bit cone shaped


shrub_model <- lm(RatDnsty ~ ShrbCvr, data = data)
qqp(residuals(shrub_model), "norm")
hist(residuals(shrub_model))
shapiro.test(residuals(shrub_model)) #P > 0.05. This is better
fitted(shrub_model)
plot(residuals(shrub_model) ~ fitted(shrub_model), ylab = "Residual", xlab = "Predicted y")

ggplot(data, aes(x = ShrbCvr, y = RatDnsty)) +
  geom_point() +
  xlab("Shrub cover") +
  ylab("Rat density") +
  geom_smooth(method = "lm", formula = y ~ x, 
              se = TRUE, level = 0.95) +
  labs(x = "Shrub cover", y = "Rat density")

seed_model <- lm(RatDnsty ~ SeedPrd, data = data)
qqp(residuals(seed_model), "norm")
hist(residuals(seed_model))
shapiro.test(residuals(seed_model)) #P > 0.05. This is better
fitted(seed_model)
plot(residuals(seed_model) ~ fitted(seed_model), ylab = "Residual", xlab = "Predicted y")

ggplot(data, aes(x = SeedPrd, y = RatDnsty)) +
  geom_point() +
  xlab("Seed production") +
  ylab("Rat density") +
  geom_smooth(method = "lm", formula = y ~ x, 
              se = TRUE, level = 0.95) +
  labs(x = "Seed production", y = "Rat density")

snake_model <- lm(RatDnsty ~ SnkDnsty, data = data)
qqp(residuals(snake_model), "norm")
hist(residuals(snake_model))
shapiro.test(residuals(snake_model)) #P > 0.05. This is better
fitted(snake_model)
plot(residuals(snake_model) ~ fitted(snake_model), ylab = "Residual", xlab = "Predicted y")

ggplot(data, aes(x = SnkDnsty, y = RatDnsty)) +
  geom_point() +
  xlab("Snake density") +
  ylab("Rat density") +
  geom_smooth(method = "lm", formula = y ~ x, 
              se = TRUE, level = 0.95) +
  labs(x = "Snake density", y = "Rat density")


#plot(datalinres ~ data$ShrbCvr, xlab = "Shrub cover", ylab = "Residuals")


#try a log transformation
data$logRatDnsty <- log(data$RatDnsty)
data$logShrbCvr <- log(data$ShrbCvr)
data$logSeedPrd <- log(data$SeedPrd)
data$logSnkDnsty <- log(data$SnkDnsty)
data

#create a new linear model with the log data
logmodel<-lm(logNumberSpp~logBiomass, data=riffledata, na.action="na.omit")

#check the assumptions of the log model
#test normality
logmodelres <- resid(logmodel)
qqp(logmodelres, "norm")
hist(logmodelres)
shapiro.test(logmodelres) #P > 0.05. This is better
#test for homogeneity of variance
# call the fitted y values
fitted(logmodel) 
#plot the fitted values against the residuals
plot(logmodelres~fitted(logmodel), ylab = "Residuals", xlab = "Predicted log y") #it looks better

#let's perfrom regression with the transformed data
summary(logmodel) #notice it also does a t-test (but t is calculated differently than in correlation)
#make an ANOVA table
anova(logmodel)


#let's perfrom regression with the transformed data
summary(logmodel) #notice it also does a t-test (but t is calculated differently than in correlation)
#make an ANOVA table
anova(logmodel)

#make a scatterplot with trendline
ggplot(riffledata, aes(x=logBiomass, y=logNumberSpp)) +
  geom_point(size = 2) +
  xlab("Log (Biomass)") +
  ylab("Log (Species Richness)") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals
