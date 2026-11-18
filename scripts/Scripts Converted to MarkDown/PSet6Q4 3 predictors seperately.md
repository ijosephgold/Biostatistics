```r
load_biostats_packages()
data <- read.csv("data/Krat data.csv")
View(data)
inspect_df(data)

#fit three separate models
#model for Shrub Cover
shrub_model <- lm(RatDnsty ~ ShrbCvr, data = data,
                  na.action = "na.omit")

#model for Seed Production
seed_model <- lm(RatDnsty ~ SeedPrd, data = data,
                 na.action = "na.omit")

#model for Snake Density
snake_model <- lm(RatDnsty ~ SnkDnsty, data = data,
                  na.action = "na.omit")

#check Shrub Cover model for assumptions 
shrub_modelres <- resid(shrub_model)
qqp(shrub_modelres, "norm", main = "qqPlot Shrub Cover Model Residuals")
hist(shrub_modelres, main = "Histogram Shrub Cover Model Residuals", xlab = "Residuals")
shapiro.test(shrub_modelres)

#plot(shrub_modelres ~ fitted(shrub_model),
#     ylab = "Residuals",
#     xlab = "Predicted rat density")
#abline(h = 0, col = "red", lty = 2)

#make a scatterplot with trendline
ggplot(data, aes(x=RatDnsty, y=ShrbCvr)) +
  geom_point(size = 2) +
  xlab("Rat Density") +
  ylab("Shrub Cover") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals

#check Seed Production model for assumptions 
seed_modelres <- resid(seed_model)
qqp(seed_modelres, "norm", main = "qqPlot Seed Production Model Residuals")
hist(seed_modelres, main = "Histogram Seed Production Model Residuals", xlab = "Residuals")
shapiro.test(seed_modelres)

#plot(seed_modelres ~ fitted(seed_model),
#     ylab = "Residuals",
#     xlab = "Predicted rat density")
#abline(h = 0, col = "red", lty = 2)

ggplot(data, aes(x=RatDnsty, y=SeedPrd)) +
  geom_point(size = 2) +
  xlab("Rat Density") +
  ylab("Seed Production") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals

#check Snake Density model for assumptions 
snake_modelres <- resid(snake_model)
qqp(snake_modelres, "norm", main = "qqPlot Snake Density Model Residuals")
hist(snake_modelres, main = "Histogram Snake Density Model Residuals", xlab = "Residuals")
shapiro.test(snake_modelres)

#plot(snake_modelres ~ fitted(snake_model),
#     ylab = "Residuals",
#     xlab = "Predicted rat density")
#abline(h = 0, col = "red", lty = 2)

ggplot(data, aes(x=RatDnsty, y=SnkDnsty)) +
  geom_point(size = 2) +
  xlab("Rat Density") +
  ylab("Snake Density") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals

#assumptions for models relatively satisfied, examine results
summary(shrub_model)
anova(shrub_model)
summary(seed_model)
anova(seed_model)
summary(snake_model)
anova(snake_model)

assumption_report(data = shrub_modelres)
assumption_report(data = seed_modelres)
assumption_report(data = snake_modelres)

```
