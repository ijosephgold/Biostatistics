```r
#This asks whether each factor is related to rat density while accounting for the other two factors

load_biostats_packages()
data <- read.csv("data/Krat data.csv")
view(data)
analysis_guide(data, outcome = "RatDnsty", predictors = c("ShrbCvr", "SeedPrd", "SnkDnsty"))
#Fit one model with all three predictors
datalin <- lm(
  RatDnsty ~ ShrbCvr + SeedPrd + SnkDnsty, #~ means â€œas a function ofâ€
  data = data,
  na.action = "na.omit"
)

#Save this model's residuals
datalinres <- resid(datalin) #resid(datalin) extracts the differences between the observed rat densities and the values predicted by this model

#Check whether the residuals are reasonably normal
qqp(datalinres, "norm")
hist(datalinres)
shapiro.test(datalinres)

#Check residual spread against predicted rat density
#plot(datalinres ~ fitted(datalin)) plots residuals vertically against predicted values horizontally. Look for a roughly even cloud around zero, without a curve or funnel shape
plot(datalinres ~ fitted(datalin), #fitted(datalin) extracts the modelâ€™s predicted rat densities
     ylab = "Residuals",
     xlab = "Predicted rat density")
abline(h = 0, col = "red", lty = 2)

#Examine the model's results
summary(datalin)

termplot(datalin, se = TRUE, partial.resid = FALSE)

```
