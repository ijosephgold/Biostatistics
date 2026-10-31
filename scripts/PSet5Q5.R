diversity_data <- read.csv("data/Diversity data.csv")
diversity_data

#location is treated as category
diversity_data$Location <- as.factor(diversity_data$Location)

#diversity remains numeric
str(diversity_data)

#sample variance within each Location
variance_by_location <- tapply(
  diversity_data$Diversity,
  diversity_data$Location,
  var
)
variance_by_location

inspect_df(diversity_data)

#diversity_data$Diversity <- as.factor(diversity_data$Diversity)
lm_model <- lm(Diversity ~ Location, data = diversity_data)
lm_model

anova(lm_model)
summary(lm_model)

lm_model_res <- resid(lm_model)

qqPlot(lm_model_res, "norm")
hist(lm_model_res, xlab = "Diversity residuals", main = "")
boxplot(lm_model_res, ylab = "Diversity residuals")
shapiro.test(lm_model_res) #P > 0.05
skewness(lm_model_res) 
kurtosis(lm_model_res)
#next test equality of variances
bartlett.test(lm_model_res ~ Location, data = diversity_data) #P > 0.05

diversity_data$log_diversity_data <- log(diversity_data$Diversity)

log_lm_model <- lm(log_diversity_data ~ Location, data = diversity_data)
log_lm_model_res <- resid(log_lm_model)
summary(log_lm_model)

qqPlot(log_lm_model_res, "norm")
hist(log_lm_model_res, xlab = "Log Diversity residuals", main = "")
boxplot(log_lm_model_res, ylab = "Log Diversity residuals")
shapiro.test(log_lm_model_res) #P > 0.05
skewness(log_lm_model_res) 
kurtosis(log_lm_model_res)
#next test equality of variances
bartlett.test(log_lm_model_res ~ Location, data = diversity_data) #P > 0.05

diversity_data$sqrt_diversity_data <- sqrt(diversity_data$Diversity)
sqrt_lm_model <- lm(sqrt_diversity_data ~ Location, data = diversity_data)
sqrt_lm_model_res <- resid(sqrt_lm_model)
summary(sqrt_lm_model)

qqPlot(sqrt_lm_model_res, "norm")
hist(sqrt_lm_model_res, xlab = "Sqrt Diversity residuals", main = "")
boxplot(sqrt_lm_model_res, ylab = "Sqrt Diversity residuals")
shapiro.test(sqrt_lm_model_res) #P > 0.05
skewness(sqrt_lm_model_res) 
kurtosis(sqrt_lm_model_res)
#next test equality of variances
bartlett.test(sqrt_lm_model_res ~ Location, data = diversity_data) #P > 0.05

kruskal.test(Diversity ~ Location, data = diversity_data)
bartlett.test(Diversity ~ Location, data = diversity_data)
