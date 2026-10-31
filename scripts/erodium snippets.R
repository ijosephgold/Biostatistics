library(tidyverse)
library(car)
library(rstatix)

erodiumdata <- read.csv("Data/Erodium data.csv") %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("control", "shade", "wind", "3XH2O")
    ),
    plant_size = Width * Height,
    log_size = log(plant_size)
  )

# Check residual normality using the log-transformed model
log_model <- lm(log_size ~ Treatment, data = erodiumdata)
qqPlot(residuals(log_model))
shapiro.test(residuals(log_model))

# Test the treatment effect; Welch's test allows unequal variances
oneway.test(log_size ~ Treatment, data = erodiumdata, var.equal = FALSE)

# Compare each treatment pair; Games–Howell allows unequal variances
games_howell_test(erodiumdata, log_size ~ Treatment)