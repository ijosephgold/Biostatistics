data <- read.csv("data/Crayfish data.csv")
view(data)

# 1. Fit Binary Logistic Regression GLM (Binomial family, logit link)
mod <- glm(Mated ~ Claw.diff, data = data, family = binomial)
summary(mod)                                # Check Wald z-value and p-value (Pr(>|z|))

# 2. Calculate Numeric Factor by which 1-unit increase in X increases Odds
exp(coef(mod))                              # exp(b1) = Odds Ratio multiplier
exp(confint(mod))                           # 95% Confidence Interval for Odds Ratio

# 3. Assess Strength of Association (McFadden's Pseudo-R2 & LRT Chi-Sq)
1 - (mod$deviance / mod$null.deviance)      # McFadden's Pseudo-R2 (proportion deviance explained)
anova(mod, test = 'Chisq')                  # Likelihood Ratio Test (Drop-in-Deviance G2)

# 4. Plot Logistic Regression Curve with Raw Points
library(ggplot2)
ggplot(data, aes(x = Claw.diff, y = Mated)) +
  geom_jitter(height = 0, width = 0, alpha = 0.6) +
  stat_smooth(method = 'glm', method.args = list(family = 'binomial'), se = TRUE, color = 'red') +
  theme_classic()
