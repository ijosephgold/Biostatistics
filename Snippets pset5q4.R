

add_transformations(corals, "FvFm")









vcorals$temperature <- as.factor(corals$temperature)
lm <- lm(FvFm ~ temperature, data = corals) #takes the data as lm(y~x)
lm
anova(lm) #provides the ANOVA table
summary(lm) #provides extra details about the model


#check the assumptions
#recall, we want to look at the residuals to assess normality
lmres <- resid(lm)
qqplot(lmres, "norm") 
hist(lmres, xlab = "FvFm residuals", main = "")
boxplot(lmres, ylab = "FvFm residuals")
shapiro.test(lmres) #P > 0.05
skewness(lmres) 
kurtosis(lmres)
#next test equality of variances
bartlett.test(FvFm~temperature, data=corals) #P > 0.05


anova_assumption_report(data = corals, formula = FvFm ~ temperature)

corals <- add_transformations(corals, "FvFm")
corals

anova_assumption_report(data = corals, formula = FvFm_log ~ temperature)



results <- transformation_workflow(corals, variable = "FvFm")
corals <- results$data
names(corals)
results <- transformation_workflow(
  corals,
  variable = "FvFm",
  save_to = "data/Coral_data_transformed.csv"
)

corals <- results$data
inspect_df(corals)



#now to tackle Aclimation
#creating one row per coral, with one column for each day
change_data <- corals %>%
  mutate(day = as.character(day)) %>%
  select(coral, temperature, day, FvFm) %>%
  pivot_wider(
    names_from = day,
    values_from = FvFm,
    names_prefix = "day_"
  ) %>%
  mutate(
    change = day_6 - day_1,
    temperature = factor(temperature, levels = c("ambient", "hot", "warm"))
  )
print(n = 30, change_data)

#Checking average change in each group
change_data %>%
  group_by(temperature) %>%
  summarize(
    mean_change = mean(change),
    sd_change = sd(change),
    n = n()
  )

# Compare hot with ambient using a Welch t-test
hot_test <- t.test(
  change ~ temperature,
  data = droplevels(
    filter(change_data, temperature %in% c("ambient", "hot"))
  ),
  var.equal = FALSE
)

# Compare warm with ambient
warm_test <- t.test(
  change ~ temperature,
  data = droplevels(
    filter(change_data, temperature %in% c("ambient", "warm"))
  ),
  var.equal = FALSE
)

hot_test
warm_test

# Adjust the two p-values for multiple comparisons
p.adjust(
  c(hot_vs_ambient = hot_test$p.value,
    warm_vs_ambient = warm_test$p.value),
  method = "holm"
)

#creating one row per coral, then calculate its day 6 minus day 1 change
change_data <- corals %>%
  mutate(day = as.character(day)) %>%
  select(coral, temperature, day, FvFm) %>%
  pivot_wider(
    names_from = day,
    values_from = FvFm,
    names_prefix = "day_"
  ) %>%
  mutate(change = day_6 - day_1)

# Compare each heated group with ambient using Welch tests
ambient_change <- filter(change_data, temperature == "ambient")$change
hot_change <- filter(change_data, temperature == "hot")$change
warm_change <- filter(change_data, temperature == "warm")$change

hot_test <- t.test(hot_change, ambient_change, var.equal = FALSE,
                   conf.level = 0.975)
warm_test <- t.test(warm_change, ambient_change, var.equal = FALSE,
                    conf.level = 0.975)

# Adjust the two p-values for multiple comparisons
adjusted_p <- p.adjust(
  c(hot_test$p.value, warm_test$p.value),
  method = "holm"
)

contrast_data <- data.frame(
  comparison = c("Hot − ambient", "Warm − ambient"),
  estimate = c(
    hot_test$estimate[1] - hot_test$estimate[2],
    warm_test$estimate[1] - warm_test$estimate[2]
  ),
  lower = c(hot_test$conf.int[1], warm_test$conf.int[1]),
  upper = c(hot_test$conf.int[2], warm_test$conf.int[2]),
  p_holm = adjusted_p
)

contrast_data

# Plot each difference in change with its confidence interval
ggplot(contrast_data, aes(x = estimate, y = comparison)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_segment(
    aes(x = lower, xend = upper, yend = comparison),
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  labs(
    x = "Difference in day 1–to–day 6 change (Fv/Fm units)",
    y = NULL
  ) +
  theme_classic()

