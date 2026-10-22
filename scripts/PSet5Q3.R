#H_0: the variance among the gall for mean thorax size is the same for all groups (sigma)^2_gall = 0
#H_1: the variance among the gall for mean thorax size is not the same for all groups (sigma)^2_gall > 0

aphids <- read.csv("data/Aphid data.csv")

head(aphids)
str(aphids)

#make Gall categorical
aphids$Gall <- as.factor(aphids$Gall)

model <- aov(Length ~ Gall, data = aphids)

summary(model)

#because p < 0.05, we reject the null hypothesis and conclude that there is a significant difference in mean thorax size among the different gall groups.
inspect_df(aphids)
MS_gall <- 1.7029100529
MS_within <- 0.0838095238
n <- 4

variance_gall <- (MS_gall - MS_within) / n
variance_environment <- MS_within

variance_gall
variance_environment

#so (sigma)^2_genetic = 0.4048
#(sigma)^2_environmental = 0.0838

#(sigma)^2_total = (sigma)^2_genetic + (sigma)^2_environmental

total_variance <- variance_gall + variance_environment

genetic_percent <-
  variance_gall / total_variance * 100

environment_percent <-
  variance_environment / total_variance * 100

total_variance
genetic_percent
environment_percent

#conclusion: The estimated genetic variance among galls was 0.4048, while the estimated environmental variance within galls was 0.0838. 
#Thus, approximately 82.8% of the phenotypic variance in aphid thorax length was attributed to differences among genetically distinct aphid colonies, while approximately 17.2% was attributed to environmental or individual-level variation within galls.