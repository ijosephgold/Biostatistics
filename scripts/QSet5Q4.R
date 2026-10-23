corals <- read.csv("data/Coral data.csv")
inspect_df(corals)

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
