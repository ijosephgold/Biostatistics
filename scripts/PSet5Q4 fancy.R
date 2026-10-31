#Repeated measures ANOVA
#H_0: The algal symbionts of these Polyps can not acclimate to elevated temperatures
#H_A: At least one can acclimate to elevated temperatures
#day is a repeated measure, temperature is a seperate 3 group factor, and coral is a random factor. The repeated-measures ANOVA will test for the effects of temperature, day, and their interaction on FvFm while accounting for the repeated measurements on the same coral.
#10 corals per temperature treatment, measured over 6 days. The data is balanced and there are no missing values or duplicates for the same coral on the same day.

#load libraries
library(tidyverse)
library(rstatix) #for identifying outliers
library(afex) #for ANOVA
library(ggpubr) #for plotting
library(emmeans)
library(ggplot2)
load_biostats_packages(install_missing = TRUE)

corals <- read.csv("data/Coral data.csv") %>%
  mutate(
    temperature = factor(
    temperature, 
    levels = c("ambient", "hot", "warm")
    ),
    coral = factor(coral),
    day = factor(day, levels = 1:6)
  )

# Check that the repeated-measures structure is intact
stopifnot(!anyNA(corals))

if(anyDuplicated(corals[c("coral", "day")])) {
  stop("Duplicate measurements found for the same coral on the same day.")
}

temperature_per_coral <- tapply(
  as.character(corals$temperature), corals$coral,
  function(x) length(unique(x))
  ) 

if(any(temperature_per_coral != 1)) {
  stop("At least one coral has measurements at multiple temperatures, which violates the repeated-measures assumption.")
}

cell_counts <- table(corals$temperature, corals$day)
print(cell_counts)

if(length(unique(as.vector(cell_counts))) !=1) {
  stop("The temperature-by-day groups are not balanced")
}

#Repeated-measures ANOVA:
#temperature = between-coral factor
#day = repeated within-coral factor
rm_anova <- aov_ez(
  id = "coral",
  dv = "FvFm",
  data = corals,
  between = "temperature",
  within = "day",
  anova_table = list(correction = "GG", es = "pes")
)

#Show ANOVA with Greenhouse-Geisser correction
nice(rm_anova, correction = "GG", es = "pes")

#Follow up interaction
#Compare temperatures seperately on each day
emm_temperature_by_day  <- emmeans(rm_anova, ~ temperature | day)
pairs(emm_temperature_by_day, adjust = "tukey")

#Repeated measure ANOVA with "Greenhouse-Geisser correction; needed to adjust day and temperature by day tests for potential violations of "Sphericity"
#the following comparison help show which days or temperature groups differ; the overall interaction alone does not identify every specific difference.
#the graph plots the mean FvFm for each temperature group across the 6 days, with error bars representing the standard error of the mean. The different line types and colors help distinguish between the temperature treatments.
graphdata <- corals %>%
  group_by(temperature, day) %>%
  summarise(
    mean = mean(FvFm),
    se = sd(FvFm) / sqrt(n())
  ) %>%
  mutate(day_num = as.numeric(as.character(day)))

fvfm_plot <- ggplot(
  graphdata,
  aes(
    x = day_num,
    y = mean,
    color = temperature,
    linetype = temperature, 
    group = temperature
  )
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.12,
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = c(ambient = "blue", "hot" = "red", "warm" = "orange")
  ) +
  scale_linetype_manual(
    values = c(ambient = "solid", hot = "dotted", warm = "dashed")
  ) +
  scale_x_continuous(breaks = 1:6) +
  coord_cartesian(ylim = c(0.30, 0.46)) +
  labs(
    x = "Day",
    y = "PS efficiency (Fv/Fm)",
    color = "Temperature",
    linetype = "Temperature"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
fvfm_plot
#y-axis is limited to 0.30-0.46 to better visualize the differences between temperature treatments, as the FvFm values are relatively close together and this range allows for a clearer comparison of the trends over time.
#the error bars represent the standard error of the mean for each temperature group on each day, providing a visual indication of the variability and reliability of the mean estimates.

#Results
#temperature by day interaction: F(5.37, 72.49) ≈ 50.38, p < 0.001 (Greenhouse-Geisser corrected). Fv/Fm patern differs over time among btemperature groups
#Day: F(2.68, 72.49) ≈ 26.98, p < 0.001 Fv/Fm changes across days when averaged over temperatures
#temperature averaged across days: F(2, 27) ≈ 14.79, p < 0.001, this overall comaprison is not significant at a 0.05 level
#Interpretation: the key result is the Interaction, where the temperatures do not show the same pattern throughout the 6 days. This supports a temperature dependent change over time. 
#Conclusion: Repeated Measure ANOVA found significant affects of temperature, day, and their interaction (all p < 0.05). interaction means Fv/Fm changed differently over time across temperatures: initial decreases followed by increases over time, suggestingng they are acclimating! 