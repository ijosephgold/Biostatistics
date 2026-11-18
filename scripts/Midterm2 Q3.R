install_biostats_packages()
load_biostats_packages()
herp_data <- read.csv("data/Herp data.csv")
view(herp_data)
biostat_wizard(herp_data)
biostat_inspect(herp_data)
analysis_start(herp_data, question = "Is habitat area a good predictor of herp species richness? Perform an appropriate test to assess this hypothesis. If the results are significant, report the strength of the effect with an informative metric. You may presume that the data meet the assumptions of normality and equal variances")

lmherp <- lm(SppRichness ~ Hectares, data=herp_data, na.action="na.omit")

lmherpres <- resid(lmherp)
# call the fitted y values
fitted(lmherp) 
#plot the fitted values against the residuals
plot(lmherpres ~ fitted(lmherp), ylab = "Residuals", xlab = "Predicted log y") #it looks better
abline(h=0, col="red") #adds a horizontal line at 0
summary(lmherp)
anova(lmherp)
Anova(lmherp)
cor.test(herp_data$Hectares, herp_data$SppRichness, method="pearson", na.rm=TRUE)
cor.test(herp_data$Hectares, herp_data$SppRichness, method="spearman", na.rm=TRUE)

#make a scatterplot with trendline

linear_regression_workflow(herp_data, SppRichness ~ Hectares)
ggplot(herp_data, aes(x=Hectares, y=SppRichness)) +
  geom_point(size = 2) +
  xlab("Hectares") +
  ylab("Species Richness") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  geom_smooth(method = "lm") #adds regression line with 95% confidence intervals
