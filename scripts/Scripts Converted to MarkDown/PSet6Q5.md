```r
load_biostats_packages()
data <- read.csv("data/Hawk data.csv")
view(data)
str(data)
inspect_df(data)

#check assumptions of normality
hist(data$Jawforce)
shapiro.test(data$Jawforce)
qqPlot(data$Jawforce)
assumption_report(data = data)

hist(data$Footforce)
shapiro.test(data$Footforce)
qqPlot(data$Footforce)

cor.test(data$Jawforce, data$Footforce, method = "pearson")

#make a scatterplot
ggplot(data, aes(x=Jawforce, y=Footforce)) +
  geom_point(size = 2) +
  xlab("Footforce") +
  ylab("Jawforce") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=13), 
        axis.text.y=element_text(color="black", size=13), 
        axis.title.x = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())
abline(lm(Footforce ~ Jawforce, data = data), col = "red")

#scatterplot the relationship between jawforce and footforce 
plot(Jawforce ~ Footforce, data = data,
     xlab = "Footforce",
     ylab = "Jawforce",
     main = "Footforce by Jawforce Hawks",
     pch = 19)
abline(lm(Jawforce ~ Footforce, data = data), col = "red")

correlation_workflow(data, x = "Jawforce", y = "Footforce", method = "pearson", plots = TRUE)

```
