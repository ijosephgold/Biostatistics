load_biostats_packages(install_missing = FALSE)

#H_0: Beak width and beak depth have no linear association in the population; the correlation is zero (ρ = 0)
#H_A: Beak width and beak depth have a linear association; the correlation is not zero (ρ ≠ 0)
penguins <- read.csv("data/penguin data - Sheet1.csv")
View(penguins)

cor.test(penguins$Beak.width, penguins$Beak.depth, method = "pearson")

plot(penguins$Beak.width, penguins$Beak.depth,
     xlab = "Beak width (mm)",
     ylab = "Beak depth (mm)",
     main = "Beak width and depth in African penguins",
     pch = 19)
abline(lm(penguins$Beak.depth ~ penguins$Beak.width), col = "blue", lwd = 2)

model <- lm(penguins$Beak.depth ~ penguins$Beak.width)

plot(model, which = 1)  
#Plot 1: Residuals vs. Fitted. residual is (depth y_i) - (depth y-hat), looking for residuals scattered random around 0, curve ~suggest relationship isn’t linear; funnel shape ~suggest residual spread changes across fitted values
plot(model, which = 2)  
#Plot 2: Normal Q-Q. compares residuals with what we’d expect if they were normally distributed. Points roughly following line suggest normality is reasonable, visually complements Shapiro–Wilk test ran

shapiro.test(resid(model))


#penguins$differences <- penguins$Beak.depth - penguins$Beak.width
#var(penguins$differences)
#qqPlot(penguins$differences)
#shapiro.test(penguins$differences)
