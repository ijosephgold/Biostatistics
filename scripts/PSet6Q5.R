load_biostats_packages()
data <- read.csv("data/Hawk data.csv")
view(data)
inspect_df(data)

#scatterplot the relationship between jawforce and footforce 
plot(Jawforce ~ Footforce, data = data,
     xlab = "Jawforce",
     ylab = "Footforce",
     main = "Footforce by Jawforce Hawks",
     pch = 19)
abline(lm(Jawforce ~ Footforce, data = data), col = "red")
