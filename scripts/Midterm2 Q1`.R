pen_data <- read.csv("data/African Penguin - Sheet1 (1).csv")
view(pen_data)
biostat_wizard(pen_data)

independent_t_test(pen_data, outcome = "X5", group = "Male.weight..kg.")

qqPlot(pen_data$X5, main = "QQ Plot of X5")

