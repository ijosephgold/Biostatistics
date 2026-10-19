#H_0: individual Yew conifers produce roughly equal proportions of male and female cones (ie. they're monoecious), population mean index is 0.5
#H_A: individual Yew conifers do not produce equal proportions of male and female cones (ie. they're dioecious), population mean index is not 0.5

library(tidyverse)
library(ggplot2)
library(car)
library(readxl)

yewsdata <- read.csv("data/Yew data.csv")
inspect_df(yewsdata)
yewsdata
qqPlot(yewsdata$Male, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Male Cones Q-Q plot", col.lines = "red")

yewssqrt_Male <- sqrt(yewsdata$Male)
qqPlot(yewssqrt_Male, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Male Cones SQRT Q-Q plot", col.lines = "red")

#yewslog_Male <- log(yewsdata$Male)
#qqPlot(yewslog_Male, distribution = "norm", envelope = list(level = 0.95, col = "blue"), main = "Male Cones LOG Q-Q plot", col.lines = "red") #can't graph Infinite values
       