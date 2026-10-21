#Run these two together to get regexplain powers
install.packages("remotes")
remotes::install_github("gadenbuie/regexplain")

#
install.packages("easystats")
library(easystats)

install.packages("gtsummary")
library(gtsummary)

install.packages("esquisse")
# Open via RStudio Addins menu -> "ggplot2 builder"

#Excel files integration easily 
install.packages("readxl")
library(readxl)

# Descriptive baseline table split by treatment group:
tbl_summary(data, by = treatment) |> add_p()

install.packages("rex")
library(rex)
install.packages("unglue")
library(unglue)
library(tidyverse)
library(car)
library(moments)
library(pubr)
library(ggplot2)