edata <- read.csv("data/Environmental data.csv")
view(edata)
#{scatterplotMatrix(edata)
#scatterplot(edata$Elevation, edata$MaxTemp)
#pearsons_c(edata$Elevation, edata$MaxTemp)
#pearsons_c(edata$PctCover, edata$Elevation)
#pearsons_c(edata$PctCover, edata$TPI)
#pearsons_c(edata$PctCover, edata$TWI)
#pear}

vari <- c("PctCover", "Elevation", "TPI", "TWI", "MaxTemp")
dat <- edata[complete.cases(edata[vari]), vari] #complete.cases(...) checks each row, returns TRUE if it has a value for every selected column, FALSE if anyis missing

scatterplotMatrix(dat)

pairs <- combn(vari, 2, simplify = FALSE)

results <- do.call(rbind, lapply(pairs, function(x) {
  test <- cor.test(dat[[x[1]]], dat[[x[2]]], method = "pearson")
  data.frame(
    variable1 = x[1],
    variable2 = x[2],
    r = unname(test$estimate),
    p = test$p.value
  )
}))

results$p_adjusted <- p.adjust(results$p, method = "holm")
results