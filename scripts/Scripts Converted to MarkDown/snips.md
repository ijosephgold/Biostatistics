```r
krats <- read.csv("data/Krat data.csv")
view(krats)

vari <- c("RatDnsty", "ShrbCvr", "SeedPrd", "SnkDnsty")
dat <- krats[complete.cases(krats[vari]), vari] #complete.cases(...) checks each row, returns TRUE if it has a value for every selected column, FALSE if anyis missing
scatterplotMatrix(dat)

pairs <- combn(vari, 2, simplify = FALSE)

```
