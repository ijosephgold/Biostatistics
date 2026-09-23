logcirc <- log(data$Circumference..cm.)
hist(logcirc, breaks = 30, main = "", xlab= "Tree height cm", xlim = c(0,10), ylim = c(0,25)) 
abline(v=lowerCL, col="blue")
abline(v=upperCL, col="blue")