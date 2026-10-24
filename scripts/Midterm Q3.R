#Count data, so expected frequencies are fitted to a Poisson distribution 

#H_0: The data follows a Poisson distribution.
#Test of choice: Chi-Square Goodness of Fit Test
observed <- c(13, 26, 21, 13, 9)

expected <- c(12.69, 23.68, 22.09, 13.74, 9.79)

chisq.test(
  x = observed,
  p = expected / sum(expected),
  rescale.p = TRUE
)

#Poisson degrees of freedom = number of categories - 2, so 5 - 2 = 3
#because p-value is greater than 0.05 (~0.9831), we fail to reject the null hypothesis and conclude that the data follows a Poisson distribution.
