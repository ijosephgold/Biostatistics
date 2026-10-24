#H_0: The Abundance of Pacific Purple Sea Urchins is the same in both Rocky & Flat Substrates
#H_A: The Abundance of Pacific Purple Sea Urchins Differ Between Rocky & Flat Substrates
  
Urchins <- read.csv("data/PSet4 Q3 - Sheet1.csv")
view(Urchins)
inspect_df(Urchins)
VarianceFlat  <- var(Urchins$Flat)
VarianceRocky <- var(Urchins$Rocky)

CoeficientofVariationFlat <- sqrt(VarianceFlat) / mean(Urchins$Flat)
CoeficientofVariationRocky <- sqrt(VarianceRocky) / mean(Urchins$Rocky)

Ttest <- t.test(Urchins$Flat, Urchins$Rocky, alternative = "two.sided")
Ttest

#Conclusion: Our T-test statistic returned a P-value of 0.05622 which is just above the threshhold
#of 0.05, suggesting marginally significant difference in Pacific Purple Sea Urchin between 
#the Rocky and Flat Substrates