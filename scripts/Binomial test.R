#Binomial test

#sign test for ecological footprints with 9 pluses (k=9)
#create a binomial probability distribution of n = 20 trials
dist <- dbinom(x = 0:20, size = 20, p = 0.5)
dist

#plot the distribution
barplot(dist, ylim= c(0,0.20), names.arg = 0:20, xlab = 'x', ylab = 'p')
#add color to the values of the distribution 0:9 
barplot(dist, ylim= c(0,0.20), names.arg = 0:20, xlab = 'x', ylab = 'p',
        col = c(rep("blue", 9), rep("gray", 11)))

#summarize the first 9 values and multiply by 2 (for a two-tailed test)
#gives the probability of observing k = 9
p <- 2*dbinom(x=0:9, size = 20, p = 0.5)
sum(p) #p-value

#check answer
binom.test(9, 20, 0.5)

#sign test for birds above 1500 m (1 plus, k=1)
#create a binomial probability distribution of n = 28 trials
dist <- dbinom(x = 0:28, size = 28, p = 0.5)
dist

#plot the distribution
barplot(dist, ylim= c(0,0.20), names.arg = 0:28, xlab = 'x', ylab = 'p')
#add color to the values of the distribution 0:1
barplot(dist, ylim= c(0,0.20), names.arg = 0:28, xlab = 'x', ylab = 'p', 
        col = c(rep("blue", 1), rep("gray", 27)))

#summarize the first 9 values and multiply by 2 
#gives the probability of observing k = 1
p <- 2*dbinom(x=0:1, size = 28, p = 0.5)
sum(p) #p-value

#check answer
binom.test(1, 28, 0.5)
