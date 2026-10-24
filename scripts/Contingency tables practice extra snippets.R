

# View observed contingency table
contable_obs

# Add across each ROW
rowSums(contable_obs)

# Add down each COLUMN
colSums(contable_obs)

# Add EVERY observation = grand total
sum(contable_obs)

# Display contingency table WITH row and column totals
addmargins(contable_obs)

#calculate totals
redtotal <- sum(redrural, redindustrial)
blacktotal <- sum(blackrural, blackindustrial)
ruraltotal <- sum(redrural, blackrural)
industrialtotal <- sum(redindustrial, blackindustrial)
rowcoltotal <- redtotal+blacktotal

#make a contingency table of observed values
contable_obs <- matrix(c(redrural, redindustrial, redtotal, blackrural, blackindustrial, blacktotal,
                         ruraltotal, industrialtotal, rowcoltotal), nrow=3, ncol=3)
contable_obs
#name rows and columns
row.names(contable_obs) <- (c("Rural", "Industrial", "Total"))
colnames(contable_obs) <- (c("Red", "Black", "Total"))
contable_obs

# Add across each ROW
rowSums(contable_obs)

# Add down each COLUMN
colSums(contable_obs)

# Add every observation in the table
sum(contable_obs)

addmargins(contable_obs)

#calculate expected values
redrural_exp <- (100/300)*(155/300)*300
redindustrial_exp <- (145/300)*(155/300)*300
blackrural_exp <- (100/300)*(145/300)*300
blackindustrial_exp <- (145/300)*(145/300)*300

#make a contingency table of expected values
contable_exp <- matrix(c(redrural_exp, redindustrial_exp, redtotal, blackrural_exp, blackindustrial_exp, 
                         blacktotal, ruraltotal, industrialtotal, rowcoltotal), nrow=3, ncol=3)
contable_exp
#name rows and columns
row.names(contable_exp) <- (c("Rural", "Industrial", "Total"))
colnames(contable_exp) <- (c("Red", "Black", "Total"))
contable_exp

#test of independence
#reassemble contingency table of observed values without totals
contable <- matrix(c(redrural, redindustrial, blackrural, blackindustrial), nrow=2, ncol=2)
row.names(contable) <- (c("Rural", "Industrial"))
colnames(contable) <- (c("Red", "Black"))
#perform G test
GTest(contable)

#make plot
names <- matrix(c("Rural", "Ind"))
barplot(contable, beside = TRUE, names.arg = names(contable), legend.text = (names), 
        bty = "n", ylim = c(0, 40), ylab = "Counts", args.legend = list(x = "topright", 
                                                                        inset = c(-0.1, -0.2),bty = "n", cex = 1))