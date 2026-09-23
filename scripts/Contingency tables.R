#Contingency tables

#load libraries 
library(tidyverse)
library(DescTools) #for G test
library(ggplot2) #for plotting

#analyze dispersal and reproduction data

#create data
antseeds <- 25
antveg <- 36
vertseeds <- 6
vertveg <- 21

#calculate totals
anttotal <- sum(antseeds, antveg)
verttotal <- sum(vertseeds, vertveg)
seedtotal <- sum(antseeds, vertseeds)
vegtotal <- sum(antveg, vertveg)
rowcoltotal <- anttotal+verttotal

#make a contingency table of observed values
contable_obs <- matrix(c(antseeds, antveg, anttotal, vertseeds, vertveg, verttotal,
                   seedtotal, vegtotal, rowcoltotal), nrow=3, ncol=3)
contable_obs
#name rows and columns
row.names(contable_obs) <- (c("Seed only", "Vegetative", "Total"))
colnames(contable_obs) <- (c("Ants", "Vertebrates", "Total"))
contable_obs

#calculate expected values
antseeds_exp <- (31/88)*(61/88)*88
antveg_exp <- (57/88)*(61/88)*88
vertseeds_exp <- (31/88)*(27/88)*88
vertveg_exp <- (57/88)*(27/88)*88

#make a contingency table of expected values
contable_exp <- matrix(c(antseeds_exp, antveg_exp, anttotal, vertseeds_exp, vertveg_exp, 
                       verttotal, seedtotal, vegtotal, rowcoltotal), nrow=3, ncol=3)
contable_exp
#name rows and columns
row.names(contable_exp) <- (c("Seed only", "Vegetative", "Total"))
colnames(contable_exp) <- (c("Ants", "Vertebrates", "Total"))
contable_exp

#test of independence
#reassemble contingency table of observed values without totals
contable <- matrix(c(antseeds, antveg, vertseeds, vertveg), nrow=2, ncol=2)
row.names(contable) <- (c("Seed only", "Vegetative"))
colnames(contable) <- (c("Ants", "Vertebrates"))
#perform G test
GTest(contable)

#make plot
names <- matrix(c("Seed", "Veg"))
barplot(contable, beside = TRUE, names.arg = names(contable), legend.text = (names), 
        bty = "n", ylim = c(0, 40), ylab = "Counts", args.legend = list(x = "topright", 
        inset = c(-0.1, -0.2),bty = "n", cex = 1))

#analyze bird foraging guild counts before and after wildfire

#load data
#first find out the working directory of the R project
getwd()
#copy this into read.csv and add /Data/data.csv to the end of the address
Bird_guild_data <- safe_read_csv("data/Bird guild data.csv")
#or do this 
guilds <- read.csv("Data/Bird guild data.csv")

#sum count data for each guild category - this is tedious, but do it anyway
#the %>% operator is called a pipe - it tells R to take the data and draw pieces from it
carnivore <- guilds %>%
  filter(Guild_type == "carnivore") 
precarn <- sum(carnivore$Prefire_counts)
postcarn <- sum(carnivore$Postfire_counts)

nectar <- guilds %>%
  filter(Guild_type == "nectar") 
prenec <- sum(nectar$Prefire_counts)
postnec <- sum(nectar$Postfire_counts)

granivore <- guilds %>%
  filter(Guild_type == "granivore") 
pregran <- sum(granivore$Prefire_counts)
postgran <- sum(granivore$Postfire_counts)

insectivore <- guilds %>%
  filter(Guild_type == "insectivore") 
preins <- sum(insectivore$Prefire_counts)
postins <- sum(insectivore$Postfire_counts)

omnivore <- guilds %>%
  filter(Guild_type == "omnivore") 
preomn <- sum(omnivore$Prefire_counts)
postomn <- sum(omnivore$Postfire_counts)

#create matrix
guildmatrix <- matrix(c(precarn, prenec, pregran, preins, preomn,
                        postcarn, postnec, postgran, postins, postomn), 
                      nrow=5, ncol=2)

#name rows and columns
rownames(guildmatrix) <- c("Carnivore", "Nectar", "Granivore", "Insectivore", "Omnivore")
colnames(guildmatrix) <- c("Pre-fire", "Post-fire")

#check data frame
guildmatrix 

#perform G-test
GTest(guildmatrix)

#make a plot for better visualization
#first need new data frame in long format
#save some work by putting it together in Excel first
gtestdata <- read.csv("Data/Gtest data.csv")
gtestdata

#make a plot with ggplot
ggplot(gtestdata, aes(x = Period, y = Count, fill = Guild, group = Guild))+
  geom_bar(aes(x = factor(Period, level = c ('Pre-fire', 'Post-fire'))),
    stat="identity", position="dodge", size=0.6)+
  guides(fill=guide_legend(title=NULL))+
  xlab("") + ylab("Bird counts")+
  ylim(0, 200)+
  theme_bw()+
  theme(panel.grid.major=element_blank())+
  theme(axis.text.x = element_text(size = 11))

