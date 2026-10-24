#Create Data
ConeFlying <- 
  BarkFlying <-
  LeavesFlying <-
  BranchFlying  <-
  
  flyingcone <- data %>%
  filter(Behavior == "Flying")

Countdata <- function(data, x, y){
  countdata <- data %>%
    filter(behavior == x) %>%
    filter(Substrate == y) %>%
    summarise(Count == sum(Count))
  unlist(Countdata)
}

make_contingency_table <- function(data) {
  xtabs(Count ~ Substrate + Behavior, data = data)
}

contable_obs <- make_contingency_table(data)
contable_obs

contable_obs <- matrix(c())

library(ggplot2)

ggplot(data, aes(x = Substrate, y = Count, fill = Behavior)) + 
  geom_col(position = "dodge") + 
  labs(
    title = "Foraging Behaviors of Mountain Chickadees by Substrate",
    x = "Substrate",
    y = "Observed Count",
    fill = "Behavior"
  ) + 
  theme_minimal()





fisher.test(contable_obs, simulate.p.value=TRUE)

library(ggplot2)


ConeFlying
2	Bark	Flying
1	Leaves	Flying
1	Branch	Flying
0	Cone	Hopping
0	Bark	Hopping
13	Leaves	Hopping
23	Branch	Hopping
0	Cone	Climbing
1	Bark	Climbing
0	Leaves	Climbing
0	Branch	Climbing
3	Cone	Gleaning
0	Bark	Gleaning
10	Leaves	Gleaning
5	Branch	Gleaning
2	Cone	Probing
0	Bark	Probing
0	Leaves	Probing
10	Branch	Probing
ConeHanging
BarkHanging
LeavesHanging
BranchHanging