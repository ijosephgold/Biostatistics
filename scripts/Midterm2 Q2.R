install_biostats_packages()
load_biostats_packages()
argasm_data <- read.csv("data/Sargassum data - Sargassum data.csv")
view(argasm_data)
analysis_start(argasm_data)
inspect_df(argasm_data)
biostat_inspect(argasm_data)
analysis_start(argasm_data, question = "")
view(argasm_data, wide)
biostat_wizard(argasm_data)

graphdata <- argasm_data %>%
  group_by(treatment) %>%
  summarize(meanbiomass=mean(biomass, na.rm=TRUE), se=sd(biomass, na.rm=TRUE)/sqrt(30))
graphdata
ggplot(graphdata, aes(x=treatment, y=meanbiomass)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_bar(stat="identity", color="dodgerblue2", fill="dodgerblue2", position="dodge", size=0.6) +
  labs(x="treatment", y="Biomass", size = 1) +
  geom_errorbar(aes(ymax=meanbiomass+se, ymin=meanbiomass-se), position=position_dodge(0.9), width=0.1) 

