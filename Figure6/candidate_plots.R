library(tidyverse)

volatiles = read.csv("Figure5/Figure_5_candidates_plots/20180905_Wild_collection_leafwash.csv", header = T, stringsAsFactors = TRUE, check.names = F)

volatiles.long = gather(volatiles, 
                        key = "metabolite",
                        value = "abundance",
                        -Sample, -Accession)

volatiles.long$Accession = factor(volatiles.long$Accession, levels = c("MM", "LA4024", "LA2133", "LA0735", "LA1840", "LA1364", "LA1578",
                                                                       "LA1278", "LA1401", "LA2172", "LA0407",
                                                                       "LA1718", "LA1954", "PI127826",
                                                                       "LA1777", "PI134418", "LYC4", "LA0716", "LA2695"), 
                                  ordered = TRUE)

volatiles.candidates = volatiles.long %>% filter(., metabolite %in% c('11.844_91.0573', '25.356_105.0726', '26.164_119.0865', '25.421_161.1340','25.968_119.0881'))



#########################################
# Create barplots of selected volatiles #
#########################################

# Theme for plotting
my.theme = theme(axis.text.x = element_text(color = "black", size = 6, angle = 45, hjust = 1),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8)
)


#Boxplot
ggplot(volatiles.candidates)+
  geom_boxplot(aes(x = Accession, y = abundance))+
  facet_wrap(~metabolite, scale = "free" , ncol = 1)+
  theme_bw()+
  my.theme

#Barplot
volatiles.candates.avg = volatiles.candidates %>% dplyr::group_by(Accession, metabolite) %>% 
  summarise(., mean_abundance = mean(abundance), n = n(), se = (sd(abundance)/sqrt(n)))


ggplot(volatiles.candates.avg)+
geom_bar(aes(x = Accession, y = mean_abundance), stat = "identity")+
  geom_errorbar(aes(x = Accession, ymin = mean_abundance - se, ymax = mean_abundance + se))+
  facet_wrap(~metabolite, scale = "free", ncol = 1)+
  theme_bw()+
  my.theme
  





