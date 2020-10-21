library(tidyverse)
library(Rmisc)
library(pals)
library(reshape2)

#######################
# Load acylsugar data #
#######################


acylsugars = volatiles = read.csv("Figure_S6/20190904_acylsugars_peak_area_all_samples.csv", header = T, stringsAsFactors = TRUE, check.names = F)
acylsugars.long = gather(acylsugars, 
                         key = "metabolite",
                         value = "abundance",
                         -sample, -accession)

#######################
# Add metabolite info #
#######################

acylsugars.long.structure <- acylsugars.long %>%
  filter(!grepl("unknown", metabolite, ignore.case = TRUE)) %>%
  mutate(backbone = case_when(grepl("S", metabolite) ~ "Sucrose",
                              grepl("G", metabolite) ~ "Glucose")) %>%
  mutate(tails = case_when(grepl("G1", metabolite) ~ "1",
                           grepl("G2", metabolite) ~ "2",
                          grepl("G3", metabolite) ~ "3",
                          grepl("G4", metabolite) ~ "4",
                          grepl("G5", metabolite) ~ "5",
                          grepl("G4", metabolite) ~ "4",
                          grepl("S1", metabolite) ~ "1",
                          grepl("S2", metabolite) ~ "2",
                          grepl("S3", metabolite) ~ "3",
                          grepl("S4", metabolite) ~ "4",
                          grepl("S5", metabolite) ~ "5",
                          grepl("S6", metabolite) ~ "6")
  ) %>%
  mutate(carbons = case_when(grepl("C4", metabolite) ~ "C4",
                             grepl("C5", metabolite) ~ "C10",
                             grepl("C6", metabolite) ~ "C10",
                             grepl("C7", metabolite) ~ "C10",
                             grepl("C8", metabolite) ~ "C10",
                             grepl("C9", metabolite) ~ "C9",
                           grepl("C10", metabolite) ~ "C10",
                           grepl("C11", metabolite) ~ "C11",
                           grepl("C12", metabolite) ~ "C12",
                           grepl("C13", metabolite) ~ "C13",
                           grepl("C14", metabolite) ~ "C14",
                           grepl("C15", metabolite) ~ "C15",
                           grepl("C16", metabolite) ~ "C16",
                           grepl("C17", metabolite) ~ "C17",
                           grepl("C18", metabolite) ~ "C18",
                           grepl("C19", metabolite) ~ "C19",
                           grepl("C20", metabolite) ~ "C20",
                           grepl("C21", metabolite) ~ "C21",
                           grepl("C22", metabolite) ~ "C22",
                           grepl("C23", metabolite) ~ "C23",
                           grepl("C24", metabolite) ~ "C24",
                           grepl("C25", metabolite) ~ "C25",
                           grepl("C26", metabolite) ~ "C26",
                           grepl("C27", metabolite) ~ "C27",
                           grepl("C32", metabolite) ~ "C32")
  )
  
# Add species info
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
acylsugars.with.species = left_join(acylsugars.long.structure,accession2species,by="accession")

#######################
# Ordering of the plots
#######################

genotype_order_whiteflies = c("LA0716","PI127826","LYC4", "LA1777", 
                              "PI134418", "LA1718", "LA1954","LA2695",
                              "LA1401","LA0407","LA1364","LA4024", "LA2172", 
                              "LA2133","LA1578","LA0735", "LA1278","MM","LA1840")

genotype_order_thrips = c("LYC4","LA0407", "LA1777", "PI134418",
                          "LA1401", "LA0716", "LA1278",  "LA2172", 
                          "LA2695","LA0735","LA1718","LA2133","PI127826",
                          "LA1578", "LA1954", "MM",  "LA1840", "LA4024","LA1364")

########
# plots #
########

# Set to accessions to WF survival order
acylsugars.with.species$accession = factor(acylsugars.with.species$accession,
                                              levels = genotype_order_whiteflies,
                                              ordered = TRUE)

#############################
# Plot acyl-chain proportions #
#############################

# custom colours 
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#0072B2", "#D55E00", "#CC79A7")

p.acyl.proportion = 
acylsugars.with.species %>%
  dplyr::group_by(accession, backbone, tails, species, color) %>% 
  dplyr::summarise(mean_abundance = mean(abundance), 
            n = n(), 
            se = (sd(abundance)/sqrt(n))
            ) %>%
  ungroup %>%
  ggplot(aes(x = accession, y = mean_abundance, fill = tails))+
  geom_bar(position = "fill", stat = "identity")+
  facet_wrap(~backbone, ncol = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1))+
  scale_fill_manual(values = cbp1)+
  labs(fill = "No of acyl chains")+
  ylab("Proportion of the acylsugar")

# Plot backbones

summary = summarySE(acylsugars.with.species, measurevar = "abundance", groupvars = c("accession","backbone"))
summary = within(summary, cum_abundance <- ave(abundance, accession, FUN = cumsum))

p.backbones = 
  summary %>%
  ggplot( aes(x = accession, y = abundance, fill = backbone))+
  geom_bar(position = "stack", stat = "identity")+
  geom_errorbar(aes(ymin = cum_abundance -se, 
                    ymax = cum_abundance+se,
                    width = 0.5))+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1))+
  scale_fill_manual(values = c("black","grey"))+
  labs(fill = "Sugar moiety")+
  ylab("Summed peak area (ion counts)")


# Plot number of carbons
p.carbon =
acylsugars.with.species %>%
  dplyr::group_by(accession, backbone, carbons, species, color) %>% 
  dplyr::summarise(mean_abundance = mean(abundance), 
            n = n(), 
            se = (sd(abundance)/sqrt(n))
  ) %>%
  ggplot(aes(x = accession, y = mean_abundance, fill = carbons))+
  geom_bar(position = "stack", stat = "identity")+
  facet_wrap(~backbone, ncol = 1)+
  theme_bw()+
  theme(legend.position  = "bottom",
        axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1))+
  scale_fill_manual(values=stepped())

# Number of carbons proportion

p.carbon.proportion =
  acylsugars.with.species %>%
  dplyr::group_by(accession, backbone, carbons, species, color) %>% 
  dplyr::summarise(mean_abundance = mean(abundance), 
                   n = n(), 
                   se = (sd(abundance)/sqrt(n))
  ) %>%
  ggplot(aes(x = accession, y = mean_abundance, fill = carbons))+
  geom_bar(position = "fill", stat = "identity")+
  facet_wrap(~backbone, ncol = 1)+
  theme_bw()+
  theme(legend.position  = "bottom",
        axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1))+
  scale_fill_manual(values=stepped())

##############
# Save plots #
##############

ggsave(filename = "Figure_3/acylsugar_analytics/Sugar_moiety.pdf", plot = p.backbones, height = 3.5, width = 6)
ggsave(filename = "Figure_3/acylsugar_analytics/Acyl_chains.pdf", plot = p.acyl.proportion, height = 3.5, width = 6)
ggsave(filename = "Figure_3/acylsugar_analytics/Carbons.pdf", plot = p.carbon, height = 5.5, width = 6)
ggsave(filename = "Figure_3/acylsugar_analytics/Carbons_proportion.pdf", plot = p.carbon.proportion, height = 5.5, width = 6)
