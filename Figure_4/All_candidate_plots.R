########
# Setup 
#######

if (is.element('checkpoint', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('checkpoint'));
} else
{
  install.packages('checkpoint');
  suppressPackageStartupMessages(library('checkpoint'));
}

# Checkpoint ensures that the same version of the packages are being used
# For more info on checkpoint: https://cran.r-project.org/web/packages/checkpoint/vignettes/checkpoint.html
checkpoint("2019-10-01", checkpointLocation = tempdir())

library(tidyverse)
library(ggpubr)

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

################################################
# Load and transform individual data measurement
###############################################

#############
# volatiles #
#############

volatiles = read.delim("Figure_4/leaf_terpenoids_normalised_peak_area.tsv", header = T, 
                     stringsAsFactors = TRUE, check.names = F)


volatiles.long = gather(volatiles, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, 
                        -accession)

### Filter to keep volatiles toxic to either whitefly or thrips
candidates = read.delim(file = "Figure_4/all_candidate_names.txt",header = T,stringsAsFactors = F,check.names = F)

volatiles.long.candidates = inner_join(
                                       candidates,volatiles.long,
                                       by = "metabolite") 

### Read and add species and color information
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
volatiles.candidates.with.species = left_join(volatiles.long.candidates,accession2species,by="accession")

volatiles.candidates.with.species$accession = factor(volatiles.candidates.with.species$accession, 
                                                     levels = genotype_order_whiteflies, 
                                                     ordered = TRUE)

##############
# acylsugars #
##############

acylsugars = volatiles = read.csv("Figure_4/20190904_acylsugars_peak_area_all_samples.csv", header = T, stringsAsFactors = TRUE, check.names = F)
acylsugars.long = gather(acylsugars, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, -accession)


acylsugar.long.candidates = inner_join(candidates,acylsugars.long,by="metabolite") 


### Read and add species and color information
acylsugar.candidates.with.species = left_join(acylsugar.long.candidates,accession2species,by="accession")
acylsugar.candidates.with.species$accession = factor(acylsugar.candidates.with.species$accession, 
                                                          levels = genotype_order_whiteflies, 
                                                          ordered = TRUE)

###########################################################
# Plot 1 = barplot of selected acylsugars toxic to whitefly
###########################################################

# Theme for plotting
my.theme = theme(axis.text.x = element_text(color = "black", size = 6, angle = 45, hjust = 1),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8),
                 strip.text.x = element_text(size = 8, colour = "black"),
                 legend.text = element_text(size = 8, colour = "black")
)


# Barplot
g1 = acylsugar.candidates.with.species %>%
  filter(toxic_to == "whitefly") %>% 
  dplyr::group_by(accession, name, species, color) %>% 
  summarise(mean_abundance = mean(abundance), 
            n = n(), 
            se = (sd(abundance)/sqrt(n))
  ) %>% 
  ggplot(.) +
  geom_bar(aes(x = accession, y = mean_abundance,fill=species), stat = "identity",color="black") + 
  geom_errorbar(
    aes(x = accession, 
        ymin = mean_abundance - se, 
        ymax = mean_abundance + se)
  ) +
  facet_wrap(~ name, scale = "free", ncol = 3) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  theme_bw() +
  my.theme


###########################################################
# Plot 2 = barplot of selected volatiles toxic to whiteflies
############################################################

# Volatiles - whitefly
g2 = volatiles.candidates.with.species %>%
  filter(toxic_to == "whitefly") %>% 
  dplyr::group_by(accession, name, species, color) %>% 
  summarise(mean_abundance = mean(abundance), 
            n = n(), 
            se = (sd(abundance)/sqrt(n))
  ) %>% 
  ggplot(.) +
  geom_bar(aes(x = accession, y = mean_abundance,fill=species), stat = "identity",color="black") + 
  geom_errorbar(
    aes(x = accession, 
        ymin = mean_abundance - se, 
        ymax = mean_abundance + se)
  ) +
  facet_wrap(~ name, scale = "free", ncol = 3) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  theme_bw() +
  my.theme


# Volatiles - thrips
g3 = volatiles.candidates.with.species %>%
  filter(toxic_to == "thrips") %>% 
  dplyr::group_by(accession, name, species, color) %>% 
  summarise(mean_abundance = mean(abundance), 
            n = n(), 
            se = (sd(abundance)/sqrt(n))
  ) %>% 
  ggplot(.) +
  geom_bar(aes(x = accession, y = mean_abundance,fill=species), stat = "identity",color="black") + 
  geom_errorbar(
    aes(x = accession, 
        ymin = mean_abundance - se, 
        ymax = mean_abundance + se)
  ) +
  facet_wrap(~ name, scale = "free", ncol = 3) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  scale_x_discrete("accession", labels = genotype_order_thrips) +
  theme_bw()+
  my.theme






############
# Save plots
############
ggsave("Figure_4/Figure_4.png",plot=g,width = 8,height = 10)
ggsave("FigureS6/Figure_S6A_acylsugars_wf",plot=g1,width = 10,height = 10)
ggsave("FigureS6/Figure_S6B_volatiles_wf.pdf",plot=g2,width = 10,height = 12)
ggsave("FigureS6/Figure_S6C_volatiles_thrips.pdf",plot=g3,width = 10,height = 12)
