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

genotype_order_for_plots = c("MM", "LA4024", "LA2133", "LA0735",
                             "LA1840", "LA1364", "LA1578",
                             "LA1278", "LA1401", "LA2172", "LA0407",
                             "LA1718", "LA1954", "PI127826",
                             "LA1777", "PI134418", "LYC4", "LA0716", "LA2695")

################################################
# Load and transform individual data measurement
###############################################

#############
# volatiles #
#############

volatiles = read.csv("Figure_S6/20180905_Wild_collection_leafwash.csv", header = T, 
                     stringsAsFactors = TRUE, check.names = F)

### Load data from tsv file
volatiles = read.delim("Figure_S6/leaf_terpenoids_normalised_peak_area.tsv", header = T, 
                       stringsAsFactors = TRUE, check.names = F)

volatiles.long = gather(volatiles, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, 
                        -accession)

### Filter to keep volatiles toxic to either whitefly or thrips
candidates = read.delim(file = "Figure_S6/all_candidate_names.txt",header = T,stringsAsFactors = F,check.names = F)

volatiles.long.candidates = inner_join(
                                       candidates,volatiles.long,
                                       by = "metabolite") 

### Read and add species and color information
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
volatiles.candidates.with.species = left_join(volatiles.long.candidates,accession2species,by="accession")

volatiles.candidates.with.species$accession = factor(volatiles.candidates.with.species$accession, 
                                                     levels = genotype_order_for_plots, 
                                                     ordered = TRUE)

##############
# acylsugars #
##############

acylsugars = volatiles = read.csv("Figure_S6/20190904_acylsugars_peak_area_all_samples.csv", header = T, stringsAsFactors = TRUE, check.names = F)
acylsugars.long = gather(acylsugars, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, -accession)


acylsugar.long.candidates = inner_join(candidates,acylsugars.long,by="metabolite") 


### Read and add species and color information
acylsugar.candidates.with.species = left_join(acylsugar.long.candidates,accession2species,by="accession")
acylsugar.candidates.with.species$accession = factor(acylsugar.candidates.with.species$accession, 
                                                          levels = genotype_order_for_plots, 
                                                          ordered = TRUE)


###########################################################
# Plot 1 = barplot of selected volatiles toxic to whiteflies
############################################################

# Theme for plotting
my.theme = theme(axis.text.x = element_text(color = "black", size = 6, angle = 45, hjust = 1),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8),
                 strip.text.x = element_text(size = 8, colour = "black"),
                 legend.text = element_text(size = 8, colour = "black")
)


# Volatiles - whitefly
g1 = volatiles.candidates.with.species %>%
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
g2 = volatiles.candidates.with.species %>%
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
  theme_bw() +
  my.theme



###########################################################
# Plot 3 = barplot of selected acylsugars toxic to whitefly
###########################################################

# Barplot
g3 = acylsugar.candidates.with.species %>%
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
  #scale_colour_manual(values=volatiles.candidates.with.species$color) +
  theme_bw() +
  my.theme

############################

############
# Save plots
############
ggsave("Figure_S6/Figure_S6A_wf_volatiles.pdf",plot=g1,width = 10,height = 10)
ggsave("Figure_S6/Figure_S6B_thrips_volatiles.pdf",plot=g2,width = 10,height = 12)
ggsave("Figure_S6/Figure_S6C_wf_acylsugars.pdf",plot=g3,width = 10,height = 10)
