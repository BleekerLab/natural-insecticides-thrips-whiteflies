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

################################################
# Load and transform individual data measurement
###############################################

############
# volatiles #
############

volatiles = read.csv("Figure6/20180905_Wild_collection_leafwash.csv", header = T, stringsAsFactors = TRUE, check.names = F)

### Load data from tsv file
volatiles = read.delim("Figure6/leaf_terpenoids_normalised_peak_area.tsv", header = T, stringsAsFactors = TRUE, check.names = F)
volatiles.long = gather(volatiles, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, 
                        -accession)
volatiles.long$accession = with(volatiles.long,
                                factor(accession,levels = unique(accession)),
                                ordered= TRUE)

### Filter to keep volatiles toxic to either whitefly or thrips
candidates = read.delim("Figure6/toxic_candidate_names.txt",stringsAsFactors = F,check.names = F)
volatiles.long.candidates = inner_join(candidates,volatiles.long,by="metabolite") 

### Read and add species and color information
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
volatiles.candidates.with.species = left_join(volatiles.long.candidates,accession2species,by="accession")

##############
# acylsugars #
##############

acylsugars = volatiles = read.csv("Figure6/20190904_acylsugars_peak_area_all_samples.csv", header = T, stringsAsFactors = TRUE, check.names = F)
acylsugars.long = gather(acylsugars, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, -accession)


acylsugar.long.candidates = inner_join(candidates,acylsugars.long,by="metabolite") 


### Read and add species and color information
acylsugar.candidates.with.species = left_join(acylsugar.long.candidates,accession2species,by="accession")
acylsugar.candidates.with.species$accession = with(acylsugar.candidates.with.species,
                                                   factor(accession, 
                                                          levels = unique(accession)
                                                          ),
                                                       ordered = TRUE)


###########################################################
# Plot 1 = barplot of selected volatiles toxic to whiteflies
############################################################

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
  facet_wrap(~ name, scale = "free", ncol = 1) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  theme_bw() +
  my.theme


########################################################
# Plot 2 = barplot of selected volatiles toxic to thrips
########################################################

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
  facet_wrap(~ name, scale = "free", ncol = 1) +
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
  facet_wrap(~ name, scale = "free", ncol = 1) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  theme_bw() +
  my.theme

#############################
# Arrange the plots together
############################
ggarrange(g1,g2,g3,ncol = 3,nrow = 1,common.legend = TRUE)

g <- ggarrange(g1,g2,g3,ncol = 3,nrow = 1,common.legend = TRUE)

############
# Save plots
############
ggsave("Figure6/Figure6.png",plot=g,width = 12,height = 8)
ggsave("Figure6/Figure6.pdf",plot=g,width = 12,height = 8)



