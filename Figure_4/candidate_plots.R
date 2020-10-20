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
library(patchwork)

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

volatiles = read.csv("Figure_4/20180905_Wild_collection_leafwash.csv", header = T, 
                     stringsAsFactors = TRUE, check.names = F)

### Load data from tsv file
volatiles = read.delim("Figure_4/leaf_terpenoids_normalised_peak_area.tsv", header = T, 
                       stringsAsFactors = TRUE, check.names = F)

volatiles.long = gather(volatiles, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, 
                        -accession)

### Filter to keep volatiles toxic to either whitefly or thrips
candidates = read.delim(file = "Figure_4/toxic_candidate_names.tsv",header = T,stringsAsFactors = F,check.names = F)

volatiles.long.candidates = inner_join(volatiles.long,
                                       candidates,
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

acylsugars = read.csv("Figure_4/20190904_acylsugars_peak_area_all_samples.csv", header = T, stringsAsFactors = TRUE, check.names = F)
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
                 legend.text = element_text(size = 8, colour = "black"),
                 legend.position = "none"
)


# Barplot
g1 = acylsugar.candidates.with.species %>%
  filter(., name %in% c("S3:15", "S3:21", "S4:17-2", "S4:22")) %>%
  dplyr::group_by(accession, name, species, color) %>% 
  summarise(mean_abundance = mean(log(abundance+1)), 
            n = n(), 
            se = (sd(log(abundance+1))/sqrt(n))
  ) %>% 
  ggplot(.) +
  geom_bar(aes(x = accession, y = mean_abundance,fill=species), stat = "identity",color="black") + 
  geom_errorbar(
    aes(x = accession, 
        ymin = mean_abundance - se, 
        ymax = mean_abundance + se,
        width = 0.4)
  ) +
  facet_wrap(~ name, scale = "free", ncol = 1) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  scale_x_discrete("accession", labels = genotype_order_whiteflies)+
  theme_bw() +
  my.theme


###########################################################
# Plot 2 = barplot of selected volatiles toxic to whiteflies
############################################################


 g2 = volatiles.candidates.with.species %>%
  dplyr::group_by(accession, name, species, color) %>% 
  summarise(mean_abundance = mean(log(abundance+1)), 
            n = n(), 
            se = (sd(log(abundance+1))/sqrt(n))
  ) %>% 
  filter(., name %in% c("10-epi-italicene ether", "y-cuprenene")) %>%
  ggplot(.) +
  geom_bar(aes(x = accession, y = mean_abundance,fill=species), stat = "identity",color="black") + 
  geom_errorbar(
    aes(x = accession, 
        ymin = mean_abundance - se, 
        ymax = mean_abundance + se,
        width = 0.4)
  ) +
  facet_wrap(~ name, scale = "free", ncol = 1) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  scale_x_discrete("accession", labels = genotype_order_whiteflies)+
  theme_bw() +
  my.theme+
  theme(axis.title.x = element_blank())
 
 ###########################################################
 # Plot 3 = barplot of selected volatiles toxic to Thrips
 ############################################################

g3 = volatiles.candidates.with.species %>%
  dplyr::group_by(accession, name, species, color) %>% 
  summarise(mean_abundance = mean(log(abundance+1)), 
            n = n(), 
            se = (sd(log(abundance+1))/sqrt(n))
  ) %>% 
  filter(., name %in% c("d-selinene", "b-cadinene")) %>%
  ggplot(.) +
  geom_bar(aes(x = accession, y = mean_abundance,fill=species), stat = "identity",color="black") + 
  geom_errorbar(
    aes(x = accession, 
        ymin = mean_abundance - se, 
        ymax = mean_abundance + se,
        width = 0.4)
  ) +
  facet_wrap(~ name, scale = "free", ncol = 1) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  scale_x_discrete("accession", labels = genotype_order_thrips)+
  theme_bw() +
  my.theme




#############################
# Arrange the plots together
############################

plot.volatiles = ggarrange(g2,g3, ncol = 1,nrow = 2)

############
# Save plots
############

ggsave("Figure_4/Figure_4_acylsugars.pdf", plot = g1, width = 3.5,height = 7)
ggsave("Figure_4/Figure_4_volatiles.pdf", plot = plot.volatiles, width = 3.5, height = 7)



