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

################################################
# Load and transform individual data measurement
###############################################
volatiles = read.csv("Figure6/20180905_Wild_collection_leafwash.csv", header = T, stringsAsFactors = TRUE, check.names = F)

volatiles.long = gather(volatiles, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, -accession)

volatiles.long$Accession = factor(volatiles.long$accession, levels = c("MM", "LA4024", "LA2133", "LA0735", "LA1840", "LA1364", "LA1578",
                                                                       "LA1278", "LA1401", "LA2172", "LA0407",
                                                                       "LA1718", "LA1954", "PI127826",
                                                                       "LA1777", "PI134418", "LYC4", "LA0716", "LA2695"), 
                                  ordered = TRUE)

volatiles.candidates = volatiles.long %>% filter(., metabolite %in% c('11.844_91.0573', '25.356_105.0726', '26.164_119.0865', '25.421_161.1340','25.968_119.0881'))


# Read and add species and color information
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
volatiles.candidates.with.species = left_join(volatiles.candidates,accession2species,by="accession")
volatiles.candidates.with.species$accession = factor(volatiles.candidates.with.species$accession, levels = c("MM", "LA4024", "LA2133", "LA0735", "LA1840", "LA1364", "LA1578",
                                                                       "LA1278", "LA1401", "LA2172", "LA0407",
                                                                       "LA1718", "LA1954", "PI127826",
                                                                       "LA1777", "PI134418", "LYC4", "LA0716", "LA2695"), 
                                  ordered = TRUE)

#########################################
# Create barplots of selected volatiles #
#########################################

# Theme for plotting
my.theme = theme(axis.text.x = element_text(color = "black", size = 6, angle = 45, hjust = 1),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8)
)


# Boxplot
g1 <- ggplot(volatiles.candidates.with.species) +
  geom_boxplot(aes(x = accession, y = abundance,fill=species)) +
  facet_wrap(~ metabolite, scale = "free" , ncol = 1) +
  theme_bw() +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  my.theme
g1

# Barplot
volatiles.candidates.avg = volatiles.candidates.with.species %>%
  dplyr::group_by(accession, metabolite, species, color) %>% 
  summarise(mean_abundance = mean(abundance), 
            n = n(), 
            se = (sd(abundance)/sqrt(n))
            )

g2 = ggplot(volatiles.candidates.avg) +
  geom_bar(aes(x = accession, y = mean_abundance,fill=species), stat = "identity",color="black") + 
  geom_errorbar(
    aes(x = accession, 
        ymin = mean_abundance - se, 
        ymax = mean_abundance + se)
    ) +
  facet_wrap(~ metabolite, scale = "free", ncol = 1) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_colour_manual(values=volatiles.candidates.with.species$color) +
  theme_bw() +
  my.theme
g2

############
# Save plots
############
ggsave("Figure6/Figure6_boxplot.png",plot=g1,width = 12,height = 8)
ggsave("Figure6/Figure6_boxplot.pdf",plot=g1,width = 12,height = 8)


ggsave("Figure6/Figure6_barplot.png",plot=g2,width = 12,height = 8)
ggsave("Figure6/Figure6_barplot.pdf",plot=g2,width = 12,height = 8)
  




