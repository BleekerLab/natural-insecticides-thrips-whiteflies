library(tidyverse)
library(Rmisc)
library(pals)
library(reshape2)
library(ggrepel)

#############################
# Custom theme for plotting #
#############################

my.theme = theme(axis.text.x = element_text(color = "black", size = 6),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8)
)


###############
# Import data #
###############
volatiles <- read_delim(file = "Figure_3/volatile_analytics/GCMS_normalised_counts.txt", delim = "\t") %>%
  mutate(sum_volatiles = rowSums(volatiles[3:ncol(volatiles)])
  ) %>%
  select(accession, sum_volatiles)

volatiles[is.na(volatiles)] <- 0

accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
volatiles.with.species = left_join(volatiles,accession2species,by="accession") %>% select(-genotype)


# Load insect survival data
pheno = read.delim("Figure_1/dataframe_for_relative_scatterplots.tsv", header = T) %>% 
  select(., c("accession", "species", "color", "wf_relative_survival", "thrips_relative_survival")) %>%
  select(-c(species, color))

volatiles.pheno <- left_join(volatiles.with.species, pheno, by = "accession")


#########
# plots #
#########

p.wf = 
volatiles.pheno %>%
  dplyr::group_by(accession, species, color, wf_relative_survival, thrips_relative_survival) %>%
  dplyr::summarise(mean_abundance = mean(sum_volatiles)) %>%
  ggplot() +
  geom_point(aes(x = wf_relative_survival, y = mean_abundance),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = wf_relative_survival, y=mean_abundance,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  geom_smooth(aes(x=wf_relative_survival, y = mean_abundance), method = "lm", alpha = 0.2, color = "black", size = 0.2)+
  labs(x = "Relative whitefly survival (%)",y = "Total volatiles (ion counts / mg FW)") +
  my.theme+
  ylim(0, NA)

p.thrips = 
volatiles.pheno %>%
  dplyr::group_by(accession, species, color, thrips_relative_survival, thrips_relative_survival) %>%
  dplyr::summarise(mean_abundance = mean(sum_volatiles)) %>%
  ggplot() +
  geom_point(aes(x = thrips_relative_survival, y = mean_abundance),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = thrips_relative_survival, y=mean_abundance,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  geom_smooth(aes(x=thrips_relative_survival, y = mean_abundance), method = "lm", alpha = 0.2, color = "black", size = 0.2)+
  labs(x = "Relative thrips survival (%)",y = "Total volatiles (ion counts / mg FW)") +
  my.theme+
  ylim(0, NA)

#################
# Linear models #
#################

lm.wf.volatiles = summary( 
  lm(data = volatiles.pheno %>%
       dplyr::group_by(accession, species, color, wf_relative_survival, thrips_relative_survival) %>%
       dplyr::summarise(mean_abundance = mean(sum_volatiles)),
     wf_relative_survival~ mean_abundance)
)

lm.thrips.volatiles = summary( 
  lm(data = volatiles.pheno %>%
       dplyr::group_by(accession, species, color, wf_relative_survival, thrips_relative_survival) %>%
       dplyr::summarise(mean_abundance = mean(sum_volatiles)),
     thrips_relative_survival~ mean_abundance)
)

########
# Save #
########

ggsave(filename = file.path("Figure_3/volatile_analytics/wf_vs_total_volatiles.pdf"),plot = p.wf, width = 8,height = 6, units = "cm")
ggsave(filename = file.path("Figure_3/volatile_analytics/thrips_vs_total_volatiles.pdf"),plot = p.thrips, width = 8,height = 6, units = "cm")

capture.output(lm.wf.volatiles, file = "Figure_3/volatile_analytics/linear_model_wf_volatiles")
capture.output(lm.thrips.volatiles, file = "Figure_3/volatile_analytics/linear_model_thrips_volatiles")

