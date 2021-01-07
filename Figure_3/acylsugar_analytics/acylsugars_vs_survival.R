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

#######################
# Load acylsugar data #
#######################


acylsugars = read.csv("Figure_S6/20190904_acylsugars_peak_area_all_samples.csv", header = T, stringsAsFactors = TRUE, check.names = F) %>%
  mutate(sum_acylsugars = rowSums(acylsugars[3:ncol(acylsugars)]))

acylsugars.long = gather(acylsugars, 
                         key = "metabolite",
                         value = "abundance",
                         -sample, -accession) %>%
  filter(metabolite == "sum_acylsugars")

# Load species information
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
acylsugars.with.species = left_join(acylsugars.long,accession2species,by="accession") %>% select(-genotype)

# Load insect phenotypes
pheno = read.delim("Figure_1/dataframe_for_relative_scatterplots.tsv", header = T) %>% 
  select(., c("accession", "species", "color", "wf_relative_survival", "thrips_relative_survival")) %>%
  select(-c(species, color))

acylsugars.pheno <- left_join(acylsugars.with.species, pheno, by = "accession")


# Group according to acylsugar abundance, summarise and plot data

p.wf = 
acylsugars.pheno %>%
  dplyr::group_by(accession, species, color, wf_relative_survival, thrips_relative_survival) %>%
  dplyr::summarise(mean_abundance = mean(abundance)) %>%
  ggplot() +
  geom_point(aes(x = wf_relative_survival, y = mean_abundance),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = wf_relative_survival, y=mean_abundance,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  geom_smooth(aes(x=wf_relative_survival, y = mean_abundance), method = "lm", alpha = 0.2, color = "black", size = 0.2)+
  labs(x = "Relative whitefly survival (%)",y = "Total acylsugars (ion counts / mg FW)") +
  my.theme+
  ylim(0, NA)

# thrips

p.thrips =
acylsugars.pheno %>%
  dplyr::group_by(accession, species, color, wf_relative_survival, thrips_relative_survival) %>%
  dplyr::summarise(mean_abundance = mean(abundance)) %>%
  ggplot() +
  geom_point(aes(x = thrips_relative_survival, y = mean_abundance),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = thrips_relative_survival, y=mean_abundance,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  geom_smooth(aes(x=thrips_relative_survival, y = mean_abundance), method = "lm", alpha = 0.2, color = "black", size = 0.2)+
  labs(x = "Relative thrips survival (%)",y = "Total acylsugars (ion counts / mg FW)") +
  my.theme+
  ylim(0, NA)

#################
# Linear models #
#################

lm.wf.acylsugars = summary( 
  lm(data = acylsugars.pheno %>%
       dplyr::group_by(accession, species, color, wf_relative_survival, thrips_relative_survival) %>%
       dplyr::summarise(mean_abundance = mean(abundance)),
     wf_relative_survival~ mean_abundance)
)

lm.thrips.acylsugars = summary( 
  lm(data = acylsugars.pheno %>%
       dplyr::group_by(accession, species, color, wf_relative_survival, thrips_relative_survival) %>%
       dplyr::summarise(mean_abundance = mean(abundance)),
     thrips_relative_survival~ mean_abundance)
)

########
# Save #
########

ggsave(filename = file.path("Figure_3/acylsugar_analytics/wf_vs_total_acylsugars.pdf"),plot = p.wf, width = 8,height = 6, units = "cm")
ggsave(filename = file.path("Figure_3/acylsugar_analytics/thrips_vs_total_acylsugars.pdf"),plot = p.thrips, width = 8,height = 6, units = "cm")

capture.output(lm.wf.acylsugars, file = "Figure_3/acylsugar_analytics/linear_model_wf_acylsugars")
capture.output(lm.thrips.acylsugars, file = "Figure_3/acylsugar_analytics/linear_model_thrips_acylsugars")
