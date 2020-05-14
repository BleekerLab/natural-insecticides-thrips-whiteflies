library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(Rmisc)
library(survival)
library(survminer)

#############################
# Custom theme for plotting #
#############################

my.theme = theme(axis.text.x = element_text(color = "black", size = 6),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8)
)

#################
# Load datasets #
#################

# trichome density data
density = read.delim("FigureS4/trichome.counts.processed.tsv", header =T)

# Calculate the average trichome densities per genotype (abaxial and adaxial are taken together)
density.avg = summarySE(
  density, 
  measurevar = "density", 
  groupvars = c("genotype", "accession", "species" ,"type", "color")
)


# whitefly and thrips ranks based on toxicity
pheno = read.delim("Figure2/data4scatterplot_ranks.tsv", header = T)
pheno = separate(pheno, col = sample, into = c("x", "y", "accession")) %>% #separate "sample" into accession names only so it can be used for fusing the datasets
  select(., c("accession", "wf", "thrips"))

# Thrips median survival (from Fig 1B)
survData = read.delim("Figure1/thrips_survival_data.tsv",header=T,stringsAsFactors = F)
accession2species = read.delim("genotype2species.txt",header = T,sep = "\t",stringsAsFactors = T)
fit <- with(survData,survfit(formula = Surv(time,status) ~ accession))
df.medians = surv_median(fit)
df.medians = mutate(df.medians,strata = gsub("accession=",replacement = "",strata))
df.medians = df.medians[order(df.medians$median,decreasing = F),]
colnames(df.medians)[1]="accession"
df.medians = dplyr::left_join(df.medians,accession2species,by="accession")
df.medians$accession = factor(df.medians$accession,levels = df.medians$accession)

#whitefly phenotype (From Fig 1A)
wf = pheno = read.delim("Figure1/whitefly_no-choice_19_accessions.tsv", header = T)
wf$alive = NULL
wf$dead = NULL
wf$total = NULL
wf = wf %>% dplyr::group_by(accession) %>% dplyr::summarise(wf_average = mean(percentage,na.rm = T))


#####################################
# merge the datasets to a new one
#####################################
density.pheno  = left_join(density.avg, df.medians %>% select(., accession, median), by = "accession") %>% dplyr::rename(median_trips = median)
density.pheno  = left_join(density.pheno, wf, by = "accession")

# Label the trichome classes nicely
density.pheno$type = factor(density.pheno$type, levels = c("non.glandular", "typeIandIV", "typeVI"),
                            labels = c("Non glandular", "Type I/IV", "Type VI"))

#########
# plots #
#########

# Plot whiteflies vs. trichome densities
p.whitelfies = 
density.pheno %>%
ggplot() +
  geom_point(aes(x = (wf_average/max(density.pheno$wf_average))*100, y = density),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = (wf_average/max(density.pheno$wf_average))*100, y=density,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  facet_wrap(~type, scale = "free", ncol = 1)+
  geom_smooth(aes(x=wf_average/max(density.pheno$wf_average)*100, y = density), method = "lm", alpha = 0.2)+
  labs(x = "Relative whitefly survival",y = "Trichome density (trichomes / mm2 leaf)") +
  my.theme+
  ylim(0, NA)

# Plot thrips vs. trichome densities
p.thrips = 
  density.pheno %>%
  ggplot() +
  geom_point(aes(x = (median_trips/max(density.pheno$median_trips))*100, y = density),fill="grey",color="black",shape=21,size=2) +
  theme_bw() +
  geom_label_repel(aes(x = (median_trips/max(density.pheno$median_trips))*100, y=density,label=accession, fill=species),
                   size = 2,
                   label.size = 0.05,
                   label.padding = 0.1,
                   show.legend = FALSE) +
  facet_wrap(~type, scale = "free", ncol = 1)+
  geom_smooth(aes(x =(median_trips/max(density.pheno$median_trips))*100, y = density), method = "lm", alpha = 0.2)+
  labs(x = "Relative thrips survival",y = "Trichome density (trichomes / mm2 leaf)") +
  my.theme+
  ylim(0, NA)


#################
# Linear models #
#################

# whitefly vs (relative) type I/IV trichome density
lm.wf.typeIIV = summary( 
  lm(data = density.pheno %>% filter(., type == "Type I/IV"),
          (wf_average/max(density.pheno$wf_average)*100) ~ density)
)

# whitefly vs (relative) type VI trichome density
lm.wf.typeVI = summary(
  lm(data = density.pheno %>% filter(., type == "Type VI"),
          (wf_average/max(density.pheno$wf_average)*100) ~ density)
)

# whitefly vs (relative) type I/IV trichome density
lm.thrips.typeIIV = summary(
  lm(data = density.pheno %>% filter(., type == "Type I/IV"),
                   (median_trips/max(density.pheno$median_trips)*100) ~ density)
)

# whitefly vs (relative) type VI trichome density
lm.thrips.typeVI = summary(
  lm(data = density.pheno %>% filter(., type == "Type VI"),
                  (median_trips/max(density.pheno$median_trips)*100) ~ density)
)


capture.output(lm.wf.typeIIV, file = "Figure3/Linear models/linear_model_wf_typeI_IV")
capture.output(lm.wf.typeVI, file = "Figure3/Linear models/linear_model_wf_type_VI")
capture.output(lm.thrips.typeIIV, file = "Figure3/Linear models/linear_model_thrips_typeI_IV")
capture.output(lm.thrips.typeVI, file = "Figure3/Linear models/linear_model_thrips_type_VI")



################
# Save the plots
################
ggsave(filename = file.path("Figure3/whitefly_toxicity_vs_trichome_densities.pdf"),plot = p.whitelfies, width = 9,height = 21, units = "cm")
ggsave(filename = file.path("Figure3/thrips_toxicity_vs_trichome_densities.pdf"),plot = p.thrips, width = 9,height = 21, units = "cm")

#################
# session info
##################
writeLines(capture.output(sessionInfo()), "Figure3/sessionInfo.txt")
