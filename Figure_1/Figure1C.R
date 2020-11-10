################
# load libraries
################
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
checkpoint("2019-10-01")

library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(survival)
library(survminer)

#########################
## Scatterplot (relative)
########################

#####
# WF
#####

# Data import and wrangling 
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors=F)

# import whitefly no-choice data
# average clip-cages results
# then average survival over the plant replicates
# add accession to species
# Calculate relative survival based on the most susceptible accession
df = read.delim("Figure_1/whitefly_no-choice_19_accessions.tsv",
                header = T,
                stringsAsFactors = F) %>% 
  select(- alive, - dead, - total) %>%  
  dplyr::group_by(accession, plant) %>%  
  dplyr::summarise(plant_average = mean(percentage,na.rm = T)) %>% 
  dplyr::group_by(accession) %>%  
  dplyr::summarise(wf_average = mean(plant_average,na.rm = T)) %>%
  dplyr::left_join(., y = accession2species) %>%
  select(- genotype) %>% 
  mutate(wf_relative_survival = wf_average/max(wf_average) * 100) # Calculate relative WF survival

########
# Thrips
########

survData = read.delim("Figure_1/thrips_survival_data.tsv",header=T,stringsAsFactors = F)

fit <- with(survData,survfit(formula = Surv(time,status) ~ accession))

# extract accession order by increasing survival time to reorder factor
df.medians = surv_median(fit)
df.medians = mutate(df.medians,strata = gsub("accession=",replacement = "",strata))
colnames(df.medians)[1]="accession"

df = dplyr::left_join(df,df.medians,by="accession") 
df$thrips_relative_survival = (df$median/max(df$median))*100
df$accession = factor(df.medians$accession,levels = df.medians$accession)

df.medians = surv_median(fit)
df.medians = mutate(df.medians,strata = gsub("accession=",replacement = "",strata))
colnames(df.medians)[1]="accession"

########
# Plot #
########

# Relative survival
g1 <- ggplot(df) +
  geom_point(
    aes(x = wf_relative_survival, y = thrips_relative_survival),
    fill = "grey", color = "black", shape = 21, size = 4) +
  theme_bw() +
  geom_label_repel(aes(x = wf_relative_survival,
                       y = thrips_relative_survival,
                       label = accession,
                       fill = species)) +
  labs(x = "Tomato genotype whitefly relative survival", 
       y = "Tomato genotype thrips relative survival")  + 
  scale_x_continuous(breaks = seq(0,100,20)) +  
  scale_y_continuous(breaks = seq(0,100,20))+
  theme(legend.position = "none")


### save plots
ggsave(filename = "Figure_1/Figure_1C_scatterplot_relative.pdf", 
       plot = g1, 
       width = 4.5,
       height = 4.5)

ggsave(filename = "Figure_1/Figure_1C_scatterplot_relative.png", 
       plot = g1, 
       width = 7, 
       height = 5)


###################################
# save dataframe used for the plots
###################################

write.table(df,
            file = "Figure_1/dataframe_for_relative_scatterplots.tsv",
            quote = F,
            row.names = F,
            sep = "\t")

##############
# Session Info
##############
writeLines(capture.output(sessionInfo()), "Figure_1/sessionInfo.txt")
