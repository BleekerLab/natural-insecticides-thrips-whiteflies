################
# load libraries
################
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(survival)
library(survminer)

# custom functions for Restricted mean survival time (RMST) calculation
source("Figure2/rmst_functions.R")


#####################
## Scatterplot (rank)
#####################
df = read.delim("Figure2/data4scatterplot_ranks.tsv",header=T,stringsAsFactors = F)

g1 <- ggplot(df) +
  geom_point(aes(wf,thrips),fill="grey",color="black",shape=21,size=4) +
  theme_bw() +
  geom_label_repel(aes(x = wf,y=thrips,label=sample,fill=species)) +
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Tomato genotype rank for thrips survival (low to high survival)") +
  scale_x_continuous(breaks=seq(0,19,1)) +
  scale_y_continuous(breaks=seq(0,19,1))

ggsave(filename = file.path("Figure2/","Version1_scatterplot_rank.svg"),plot = g1,width = 7,height = 5)
ggsave(filename = file.path("Figure2/","Version1_scatterplot_rank.png"),plot = g1,width = 7,height = 5)

#########################
## Scatterplot (relative)
########################

#####
# WF
#####

### Data import and wrangling 
# import whitefly no-choice data
df = read.delim("./Figure1/whitefly_no-choice_19_accessions.tsv",header = T,stringsAsFactors = F)

# remove unecessary variables
df$alive = NULL
df$dead = NULL
df$total = NULL

# average clip-cages results
df = df %>% group_by(accession,plant) %>% summarise(average = mean(percentage,na.rm = T))

# import accession to species
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors=F)
df = dplyr::left_join(x = df,y = accession2species)

### Extract the survival percentages
wf = df %>% select(accession,plant,species,color,average) %>% group_by(accession) %>% summarise(perc_survival=mean(average)) %>% arrange(perc_survival)

### Scale from 0 to 100% (the high survival percentage becomes 100%)
multiplication_coefficient = 100 / max(wf$perc_survival)
wf_relative = mutate(wf,scaled_wf_survival = perc_survival * multiplication_coefficient)

########
# Thrips
########
survData = read.delim("Figure1/thrips_survival_data.tsv",header=T,stringsAsFactors = F)

# import accessions to species correspondence
accession2species = read.delim("genotype2species.txt",header = T,sep = "\t",stringsAsFactors = T)


###################
# Survival analysis
###################
fit <- with(survData,survfit(formula = Surv(time,status) ~ accession))

# extract accession order by increasing survival time to reorder factor
df.medians = surv_median(fit)
df.medians = mutate(df.medians,strata = gsub("accession=",replacement = "",strata))
df.medians = df.medians[order(df.medians$median,decreasing = F),]
colnames(df.medians)[1]="accession"
df.medians = dplyr::left_join(df.medians,accession2species,by="accession")
df.medians$accession = factor(df.medians$accession,levels = df.medians$accession)

### Extract the median survival time
thrips = df.medians %>% select(accession,median,species,color) 

### Scale from 0 to 100% (the high survival percentage becomes 100%)
multiplication_coefficient = 100 / max(thrips$median)
thrips_relative_median_survival = mutate(thrips,
                                         scaled_thrips_survival_median_time = median * multiplication_coefficient)

### Scatterplot using relative survival numbers
df_for_relative_scatterplot = inner_join(wf_relative,thrips_relative_median_survival,by="accession")

g2 <- ggplot(df_for_relative_scatterplot) +
  geom_point(
    aes(x = scaled_wf_survival,y = scaled_thrips_survival_median_time),
    fill="grey",color="black",shape=21,size=4) +
  theme_bw() +
  geom_label_repel(aes(x=scaled_wf_survival,y=scaled_thrips_survival_median_time,label=accession,fill=species)) +
  labs(x = "Tomato genotype whitefly survival (relative survival,%)",y = "Tomato genotype thrips survival (relative survival,%)")  + 
  scale_x_continuous(breaks=seq(0,100,20)) +  scale_y_continuous(breaks=seq(0,100,20))
g2

### save plots
ggsave(filename = file.path("Figure2/","Version2_scatterplot_relative.svg"),plot = g2,width = 7,height = 5)
ggsave(filename = file.path("Figure2/","Version2_scatterplot_relative.png"),plot = g2,width = 7,height = 5)

#####################################
# RMST: Restricted Mean Survival Time 
# The RMST is the mean survival time of all subjects in the study population followed up to t, and is simply the area under the survival curve up to t. (Zhao et al., 2016 Biometrics)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5114026/
#####################################

### calculate RMST based on the Kaplan-Meier survival curve
rmst.results = rmstKM(formula = Surv(time, status) ~ accession, 
       data=survData, 
       trunc=5, # end point integration.  The truncation time must be shorter than the minimum of the largest observed time in each group: 5.000
       alpha=0.05) # The default is 0.05. (1-alpha) confidence intervals are reported. https://cran.r-project.org/web/packages/survRM2/survRM2.pdf

rmst.df = as.data.frame(rmst.results$RMST)
genotypes = gsub(pattern = "RMST ",x = row.names(rmst.df),replacement = "")  
genotypes = gsub(pattern = ":",x = genotypes,replacement = "")  
rmst.df$accession = genotypes
colnames(rmst.df)[1]="rmst"


### Scale from 0 to 100% (the high survival percentage becomes 100%)
multiplication_coefficient = 100 / max(rmst.df$rmst)
thrips_relative_rmst = mutate(rmst.df,scaled_thrips_survival_rmst = rmst * multiplication_coefficient)

## Add to dataframe with relative data
thrips_relative_rmst = 
  mutate(rmst.df,
         scaled_thrips_survival_rmst = rmst * multiplication_coefficient) %>%
  select(accession,
         rmst,
         scaled_thrips_survival_rmst)

df_for_relative_scatterplot_all = inner_join(
  df_for_relative_scatterplot,
  thrips_relative_rmst,
  by="accession")

# scatterplot with RMST relative %
g3 <- ggplot(df_for_relative_scatterplot_all) +
  geom_point(
    aes(x = scaled_wf_survival,y = scaled_thrips_survival_rmst),
    fill="grey",color="black",shape=21,size=4) +
  theme_bw() +
  geom_label_repel(aes(x=scaled_wf_survival,y=scaled_thrips_survival_rmst,label=accession,fill=species)) +
  labs(x = "Tomato genotype whitefly survival (relative survival,%)",y = "Tomato genotype thrips survival (relative survival,%)")  + 
  scale_x_continuous(breaks=seq(0,100,20)) +  scale_y_continuous(breaks=seq(0,100,20))
g3

### save plots
ggsave(filename = file.path("Figure2/","Version3_scatterplot_relative_rmst.svg"),plot = g3,width = 7,height = 5)
ggsave(filename = file.path("Figure2/","Version3_scatterplot_relative_rmst.png"),plot = g3,width = 7,height = 5)


###################################
# save dataframe used for the plots
###################################

# rename some columns
colnames(df_for_relative_scatterplot_all)[1] = "wf_survival"
colnames(df_for_relative_scatterplot_all)[3] = "thrips_median_survival_time"

write.table(df_for_relative_scatterplot_all,
            file = "Figure2/dataframe_for_relative_scatterplots.tsv",
            quote = F,row.names = F,sep = "\t")

##############
# Session Info
##############
writeLines(capture.output(sessionInfo()), "Figure2/sessionInfo.txt")