################
# load libraries
################
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(survival)
library(survminer)

# custom functions for Restricted mean survival time (RMST) calculation
source("Figure_1C/rmst_functions.R")


#####################
## Scatterplot (rank)
#####################
df = read.delim("Figure_1C/data4scatterplot_ranks.tsv",header=T,stringsAsFactors = F)

g1 <- ggplot(df) +
  geom_point(aes(wf,thrips),fill="grey",color="black",shape=21,size=4) +
  theme_bw() +
  geom_label_repel(aes(x = wf,y=thrips,label=sample,fill=species)) +
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Tomato genotype rank for thrips survival (low to high survival)") +
  scale_x_continuous(breaks=seq(0,19,1)) +
  scale_y_continuous(breaks=seq(0,19,1))

ggsave(filename = file.path("Figure_1C/","Version1_scatterplot_rank.svg"),plot = g1,width = 7,height = 5)
ggsave(filename = file.path("Figure_1C/","Version1_scatterplot_rank.png"),plot = g1,width = 7,height = 5)

#########################
## Scatterplot (relative)
########################

#####
# WF
#####

### Data import and wrangling 
# import whitefly no-choice data
df = read.delim("./Figure_1C/whitefly_no-choice_19_accessions.tsv",header = T,stringsAsFactors = F)

# remove unecessary variables
df$alive = NULL
df$dead = NULL
df$total = NULL

# average clip-cages results
df = df %>% dplyr::group_by(accession) %>% dplyr::summarise(wf_average = mean(percentage,na.rm = T))

# import accession to species
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors=F)
df = dplyr::left_join(x = df,y = accession2species) %>% select(-genotype)

#Calculate relative WF survival
df$wf_relative_survival = (df$wf_average/max(df$wf_average))*100

########
# Thrips
########
survData = read.delim("Figure_1C/thrips_survival_data.tsv",header=T,stringsAsFactors = F)


###################
# Survival analysis
###################
fit <- with(survData,survfit(formula = Surv(time,status) ~ accession))

# extract accession order by increasing survival time to reorder factor
df.medians = surv_median(fit)
df.medians = mutate(df.medians,strata = gsub("accession=",replacement = "",strata))
colnames(df.medians)[1]="accession"

df = dplyr::left_join(df,df.medians,by="accession")
df$thrips_relative_survival = (df$median/max(df$median))*100
df$accession = factor(df.medians$accession,levels = df.medians$accession)


########
# Plot #
########

# Relative survival
g1 <- ggplot(df) +
  geom_point(
    aes(x = wf_relative_survival,y = thrips_relative_survival),
    fill="grey",color="black",shape=21,size=4) +
  theme_bw() +
  geom_label_repel(aes(x=wf_relative_survival,y=thrips_relative_survival,label=accession,fill=species)) +
  labs(x = "Tomato genotype whitefly relative survival",y = "Tomato genotype thrips relative survival")  + 
  scale_x_continuous(breaks=seq(0,100,20)) +  scale_y_continuous(breaks=seq(0,100,20))

# Survival
g2 = ggplot(df) +
  geom_point(
    aes(x = wf_average,y = median),
    fill="grey",color="black",shape=21,size=4) +
  theme_bw() +
  geom_label_repel(aes(x = wf_average,y = median,label=accession,fill=species)) +
  labs(x = "Tomato genotype whitefly survival (%)",y = "Tomato genotype thrips survival (median days)")+
  scale_x_continuous(breaks=seq(0,85,10)) +  scale_y_continuous(breaks=seq(0,20,5))

g3 = ggarrange(g1,g2,ncol = 2,nrow = 1,common.legend = TRUE)


### save plots
ggsave(filename = file.path("Figure_1C/","Version3_scatterplot_relative.svg"),plot = g1,width = 7,height = 5)
ggsave(filename = file.path("Figure_1C/","Version3_scatterplot_relative.png"),plot = g1,width = 7,height = 5)

ggsave(filename = file.path("Figure_1C/","Version3_scatterplot_survival.svg"),plot = g2,width = 7,height = 5)
ggsave(filename = file.path("Figure_1C/","Version3_scatterplot_survival.png"),plot = g2,width = 7,height = 5)

ggsave(filename = file.path("Figure_1C/","Version3_relative_vs_normal_survival.pdf"),plot = g3,width = 15,height = 6)


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
ggsave(filename = file.path("Figure_1C/","Version3_scatterplot_relative_rmst.svg"),plot = g3,width = 7,height = 5)
ggsave(filename = file.path("Figure_1C","Version3_scatterplot_relative_rmst.png"),plot = g3,width = 7,height = 5)


###################################
# save dataframe used for the plots
###################################

# rename some columns
colnames(df_for_relative_scatterplot_all)[1] = "wf_survival"
colnames(df_for_relative_scatterplot_all)[3] = "thrips_median_survival_time"

write.table(df_for_relative_scatterplot_all,
            file = "Figure_1C/dataframe_for_relative_scatterplots.tsv",
            quote = F,row.names = F,sep = "\t")

##############
# Session Info
##############
writeLines(capture.output(sessionInfo()), "Figure_1C/sessionInfo.txt")
