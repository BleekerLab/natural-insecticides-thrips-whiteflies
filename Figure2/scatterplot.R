################
# load libraries
################
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(survival)
library(survminer)


#####################
## Scatterplot (rank)
#####################
df = read.delim("Figure2/data4scatterplot_ranks.tsv",header=T,stringsAsFactors = F)

g <- ggplot(df) +
  geom_point(aes(wf,thrips),fill="grey",color="black",shape=21,size=4) +
  theme_bw() +
  geom_label_repel(aes(x = wf,y=thrips,label=sample,fill=species)) +
  labs(x = "Tomato genotype rank for whitefly survival (low to high survival)",y = "Tomato genotype rank for thrips survival (low to high survival)") +
  scale_x_continuous(breaks=seq(0,19,1)) +
  scale_y_continuous(breaks=seq(0,19,1))

ggsave(filename = file.path("Figure2/","scatterplot_rank.svg"),plot = g,width = 7,height = 5)
ggsave(filename = file.path("Figure2/","scatterplot_rank.png"),plot = g,width = 7,height = 5)

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
thrips_relative = mutate(thrips,scaled_thrips_survival = median * multiplication_coefficient)

##############
# Session Info
##############
writeLines(capture.output(sessionInfo()), "Figure2/sessionInfo.txt")