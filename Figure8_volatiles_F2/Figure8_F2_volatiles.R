# library
library(tidyverse)
library(superheat)

# load metabolites data
metabolites = read.delim("Figure8_volatiles_F2/20180625_F2metabolite_table.txt",header=T,stringsAsFactors = F,check.names = F)
head(metabolites)


# load whitefly data
wf = read.delim("Figure7_F2_whitefly_survival/20170912_wf_bioassay_F2.txt",header = T,stringsAsFactors = F)


####### plot heatmap ########

# calculate mean WF survival and order
wf = wf %>% mutate(percentage = alive/(alive+dead)*100)
wf = wf %>% group_by(line) %>% summarise(survival = mean(percentage)) 
wf =  wf[order(wf$survival),]

# reorder metabolites according to increasing wf survival values
metabolites.sorted = metabolites[match(wf$line,metabolites$line)]

# prepare matrix for heatmaps
