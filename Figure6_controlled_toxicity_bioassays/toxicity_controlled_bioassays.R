###########
# library
##########
library(plyr) # load it before dplyr (in tidyverse)
library(tidyverse)
library("survival") # to fit a survival curve
library("survminer") # to draw plots
library("RColorBrewer")


###############
# Load datasets
###############
list_of_wf_assays = list.files(path = "Figure6_controlled_toxicity_bioassays",pattern = "*wf*",full.names = T)
list_of_thrips_assays = list.files(path = "Figure6_controlled_toxicity_bioassays",pattern = "*thrips*",full.names = T)

wf.dfs = lapply(list_of_wf_assays,FUN = function(x){read.delim(x,header = T,stringsAsFactors = F,row.names = NULL)})
thrips.dfs = lapply(list_of_thrips_assays,FUN = function(x){read.delim(x,header = T,stringsAsFactors = F,row.names = NULL)})

########################################################
# Boxplots for whitefly survival (end-point measurement)
########################################################
wf = ldply(wf.dfs,data.frame)
wf$percentage = with(wf,alive / (alive + dead) * 100)

wf$dose = factor(x = wf$dose,levels = unique(wf$dose))

p1 <- ggplot(wf,mapping = aes(x = dose,y = percentage,fill=compound)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  facet_grid(cols = vars(compound)) +
  labs(x = "Dose",y="Whitefly survival after 48h (%)") +
  scale_fill_brewer(palette = "OrRd")
p1


##############################################################
# Survival curves for thrips survival (continous measurements)
##############################################################
thrips = ldply(thrips.dfs,data.frame)

# reorder levels for dose
thrips$dose = factor(thrips$dose,levels = unique(thrips$dose))

# fit a survival curve for each individual compound
compounds = unique(thrips$compound)
fits = list()
compound.dfs=list()
for (i in seq_along(compounds)){
  metabolite = compounds[i]
  compound.dfs[[i]] = filter(thrips,compound == metabolite) # a dataframe for one compound
  fits[[i]] = with(compound.dfs[[i]],survfit(formula = Surv(time,status) ~ dose,se.fit=T))
}

# plots
list_of_plots = list()
for (i in seq_along(compounds)){
  fit = fits[[i]]
  df = compound.dfs[[i]] # a dataframe for one compound
  p = ggsurvplot(fit,
                 data = df,
                 size=1,
                 palette=brewer.pal(n = length(unique(df$dose)),"Set2"),
                 conf.int=TRUE,
                 ggtheme = theme_bw()
                 ) 
  g <- p$plot + theme_bw() + facet_grid(~strata)
  list_of_plots[[i]]=g
}


