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

p <- ggplot(wf,mapping = aes(x = dose,y = percentage,fill=compound)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  facet_wrap(. ~ compound,scales="free_x",nrow = 2) +
  labs(x = "Dose",y="Whitefly survival after 48h (%)") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_text(colour = "black",size = 8),axis.text.x = element_text(angle = 0))

svg(filename = "Figure6_controlled_toxicity_bioassays/Figure6A.svg",width = 10,height = 5)
print(p)
dev.off()

pdf(file = "Figure6_controlled_toxicity_bioassays/Figure6A.pdf",width = 10,height = 5)
print(p)
dev.off()

##############################################################
# Survival curves for thrips survival (continous measurements)
##############################################################
thrips = ldply(thrips.dfs,data.frame)

# fit a survival curve for each individual compound
compounds = unique(thrips$compound)
fits = list()
compound.dfs=list()
for (i in seq_along(compounds)){
  metabolite = compounds[i]
  # create a dataframe with only one compound
  compound.dfs[[i]] = filter(thrips,compound == metabolite) 
  # reorder dose levels for the compound
  compound.dfs[[i]]$dose = factor(compound.dfs[[i]]$dose,levels = unique(compound.dfs[[i]]$dose))
  # fit a survival curve on that single dataframe
  fits[[i]] = with(compound.dfs[[i]],survfit(formula = Surv(time,status) ~ dose,se.fit=T))
}

# plots
for (i in seq_along(compounds)){
  metabolite = compounds[i] # name of the compound
  # extract one fit
  fit = fits[[i]]
  # remove the 'dose=' label
  names(fit$strata) = gsub(pattern = "dose=",replacement = "",x = names(fit$strata))
  # a dataframe with one single ccompound
  df = compound.dfs[[i]] 
  # one plot per compound with all doses indicated
  p = ggsurvplot(fit,
                 data = df,
                 size=1,
                 palette=brewer.pal(n = length(unique(df$dose)),"Set2"),
                 conf.int=TRUE,
                 ggtheme = theme_bw()
                 ) 
  p <- ggpar(p,legend = "none",ylab = "Survival probability",xlab = "Time (days)")
  g <- p$plot + 
    theme_bw() + 
    ggtitle(compounds[[i]]) +
    facet_wrap(~strata,nrow = 1) +
    theme(axis.title.y = element_blank(),legend.position = "none")
  # save svg plot
  svg(filename = paste("Figure6_controlled_toxicity_bioassays/",metabolite,".svg",sep = ''),width = 7,height = 4)
  print(g)
  dev.off()
  # save pdf plot
  pdf(file = paste("Figure6_controlled_toxicity_bioassays/",metabolite,".pdf",sep = ''),width = 7,height = 4)
  print(g)
  dev.off()
}


