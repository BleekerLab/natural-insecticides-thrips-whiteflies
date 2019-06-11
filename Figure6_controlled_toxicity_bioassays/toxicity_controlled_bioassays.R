###########
# library
##########
library(tidyverse)
library("survival")        # to fit a survival curve
library("survminer")       # to draw plots
library("RColorBrewer")
library(gridExtra)         # to arrange plots in a grid

###########################################
# creates a directory where plots are saved
###########################################
dir.create(file.path("Figure6_controlled_toxicity_bioassays/plots/"), showWarnings = FALSE)

##################
# helper functions
##################

calculate_whitefly_percentage_of_survival <- function(df){
  df = dplyr::mutate(df,percentage = alive / (alive + dead)*100)
  return(df)
}

# whitefly plotting function
plot_wf_survival <- function(df,color){
  # reorder levels for x-axis ordering
  df$dose = factor(x = df$dose,levels = unique(df$dose))
  # make the plot
  p = ggplot(df,mapping = aes(x = dose,y = percentage)) +
    geom_boxplot(fill=color) +
    labs(x = "Dose",y="Whitefly survival after 48h (%)") +
    ggtitle(unique(df$compound)) + 
    theme(axis.text = element_text(colour = "black",size = 8),axis.text.x = element_text(angle = 0))
  return(p)
}

# thrips plotting function

plot_thrips_survival_curve <- function(fit,df,color){
  # reorder dose levels for the compound
  df$dose = factor(df$dose,levels = unique(df$dose))
  # plot
  p = ggsurvplot(fit,
                 data = df,
                 color = color,
                 conf.int=TRUE)
  g = p$plot + 
    facet_wrap(~strata,nrow = 1) +
    theme_bw() + 
    ggtitle(unique(df$compound)) + 
    theme(axis.title.y = element_blank(),legend.position = "none",text = element_text(colour="black"))
}


# plot_thrips_survival_curve <- function(df="thrips",compoundName,destDirectory){
#   # filters the complete dataset and create a survival curve
#   df = filter(thrips,compound == compoundName) 
#   # reorder dose levels for the compound
#   df$dose = factor(df$dose,levels = unique(df$dose))
#   
#   p = ggsurvplot(fit,data = df,
#                  palette=brewer.pal(n = length(unique(phe$dose)),"Set2"),
#                  conf.int=TRUE)
#   g = p$plot + theme_bw() + ggtitle(compoundName) + facet_wrap(~strata,nrow=1) + theme(axis.title.y = element_blank(),legend.position = "none")
#   # save svg plot
#   svg(filename = paste(destDirectory,compoundName,".svg",sep = ''),width = 7,height = 4)
#   print(g)
#   dev.off()
#   # save pdf plot
#   pdf(file = paste(destDirectory,compoundName,".pdf",sep = ''),width = 7,height = 4)
#   print(g)
#   dev.off()
# }
###############
# Load datasets
###############

# reads the whitefly data tables and place them in a list
wf.dfs = lapply(
  X = list.files(path = "Figure6_controlled_toxicity_bioassays",pattern = "*wf*",full.names = T),
  FUN = function(x){read.delim(x,header = T,stringsAsFactors = F,row.names = NULL)}
  )

# reads the thrips data tables and place them in a list
thrips.dfs = lapply(
  X = list.files(path = "Figure6_controlled_toxicity_bioassays",pattern = "*thrips*",full.names = T),
  FUN = function(x){read.delim(x,header = T,stringsAsFactors = F,row.names = NULL)}
  )
  
########################################################
# Boxplots for whitefly survival (end-point measurement)
########################################################

# calculate percentage of survival
wf.dfs = map(.x = wf.dfs,.f = calculate_whitefly_percentage_of_survival)

# make a set of colors
brewerColors = brewer.pal(n = length(wf.dfs),name = "Set1")

# create plots by using the puur map2 function
wf_plots  = map2(
  .x = wf.dfs,
  .y = brewerColors,
  .f = plot_wf_survival)

# use gridExtra to arrange them
nCol <- 3
pdf(file = file.path("Figure6_controlled_toxicity_bioassays/plots/Figure6A.pdf"),width = 10,height = 7)
do.call("grid.arrange", c(wf_plots, ncol=nCol))
dev.off()

svg(file = file.path("Figure6_controlled_toxicity_bioassays/plots/Figure6A.svg"),width = 10,height = 7)
do.call("grid.arrange", c(wf_plots, ncol=nCol))
dev.off()


##############################################################
# Survival curves for thrips survival (continous measurements)
##############################################################
# fit a survival curve for each dataset
fits = map(.x = thrips.dfs,
           ~ survfit(data=.x,formula = Surv(time,status) ~ dose,se.fit=TRUE)
           )

# make a set of colors
brewerColors = brewer.pal(n = length(thrips.dfs),name = "Set2")

thrips_plots = pmap(.l = list(fits,thrips.dfs,brewerColors),.f = plot_thrips_survival_curve)

# use gridExtra to arrange them
n <- length(thrips_plots)
nCol <- floor(sqrt(n))
pdf(file = file.path("Figure6_controlled_toxicity_bioassays/plots/Figure6B.pdf"),width = 10,height = 7)
do.call("grid.arrange", c(thrips_plots, ncol=nCol))
dev.off()

svg(file = file.path("Figure6_controlled_toxicity_bioassays/plots/Figure6B.svg"),width = 10,height = 7)
do.call("grid.arrange", c(thrips_plots, ncol=nCol))
dev.off()


# # fi,c  t a survival curve for each individual compound
# compounds = unique(thrips$compound)
# fits = list()
# compound.dfs=list()
# for (i in seq_along(compounds)){
#   metabolite = compounds[i]
#   # create a dataframe with only one compound
#   compound.dfs[[i]] = filter(thrips,compound == metabolite) 
#   # reorder dose levels for the compound
#   compound.dfs[[i]]$dose = factor(compound.dfs[[i]]$dose,levels = unique(compound.dfs[[i]]$dose))
#   # fit a survival curve on that single dataframe
#   df = compound.dfs[[i]] 
#   fits[[i]] = survfit(formula = Surv(time,status) ~ dose,se.fit=TRUE,data = df)
# 
#   # extract one fit
#   fit = fits[[i]]
#   # remove the 'dose=' label
#   names(fit$strata) = gsub(pattern = "dose=",replacement = "",x = names(fit$strata))
#   
#   # one plot per compound with all doses indicated
#   p = ggsurvplot(fit,
#                  data = df,
#                  palette=brewer.pal(n = length(unique(df$dose)),"Set2"),
#                  conf.int=TRUE)
#   p <- ggpar(p,legend = "none",ylab = "Survival probability",xlab = "Time (days)")
#   g <- p$plot + 
#     theme_bw() + 
#     ggtitle(compounds[[i]]) +
#     facet_wrap(~strata,nrow = 1) +
#     theme(axis.title.y = element_blank(),legend.position = "none")
#   # save svg plot
#   svg(filename = paste("Figure6_controlled_toxicity_bioassays/",metabolite,".svg",sep = ''),width = 7,height = 4)
#   print(g)
#   dev.off()
#   # save pdf plot
#   pdf(file = paste("Figure6_controlled_toxicity_bioassays/",metabolite,".pdf",sep = ''),width = 7,height = 4)
#   print(g)
#   dev.off()
# }

##### to fix issue: "f(...) cannot change aethetics while plotting ######
##### maybe it has to do with the confidence intervals...   #############

# get list of compounds
#compounds = unique(thrips$compound)



# make plots
#for (i in seq_along(compounds)){
#  plot_thrips_survival_curve(df = thrips,compoundName = compounds[i],destDirectory = "Figure6_controlled_toxicity_bioassays/plots")
#}


