#########
# Library
#########

### survminer
if (is.element('survminer', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('survminer'));
} else
{
  install.packages('survminer');
  suppressPackageStartupMessages(library('survminer'));
}

### survival 
if (is.element('survival', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('survival'));
} else
{
  install.packages('survival');
  suppressPackageStartupMessages(library('survival'));
}

### tidyverse
if (is.element('tidyverse', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('tidyverse'));
} else
{
  install.packages('tidyverse');
  suppressPackageStartupMessages(library('tidyverse'));
}

### svglite
if (is.element('svglite', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('svglite'));
} else
{
  install.packages('svglite');
  suppressPackageStartupMessages(library('svglite'));
}

### ggfortify
if (is.element('ggfortify', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('ggfortify'));
} else
{
  install.packages('ggfortify');
  suppressPackageStartupMessages(library('ggfortify'));
}

#############
# Data import
#############
survData = read.delim("Figure_1AandB/thrips_survival_data.tsv",header=T,stringsAsFactors = F)

# import accessions to species correspondence
accession2species = read.delim("genotype2species.txt",header = T,sep = "\t",stringsAsFactors = T)


###################
# Survival analysis
###################
fit <- with(survData,survfit(formula = Surv(time,status) ~ accession))

# extract result into a dataframe
df = fortify(fit)

# extract accession order by increasing survival time to reorder factor
df.medians = surv_median(fit)
df.medians = mutate(df.medians,strata = gsub("accession=",replacement = "",strata))
df.medians = df.medians[order(df.medians$median,decreasing = F),]
colnames(df.medians)[1]="accession"
df.medians = dplyr::left_join(df.medians,accession2species,by="accession")
df.medians$accession = factor(df.medians$accession,levels = df.medians$accession)

##################################
# Plot barplot of median survivals
##################################
p <- ggplot(data=df.medians,aes(x=accession,y=median,fill=species)) +
  geom_bar(stat="identity",color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1),
        axis.text = element_text(colour = "black"))

ggsave(filename = "Figure_1AandB/Figure1B.pdf",plot = p,width = 7,height = 3)
ggsave(filename = "Figure_1AandB/Figure1B.svg",plot = p,width = 7,height = 3)

######################
# Plot survival curves
######################

# add species information
colnames(df)[ncol(df)]="accession"
df = dplyr::left_join(df,accession2species,by="accession")

# reorder by increasing survival time
df$genotype = factor(x = df$genotype,levels = unique(df.medians$genotype))

# plot
# ggplot(data = df, aes(x = time, y = surv, color = color)) +
#   geom_line(size=1) + 
#   # plot censor marks
#   geom_point(aes(shape = factor(ifelse(n.censor >= 1, 1, NA))),size=3) + 
#   # format censor shape as "+"
#   scale_shape_manual(values = 3) + 
#   # hide censor legend 
#   guides(shape = "none") +
#   facet_wrap(~ genotype) +
#   xlim(0,20)


ggsurvplot_group_by(fit = with(survData,survfit(formula = Surv(time,status) ~ accession)),
           data = survData,group.by = c("accession")
           )

# # define helper function to plot survival curve
# plotSurvival = function(df,confidence=0.95,xmin=0,xmax = 19){
#   # calculate survival object
#   survObject = with(df,Surv(time = time,event = status))
#   # fit model
#   fit <- survfit(formula = survObject ~ accession,data = df,conf.int=confidence)
#   # plot model
#   g <- ggsurvplot(
#     fit,                      # survfit object with calculated statistics.
#     risk.table = FALSE,        # show risk table.
#     pval = TRUE,              # show p-value of log-rank test.
#     conf.int = TRUE,          # show confidence intervals for 
#     # point estimaes of survival curves.
#     xlim = c(xmin,xmax),           # present narrower X axis, but not affect
#     # survival estimates.
#     break.time.by = 1,        # break X axis in time intervals by 500.
#     ggtheme = theme_bw(),  # customize plot and risk table with a theme.
#     risk.table.y.text.col = T,  # colour risk table text annotations.
#     risk.table.y.text = FALSE, # show bars instead of names in text annotations in legend of risk table
#     font.main = 18,              # title font size
#     font.x = 16,                 # font x axis 
#     font.y = 16                 # font y axis
#   )
#   return(g)
# }
# 
# 
# # empty list for plots
# l = list()
# 
# # make plots
# genotypes = unique(df.with.species$genotype)
# for (i in seq_along(genotypes)){
#   tmp_df = dplyr::filter(df,genotype == genotypes[i])
#   l[[i]] = plotSurvival(tmp_df)
#   plot_title = genotypes[i]
#   gg = l[[i]]$plot + ggtitle(plot_title)
#   ggsave(filename = file.path("Figure_1AandB/plots/",paste(genotypes[i],".png",sep = "")),plot = gg,width = 7,height = 5,dpi = 400)
# }

###############
# Fit Cox model
###############
# fit a Cox PH model
fit.cox = coxph(formula = Surv(time,status) ~ accession -1 ,data = survData)

# Get coefficients and p-values
coefs = as.data.frame(unclass(summary(fit.cox))$coefficients)

colnames(coefs)=c("coef","exp(coef)","se(coef)","z","pval")
coefs$accession = row.names(coefs$accession);

coefs$accession = sub("^accession",replacement = "",x = row.names(coefs))

# which conditions have the significant pvalues?
coefs$text = NA
coefs[which(coefs$pval > 0.001),]$text <- "ns"
coefs[which(coefs$pval < 0.001),]$text <- "***"

# write to table
write.table(x = coefs,file = "TableS2_Cox/TableS2_Cox_model_coefs.tsv",quote = F,sep = "\t",row.names = F)

#################
# session info
##################
writeLines(capture.output(sessionInfo()), "Figure_1AandB/sessionInfo.txt")
