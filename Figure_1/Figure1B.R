#########
# Library
#########
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

library(survminer)
library(survival)
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('svglite'))
suppressPackageStartupMessages(library('ggfortify'))

#############
# Data import
#############
survData = read.delim("Figure_1/thrips_survival_data.tsv",header=T,stringsAsFactors = F)

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

##################################
# Plot barplot of median survivals
##################################
p <- ggplot(data = df.medians,
            aes(x = accession, y = median, fill = species)) +
  geom_bar(stat="identity",color="black") +
  labs(x = "genotype", y = "Median survival time (days)")
  theme_bw()+
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1),
        axis.text = element_text(colour = "black"))

ggsave(filename = "Figure_1/Fig1B_thrips.png", plot = p, width = 7, height = 3)
ggsave(filename = "Figure_1/Fig1B_thrips.svg", plot = p, width = 7, height = 3)

#################
# session info
##################
writeLines(capture.output(sessionInfo()), "Figure_1/sessionInfo.Figure1B.txt")
