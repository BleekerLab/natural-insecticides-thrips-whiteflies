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

##################################
# Plot barplot of median survivals
##################################
p <- ggplot(data=df.medians,aes(x=accession,y=median,fill=species)) +
  geom_bar(stat="identity",color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1),
        axis.text = element_text(colour = "black"))

ggsave(filename = "Figure1/Figure1B.pdf",plot = p,width = 7,height = 3)
ggsave(filename = "Figure1/Figure1B.svg",plot = p,width = 7,height = 3)

#################
# session info
##################
writeLines(capture.output(sessionInfo()), "Figure1/sessionInfo.Figure1B.txt")
