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
survData = read.delim("Figure2_thrips/survival_data.tsv",header=T,stringsAsFactors = F)

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
  theme(axis.text.x = element_text(angle=90))

ggsave(filename = "Figure2_thrips/Figure2B.pdf",plot = p,width = 7,height = 5)
ggsave(filename = "Figure2_thrips/Figure2B.svg",plot = p,width = 7,height = 5)

######################
# Plot survival curves
######################

# add species information
colnames(df)[ncol(df)]="accession"
df = dplyr::left_join(df,accession2species,by="accession")

# reorder by increasing survival time
df$genotype = factor(x = df$genotype,levels = unique(df$genotype))

# plot
ggplot(data = df, aes(x = time, y = surv, color = color)) +
  geom_line() + 
  # plot censor marks
  geom_point(aes(shape = factor(ifelse(n.censor >= 1, 1, NA)))) + 
  # format censor shape as "+"
  scale_shape_manual(values = 3) + 
  # hide censor legend 
  guides(shape = "none") +
  facet_wrap(~ genotype) +
  xlim(0,20)


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
writeLines(capture.output(sessionInfo()), "Figure2_thrips/sessionInfo.txt")
