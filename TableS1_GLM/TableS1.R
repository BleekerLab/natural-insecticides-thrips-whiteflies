#########
# Library
#########
suppressMessages(library(lme4))

##########################
# Generalized Linear Model
##########################

#############
# Data import
#############

# import accession to species
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors=F)

# remove missing values
# convert cage to factor
# Calculate a probability to die for each observation
df = read.delim("Figure1/whitefly_no-choice_19_accessions.tsv",header = T,stringsAsFactors = F)
df = na.omit(df)

# Remove LA0716 since it has only zero values (cannot estimate the coefficient)
df = dplyr::filter(.data = df,accession != "LA0716")

# Fit models
# dead / total - dead
fit = glm(cbind(dead,total-dead) ~ 1 + accession, data=df, family = binomial(link = logit)) # model with accession

# extract coefficients
coefs = as.data.frame(summary(fit)$coefficients)
colnames(coefs)=c("coeff","stderr","zval","pval")
coefs$accession = row.names(coefs$accession)
coefs$accession = sub("^accession",replacement = "",x = row.names(coefs))
coefs$accession[1] = "Intercept"

# add species
coefs = dplyr::left_join(coefs,accession2species)

# which conditions have the significant pvalues?
coefs$text = NA
coefs[which(coefs$pval > 0.05),]$text <- "ns"
coefs[which(coefs$pval < 0.05),]$text <- "***"

##########
# Table S1
##########
write.table(x = coefs,file = "TableS1_GLM/TableS1.tsv",sep = "\t",quote = F,row.names = F)

