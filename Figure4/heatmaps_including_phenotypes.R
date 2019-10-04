# Script to plot heatmaps of the metabolites plus the annotations of the phenotypes and metabolic classes

library(tidyverse)
library(magrittr)
library(pheatmap)
library(grDevices)

##############
# Acylsugars #
##############

# Load dataset and shape accession names
acylsugars = read.delim("Figure4/acylsugars_normalised.tsv", header = T,check.names = F)

acylsugars = separate(acylsugars, col = "accession", into = c("accession", "x"))  %>% select(., -x) %>% column_to_rownames(., var = "accession")

# Change NA to 0 -> +1 and logtransfrom
acylsugars[is.na(acylsugars)] <- 0
acylsugars = acylsugars[] + 1
log.acylsugars = log10(acylsugars[,2:ncol(acylsugars)])
log.acylsugars = log.acylsugars %>% select(., -c("summed_glucose", "summed_sucrose", "summed_total")) #Remove the 'summed' columns


# Load phenotypes (i.e. toxic / non-toxic) and acylsugar annotations
phenolink = read.delim("Figure4/phenolink_thrips_and_whitefly.tsv", row.names = 1, header = T)
acylsugars.annotation = read.delim("Figure4/acylsugars_annotation.tsv", row.names = 1, header = T)
row.names(acylsugars.annotation) = colnames(log.acylsugars)

colors4annotation = list(
  "Sugar backbone" = c(
    "Sucrose" = "#7570B3",
    "Glucose" = "#E7298A"),
  Thrips.phenotype = c("Toxic" = "black","Non-toxic" ="white"),
  Whitefly.phenotype = c("Toxic" = "black","Non-toxic" ="white")
)


# Make and save the plot (heatmap including annotations)
pdf(file = "Figure4/acylsugar_heatmap.pdf",width = 10,height = 5)

pheatmap(log.acylsugars, 
         annotation_row = phenolink, 
         annotation_col = acylsugars.annotation,
         annotation_colors = colors4annotation,
         fontsize = 6, 
         fontsize_row = 6, 
         fontsize_col = 6
         )

dev.off()


############# 
# Volatiles #
#############

#Load dataset and shape accession names
volatiles = read.delim("Figure4/pheno_terpenoids.tsv", header = T,check.names = F)
volatiles = separate(volatiles, sample, into = c("x", "y", "accession")) %>% select(., -c("x","y", "wf", "thrips")) %>% column_to_rownames(., var = "accession")

# Change NA to 0 -> +1 -> logtransfrom
volatiles[is.na(volatiles)] <- 0
volatiles = volatiles[] + 1
log.volatiles = log10(volatiles[,2:ncol(volatiles)])


# Load phenotypes (i.e. toxic / non-toxic) and volatile class-annotations
volatiles.annotation = read.delim(file = "Figure4/volatiles_annotation.tsv", header = T, row.names = 1,check.names = F)

colors4annotation = list(
  metabolite_class = c(
    "methylketon" = "#7570B3",
    "monoterpene" = "#E7298A",
    "sesquiterpene" ="#66A61E",
    "other_carbohydrate" = "firebrick"),
  Thrips.phenotype = c("Toxic" = "black","Non-toxic" ="white"),
  Whitefly.phenotype = c("Toxic" = "black","Non-toxic" ="white")
)


# Make and save the plot (heatmap including annotations)
pdf(file = "Figure4/volatile_heatmap.pdf",width = 10,height = 5)

pheatmap(log.volatiles,
         annotation_row = phenolink, 
         annotation_col = volatiles.annotation,
         annotation_colors = colors4annotation,
         fontsize = 6, 
         fontsize_row = 6, 
         fontsize_col = 6
)

dev.off()
