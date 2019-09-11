# Script to plot heatmaps of the metabolites plus the annotations of the phenotypes and metabolic classes

library(tidyverse)
library(pheatmap)

##############
# Acylsugars #
##############

#Load dataset and shape accession names
acylsugars = read.delim("Figure6_heatmaps/acylsugars_normalised.tsv", header = T) %>% filter(., accession != "LA2663-chmi")
acylsugars = separate(acylsugars, col = "accession", into = c("accession", "x")) %>% select(., -x) %>% column_to_rownames(., var = "accession")

# Change NA to 0 -> +1 and logtransfrom
acylsugars[is.na(acylsugars)] <- 0
acylsugars = acylsugars[] + 1
log.acylsugars = log10(acylsugars[,2:ncol(acylsugars)])
log.acylsugars = log.acylsugars %>% select(., -c("summed_glucose", "summed_sucrose", "summed_total")) #Remove the 'summed' columns


# Load phenotypes (i.e. toxic / non-toxic) and acylsugar annotations
phenolink = read.delim("Figure6_heatmaps/plots_metabolites_vs_phenotypes/phenolink_thrips_and_whitefly.txt", row.names = 1, header = T)
acylsugars.annotation = read.delim("Figure6_heatmaps/plots_metabolites_vs_phenotypes/acylsugars_annotation.tsv", row.names = 1, header = T)
row.names(acylsugars.annotation) = colnames(log.acylsugars)

# Heatmap including annotations
pheatmap(log.acylsugars, 
         annotation_row = phenolink, 
         annotation_col = acylsugars.annotation,
         fontsize = 6, 
         fontsize_row = 6, 
         fontsize_col = 6
         )

savePlot(filename = "/Users/Ruy/Desktop/heatmap.pdf", type = "pdf")

############# 
# Volatiles #
#############

#Load dataset and shape accession names
volatiles = read.delim("Figure6_heatmaps/pheno_terpenoids.tsv", header = T)
volatiles = separate(volatiles, sample, into = c("x", "y", "accession")) %>% select(., -c("x","y", "wf", "thrips")) %>% column_to_rownames(., var = "accession")

# Change NA to 0 -> +1 -> logtransfrom
volatiles[is.na(volatiles)] <- 0
volatiles = volatiles[] + 1
log.volatiles = log10(volatiles[,2:ncol(volatiles)])


# Load phenotypes (i.e. toxic / non-toxic) and volatile class-annotations
volatiles.annotation = read.delim(file = "Figure6_heatmaps/plots_metabolites_vs_phenotypes/volatiles_annotation.tsv", header = T, row.names = 1)

# Heatmap including annotations
pheatmap(log.volatiles,
         annotation_row = phenolink, 
         annotation_col = volatiles.annotation,
         fontsize = 6, 
         fontsize_row = 6, 
         fontsize_col = 6
)
