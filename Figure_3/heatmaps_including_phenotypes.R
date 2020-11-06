# Script to plot heatmaps of the metabolites plus the annotations of the phenotypes and metabolic classes

library(tidyverse)
library(magrittr)
library(pheatmap)
library(grDevices)

##############
# Acylsugars #
##############

# Load dataset and shape accession names
acylsugars = read.delim("Figure_3/acylsugars_normalised.tsv", header = T,check.names = F)

acylsugars = separate(acylsugars, col = "accession", into = c("accession", "x"))  %>% select(., -x) %>% column_to_rownames(., var = "accession")

# Change NA to 0 -> +1 and logtransfrom
acylsugars[is.na(acylsugars)] <- 0
acylsugars = acylsugars[] + 1
log.acylsugars = log10(acylsugars[,2:ncol(acylsugars)])
log.acylsugars = log.acylsugars %>% select(., -c("summed_glucose", "summed_sucrose", "summed_total")) #Remove the 'summed' columns

# Cluster accessions to phenotypes ("make the 4 quadrants")
log.acylsugars = log.acylsugars %>% rownames_to_column()

# Order accessions according to whitefly survival (low to high survival) as on Figure 1
log.acylsugars$rowname = factor(log.acylsugars$rowname, levels = c("LA0716" ,  "PI127826", "LA1777" ,  "LYC4"   ,  "PI134418", "LA1718" ,  "LA1954" ,  "LA2695"  , "LA4024"  , "LA2172"  , "LA1401" , 
                                                                   "LA0407" ,  "LA1578"  , "LA1364" ,  "LA2133" ,  "MM"   ,    "LA1840"  , "LA0735" ,  "LA1278"),
                                                                    ordered = TRUE) 


log.acylsugars = log.acylsugars[order(log.acylsugars$rowname),]
rownames(log.acylsugars) = log.acylsugars[,1]
log.acylsugars[,1] = NULL



# Load phenotypes (i.e. toxic / non-toxic) and acylsugar annotations
# phenolink = read.delim("Figure_3/phenolink_thrips_and_whitefly.tsv", row.names = 1, header = T)
acylsugars.annotation = read.delim("Figure_3/acylsugars_annotation.tsv", row.names = 1, header = T)
row.names(acylsugars.annotation) = colnames(log.acylsugars)

colors4annotation = list(
  Sugar.backbone = c("Sucrose" = "black", "Glucose" = "grey")
)


# Make and save the plot (heatmap including annotations)
pdf(file = "Figure_3/acylsugar_heatmap.pdf",width = 10,height = 5)

pheatmap(log.acylsugars,
         cluster_rows = FALSE,
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
volatiles = read.delim("Figure_3/terpenoids_normalized.tsv", header = T,check.names = F) %>% select(., -'13.604_105.0725_blank')
volatiles = separate(volatiles, sample, into = c("x", "y", "accession")) %>% select(., -c("x","y", "wf", "thrips")) 

volatiles.id <- read.delim("Figure_3/volatile_identifications_with_KI.tsv", sep = "\t") %>% select(metabolite, tentative_id)

volatiles.long <- pivot_longer(volatiles, -accession, names_to = "metabolite", values_to = "abundance")
volatiles.long$metabolite <- factor(volatiles.long$metabolite)
volatiles.long.id <- left_join(volatiles.long, volatiles.id, by = "metabolite") %>% 
  select(tentative_id, abundance, accession)

volatiles = pivot_wider(volatiles.long.id, values_from = "abundance", names_from = "tentative_id") %>%
  column_to_rownames(var = "accession")

# Change NA to 0 -> +1 -> logtransfrom
volatiles[is.na(volatiles)] <- 0
volatiles = as.data.frame(volatiles[] + 1)
log.volatiles = log10(volatiles[,2:ncol(volatiles)])

log.volatiles = log.volatiles %>% rownames_to_column()
log.volatiles$rowname = factor(log.volatiles$rowname, levels = c("LA0716" ,  "PI127826", "LA1777" ,  "LYC4"   ,  "PI134418", "LA1718" ,  "LA1954" ,  "LA2695"  , "LA4024"  , "LA2172"  , "LA1401" , 
                                                                 "LA0407" ,  "LA1578"  , "LA1364" ,  "LA2133" ,  "MM"   ,    "LA1840"  , "LA0735" ,  "LA1278"),
                               ordered = TRUE) 

log.volatiles = log.volatiles[order(log.volatiles$rowname),]
rownames(log.volatiles) = log.volatiles[,1]
log.volatiles[,1] = NULL


# Load phenotypes (i.e. toxic / non-toxic) and volatile class-annotations
volatiles.annotation = read.delim("Figure_3/volatile_identifications_with_KI.tsv", sep = "\t") %>% select(tentative_id, class)%>% 
  select(tentative_id, class) %>% column_to_rownames(var = "tentative_id")

colors4annotation.2 = list(
  class = c("monoterpene"="#999999", 
            "sesquiterpene"="#E69F00", 
            "aromatic hydrocarbon" = "#56B4E9", 
            "unknown"= "#009E73",
            "cyclic hydrocarbon" = "#0072B2", 
            "acyclic hydrocarbon" = "#D55E00", 
            "methylketone" = "#CC79A7")
  )


# Make and save the plot (heatmap including annotations)
pdf(file = "Figure_3/volatile_heatmap.pdf",width = 10,height = 5)

pheatmap(log.volatiles,
         cluster_rows = FALSE,
         annotation_col = volatiles.annotation,
         annotation_colors = colors4annotation.2,
         fontsize = 6, 
         fontsize_row = 6, 
         fontsize_col = 6
)

dev.off()