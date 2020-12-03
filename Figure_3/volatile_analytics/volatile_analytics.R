library(tidyverse)
library(pheatmap)

volatiles <- read_delim(file = "Figure_3/volatile_analytics/GCMS_normalised_counts.txt", delim = "\t")
volatiles$accession <- factor(volatiles$accession, levels = c("LA0716" ,  "PI127826", "LA1777" ,  "LYC4"   ,  "PI134418", "LA1718" ,  "LA1954" ,  "LA2695"  , "LA4024"  , "LA2172"  , "LA1401" , 
                                                              "LA0407" ,  "LA1578"  , "LA1364" ,  "LA2133" ,  "MM"   ,    "LA1840"  , "LA0735" ,  "LA1278"),
                              ordered = TRUE)
# Order volatiles to wf survival
volatiles = volatiles[order(volatiles$accession),]

# Load annotations for heatmap
volatile.annotation <- read_delim(file = "Figure_3/volatile_identifications_with_KI.tsv", delim = "\t") %>%
  select(metabolite, class) %>% column_to_rownames(var = "metabolite")

sample.annotation <- volatiles %>% select(sample, accession) %>% column_to_rownames(var = "sample")


# Prepare matrix for heatmap
mat <- volatiles %>% select(-accession) %>% column_to_rownames(var = "sample")
mat.log <- log(mat+1)

# Plot
pheatmap(mat.log,
         annotation_row = sample.annotation,
         annotation_col = volatile.annotation, 
         cluster_rows = F,
         cluster_cols = T)


####################
# plot proportions #
####################

volatiles.long <- pivot_longer(volatiles, -c(sample, accession), names_to = "metabolite", values_to = "abundance")
volatiles.long <- left_join(volatiles.long, (volatile.annotation %>% rownames_to_column(var = "metabolite")), by = "metabolite")

#######################
# Ordering of the plots
#######################

genotype_order_whiteflies = c("LA0716" ,  "PI127826", "LA1777" ,  "LYC4"   ,  "PI134418", "LA1718" ,  "LA1954" ,  "LA2695"  , "LA4024"  , "LA2172"  , "LA1401" , 
                              "LA0407" ,  "LA1578"  , "LA1364" ,  "LA2133" ,  "MM"   ,    "LA1840"  , "LA0735" ,  "LA1278")

genotype_order_thrips = c("LYC4","LA0407", "LA1777", "PI134418",
                          "LA1401", "LA0716", "LA1278",  "LA2172", 
                          "LA2695","LA0735","LA1718","LA2133","PI127826",
                          "LA1578", "LA1954", "MM",  "LA1840", "LA4024","LA1364")

volatiles.long$accession = factor(volatiles.long$accession,
                                 levels = genotype_order_whiteflies,
                                 ordered = TRUE)



# Manual colours for plot
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#0072B2", "#D55E00", "#CC79A7")

########
# Plot #
########

p.volatiles = 
volatiles.long %>%
  dplyr::group_by(accession, metabolite, class) %>%
  dplyr::summarise(mean_abundance = mean(abundance)) %>%
  dplyr::group_by(accession, class) %>%
  dplyr::summarise(sum_abundance = sum(mean_abundance)) %>%
  ggplot(aes(x = accession, y = sum_abundance, fill = class)) +
  geom_bar(position = "fill", stat = "identity")+
  theme_bw()+
  theme(legend.position  = "bottom",
        axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 10))+
  scale_fill_manual(values = cbp1)+
  labs(class = "Chemical class")+
  ylab("Proportion of total volatiles")

################
# Plot stacked #
################

p.volatiles.abundance = 
  volatiles.long %>%
  dplyr::group_by(accession, metabolite, class) %>%
  dplyr::summarise(mean_volatiles = mean(abundance)) %>%
  dplyr::group_by(accession, class) %>%
  dplyr::summarise(mean_abundance = mean(log(mean_volatiles+1)),
                   se = sd(log(mean_volatiles+1))/sqrt(n())) %>%
  ggplot(aes(x = accession, y = mean_abundance, fill = class)) +
  geom_bar(position = "stack", stat = "identity")+
  geom_errorbar(aes(x = accession,
                    ymin = mean_abundance - se,
                    ymax = mean_abundance + se),
                width = 0.5)+
  theme_bw()+
  theme(legend.position  = "bottom",
        axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 10))+
  scale_fill_manual(values = cbp1)+
  labs(class = "Chemical class")+
  ylab("Proportion of total volatiles")

ggsave(filename = "Figure_3/volatile analytics/stacked_volatile_proportions.pdf", plot = p.volatiles, height = 5.5, width = 6)
ggsave(filename = "Figure_3/volatile analytics/stacked_volatile_abundance.pdf", plot = p.volatiles.abundance, height = 5.5, width = 6)

####################################
# Create df of chemical structures #
####################################

#calculate means of metabolites over the samples
df <- volatiles.long %>%
  dplyr::group_by(accession, metabolite, class) %>%
  dplyr::summarise(mean_abundance = mean(abundance))

# Sum the metabolites from the same class
df.sum <- df %>%
  dplyr::group_by(accession, class) %>%
  dplyr::summarise(summed_abundance = sum(mean_abundance))

# Create a wide df
df.sum.wide <- df.sum %>% pivot_wider(names_from = class,
                                      values_from = summed_abundance) %>%
  rename(sample = accession)

df.sum.wide$total_volatiles <- rowSums(df.sum.wide[,2:8])

# remove "unknown" column
df.sum.wide <- df.sum.wide %>% select(-unknown)
df.sum.wide[is.na(df.sum.wide)] <- 0

# Load dataframe with phenotypes
phenotype <- read.delim(file = "Random_Forest/phenotypes_vs_leaf_terpenoids.tsv", check.names = FALSE) %>%
  select(sample, wf, thrips) %>%
  mutate(sample = substr(sample, start = 7, stop = 15)) 

df.final <- left_join(phenotype, df.sum.wide, by = "sample")

write.table(df.final, file = "Figure_3/volatile analytics/volatiles_summary.tsv",
            col.names = TRUE, row.names = FALSE, sep = "\t")
