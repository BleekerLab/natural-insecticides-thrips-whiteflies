library(tidyverse)
library(pheatmap)

volatiles <- read_delim(file = "Figure_3/volatile analytics/GCMS_normalised_counts.txt", delim = "\t")
volatiles$accession <- factor(volatiles$accession, levels = c("LA0716","PI127826","LYC4", "LA1777", 
                                                                      "PI134418", "LA1718", "LA1954","LA2695",
                                                                      "LA1401","LA0407","LA1364","LA4024", "LA2172", 
                                                                      "LA2133","LA1578","LA0735", "LA1278","MM","LA1840"),
                              ordered = TRUE)
# Order volatiles to wf survival
volatiles = volatiles[order(volatiles$accession),]

# Load annotations for heatmap
volatile.annotation <- read_delim(file = "Figure_3/volatile analytics/volatile_identifications.txt", delim = "\t") %>%
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

genotype_order_whiteflies = c("LA0716","PI127826","LYC4", "LA1777", 
                              "PI134418", "LA1718", "LA1954","LA2695",
                              "LA1401","LA0407","LA1364","LA4024", "LA2172", 
                              "LA2133","LA1578","LA0735", "LA1278","MM","LA1840")

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

volatiles.long %>%
  dplyr::group_by(accession, class) %>%
  dplyr::summarise(mean_abundance = mean(abundance)) %>%
  ggplot(aes(x = accession, y = mean_abundance, fill = class)) +
  geom_bar(position = "fill", stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1))+
  scale_fill_manual(values = cbp1)+
  labs(class = "Chemical class")+
  ylab("Proportion of total volatiles")
