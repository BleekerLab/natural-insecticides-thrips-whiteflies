library(tidyverse)

# Check what lines have been screened
wf.lines <- read.delim("F2_Random_Forest/whitefly/20170912_wf_bioassay_F2.tsv",header = T,stringsAsFactors = F) %>% select(line)
toxic.lines <- read.delim("F2_Random_Forest/whitefly/toxic_lines.tsv", header = T, stringsAsFactors = F) %>% select(line)
metabolites.lines <- read.delim("F2_Random_Forest/volatiles/20180625_F2metabolite_table.tsv",header=T,stringsAsFactors = F,check.names = F) %>% select(line)

# load metabolites data
metabolites <- read.delim("F2_Random_Forest/volatiles/20180625_F2metabolite_table.tsv",header=T,stringsAsFactors = F,check.names = F) %>%
  filter(line %in% wf.lines$line)

# load whitefly data
wf <- read.delim("F2_Random_Forest/whitefly/20170912_wf_bioassay_F2.tsv",header = T,stringsAsFactors = F) %>%
  filter(line %in% metabolites.lines$line) %>%
  mutate(survival = ((alive/(alive+dead))*100)) %>%
  mutate(phenotype = if_else(line %in% toxic.lines$line, "Resistant", "Susceptible")) %>%
  arrange(survival) 
 

# Boxplot WF survival

ggplot(wf, aes(x = reorder(line, survival), y = survival, fill = phenotype)) +
  geom_boxplot() +
  geom_point(colour = "black", size = 1)+
  theme_bw() +
  scale_color_grey()+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust = 1, colour = "black", size = 8),
        axis.text.y = element_text(colour = "black", size = 8),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8))+
  xlab(NULL) +
  ylab("Whitefly survival (%)")

# Heatmap

# Calculate AVG WF survival for ordering the heatmap
wf.avg <- wf %>% 
  group_by(line) %>%
  dplyr::summarise(average_survival = mean(survival),
                   se = sd(survival)/sqrt(length(survival))) %>%
  arrange(average_survival)


# Order metabolites to ascending wf survival 
metabolites$line <- factor(metabolites$line, levels = wf.avg$line, ordered = TRUE)

# reorder metabolites according to increasing wf survival values
mat <- metabolites %>% column_to_rownames(var = "line")
mat.log <- log(mat+1)

pheatmap(mat.log,
         cluster_cols = T)
