library(mixOmics)
library(mdatools)
library(tidyverse)

# Read data - select the phenotype - create matrix for sPLS-DA
df = read.delim(file = "Random_Forest/phenotypes_vs_leaf_terpenoids.tsv", check.names = FALSE)

# Define the phenotype
pheno <- df$thrips

# Create numeric matrix for sPLS-DA
matrix <- df %>% select(-thrips, -wf ) %>% column_to_rownames(var = "sample")
matrix = log(matrix+1)

# Check model (e.g. number of components to include) with lowest error rate manually  
plsda.results <- splsda(matrix, pheno, ncomp = 8)
set.seed(30)
plsda.performance <- perf(plsda.results, validation = "Mfold", folds = 3, 
     progressBar = FALSE, nrepeat = 10)
plot(plsda.performance, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")


#############################
# Tune the plsda parameters #
#############################

list.keepX <- c(1:8,  seq(20, 100, 10))
tune.plsda<- tune.splsda(matrix, pheno, ncomp = 8, 
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)

error <- tune.plsda$error.rate
ncomp <- tune.plsda$choice.ncomp$ncomp

select.keepX <- tune.plsda$choice.keepX[1:2]  # optimal number of variables to select
select.keepX

plot(tune.plsda)

##################
# Do the sPLS-DA #
##################

plsda.results <- splsda(matrix, pheno, ncomp = 3) # use number of componentsa with lowest error rate as checked before

# plot the samples and loadings
plotIndiv(plsda.results, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")

plotVar(plsda.results, cutoff = 0.7)

# list the metabolites and their scores. Select the top 5 as "candidates"
scores <- selectVar(plsda.results, comp=1)$value
candidates <- rownames(scores)[1:5]



#######################
# Plot the candidates #
#######################

volatiles <- read.delim("Figure_4/leaf_terpenoids_normalised_peak_area.tsv", header = T, 
                     stringsAsFactors = TRUE, check.names = F)

volatiles.long <- gather(volatiles, 
                        key = "metabolite",
                        value = "abundance",
                        -sample, 
                        -accession)

# Filter metabolites for canditates
volatiles.long.candidates <- volatiles.long %>% filter(metabolite %in% candidates)

# Add accession species information 
accession2species <- read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
volatiles.long.candidates.species <- left_join(volatiles.long.candidates, accession2species, by = "accession")

# Add volatile identifications ("tentative_id")
volatiles.id <- read.delim(file = "Figure_3/volatile_identifications_with_KI.tsv") %>%
  select(metabolite, tentative_id)

volatiles.long.candidates.species <- left_join(volatiles.long.candidates.species, volatiles.id,
                                               by = "metabolite")


#########################
# Settings for the plot #
#########################

genotype_order_whiteflies = c("LA0716" ,  "PI127826", "LA1777" ,  "LYC4"   ,  "PI134418", "LA1718" ,  "LA1954" ,  "LA2695"  , "LA4024"  , "LA2172"  , "LA1401" , 
                              "LA0407" ,  "LA1578"  , "LA1364" ,  "LA2133" ,  "MM"   ,    "LA1840"  , "LA0735" ,  "LA1278")

genotype_order_thrips = c("LYC4","LA0407", "LA1777", "PI134418",
                          "LA1401", "LA0716", "LA1278",  "LA2172", 
                          "LA2695","LA0735","LA1718","LA2133","PI127826",
                          "LA1578", "LA1954", "MM",  "LA1840", "LA4024","LA1364")

# Set order of the accessions (from low to high survival)
volatiles.long.candidates.species$accession = factor(volatiles.long.candidates.species$accession, 
                                                     levels = genotype_order_whiteflies, 
                                                     ordered = TRUE)


# Custom theme for plot
my.theme = theme(axis.text.x = element_text(color = "black", size = 6, angle = 45, hjust = 1),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8),
                 strip.text.x = element_text(size = 8, colour = "black"),
                 legend.text = element_text(size = 8, colour = "black"),
                 legend.position = "none"
)

##########################
# Plot PLS-DA candidates #
##########################

p.volatiles = 
volatiles.long.candidates.species %>%
  dplyr::group_by(accession, tentative_id, species, color) %>% 
  summarise(mean_abundance = mean(log(abundance+1)), 
            n = n(), 
            se = (sd(log(abundance+1))/sqrt(n))
  ) %>% 
  ggplot(.) +
  geom_bar(aes(x = accession, y = mean_abundance,fill=species), stat = "identity",color="black") + 
  geom_errorbar(
    aes(x = accession, 
        ymin = mean_abundance - se, 
        ymax = mean_abundance + se,
        width = 0.4)
  ) +
  facet_wrap(~tentative_id, scale = "free", ncol = 2) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_x_discrete("accession", labels = genotype_order_whiteflies)+
  scale_colour_manual(values=volatiles.long.candidates.species$color)+
  theme_bw()+
  my.theme

ggsave(file = "PLS-DA script/volatiles_PLSDA_candidates_WF.pdf", plot = p.volatiles)
