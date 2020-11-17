library(mixOmics)
library(mdatools)
library(tidyverse)

# Load data
df = read.delim(file = "Random_Forest/phenotypes_vs_acylsugars.tsv", check.names = FALSE)

# Define phenotype of the accessions
pheno <- df$thrips

# Create the numeric matrix for sPLSDA
matrix <- df %>% select(-thrips, -wf ) %>% column_to_rownames(var = "sample")
matrix = log(matrix+1)

# Check the performace of the model 
plsda.results <- splsda(matrix, pheno, ncomp = 8)
set.seed(30)
plsda.performance <- perf(plsda.results, validation = "Mfold", folds = 3, 
     progressBar = FALSE, nrepeat = 50)
plot(plsda.performance, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")


################
# Do the PLSDA #
################

# Perform the splsda using the number of components giving the lowest error rate
plsda.results <- splsda(matrix, pheno, ncomp = 2)

# Plot samples
plotIndiv(plsda.results, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")

# plot variables
plotVar(plsda.results, cutoff = 0.8)

# Select candidates: the top 5 metabolites contributing to the phenotype
scores <- selectVar(plsda.results, comp=1)$value
candidates <- rownames(scores)[1:5]
candidates <- str_replace_all(candidates, "\\.", ":")
candidates <- str_replace_all(candidates, "C", "")
candidates <- str_replace_all(candidates, "_", "-")


#######################
# Plot the candidates #
#######################


########################################
# load original acylsugar measurements #
########################################
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
acylsugars = read.csv("Figure_4/20190904_acylsugars_peak_area_all_samples.csv", header = T, stringsAsFactors = TRUE, check.names = F)
acylsugars.long = gather(acylsugars, 
                         key = "metabolite",
                         value = "abundance",
                         -sample, -accession)


acylsugar.long.candidates = acylsugars.long %>% filter(metabolite %in% candidates)


### Read and add species and color information
acylsugar.candidates.with.species = left_join(acylsugar.long.candidates,accession2species,by="accession")


genotype_order_whiteflies = c("LA0716" ,  "PI127826", "LA1777" ,  "LYC4"   ,  "PI134418", "LA1718" ,  "LA1954" ,  "LA2695"  , "LA4024"  , "LA2172"  , "LA1401" , 
                              "LA0407" ,  "LA1578"  , "LA1364" ,  "LA2133" ,  "MM"   ,    "LA1840"  , "LA0735" ,  "LA1278")

genotype_order_thrips = c("LYC4","LA0407", "LA1777", "PI134418",
                          "LA1401", "LA0716", "LA1278",  "LA2172", 
                          "LA2695","LA0735","LA1718","LA2133","PI127826",
                          "LA1578", "LA1954", "MM",  "LA1840", "LA4024","LA1364")



my.theme = theme(axis.text.x = element_text(color = "black", size = 6, angle = 45, hjust = 1),
                 axis.text.y = element_text(color = "black", size = 6),
                 axis.title.x = element_text(color = "black", size = 8),
                 axis.title.y = element_text(color = "black", size = 8),
                 strip.text.x = element_text(size = 8, colour = "black"),
                 legend.text = element_text(size = 8, colour = "black"),
                 legend.position = "none"
)
########
# Plot #
########

# Order the accessions based on the phenotype
acylsugar.candidates.with.species$accession = factor(acylsugar.candidates.with.species$accession, 
                                                     levels = genotype_order_whiteflies, 
                                                     ordered = TRUE)

p.acylsugars = 
acylsugar.candidates.with.species %>%
  dplyr::group_by(accession, metabolite, species, color) %>% 
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
  facet_wrap(~metabolite, scale = "free", ncol = 2) +
  labs(x = "Tomato genotype", y="Mean normalised peak area (AU)") +
  scale_x_discrete("accession", labels = genotype_order_whiteflies)+
  scale_colour_manual(values=acylsugar.candidates.with.species$color)+
  theme_bw()+
  my.theme

ggsave(file = "PLS-DA script/Acylsugar_PLSDA_candidates_WF.pdf", plot = p.acylsugars)
