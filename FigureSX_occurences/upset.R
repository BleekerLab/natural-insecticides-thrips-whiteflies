library("plyr")
library("tidyverse")
library("RColorBrewer")
library("ggpubr")
library("superheat")

#############
# Functions
#############
make_classes_of_occurences <- function(matrix_of_metabolites){
  # takes a matrix of metabolites (first column contains first metabolite)
  # turns in into a binary matrix and count occurences + make classes of occurences
  #           met1             met2        met3          met4      met5
  #1          0.00             0           0.0          0.00      88210.77
  #2     234383.79             0           0.0          0.00     131824.19
  #3          0.00             0      104528.8      45335.47      27304.45
  #4          0.00             0           0.0          0.00     457891.16

  # binary matrix
  matrix_of_metabolites[matrix_of_metabolites > 0] <- 1 # all detected metabolites get a 1
  # count number of times each compound occurs
  occurences = data.frame(freq=colSums(matrix_of_metabolites))
  classes.occurences = as.data.frame(table(occurences$freq))
  colnames(classes.occurences)=c("class.occurence","freq")
  # number of metabolites
  nTotal = sum(classes.occurences$freq)
  # calculate percentages per class of occurences
  classes.occurences = classes.occurences %>% 
    mutate(percentage = round(freq / nTotal * 100,digits = 1)) %>% 
    mutate(lab.ypos = cumsum(percentage) - 0.5*percentage)
  return(classes.occurences)
}

##########
# Terpenes
##########
terpenes <- read.delim("Figure6_heatmaps/pheno_terpenoids.tsv",header = T,stringsAsFactors = F,check.names = F)
terpenes = terpenes[,4:ncol(terpenes)] # keep only compound abundances
terpenes[terpenes > 0] <- 1

# make classes of occurences
classes.occurences.terpenes = make_classes_of_occurences(terpenes)

############
# Acylsugars
############
acylsugars <- read.delim("Figure6_heatmaps/pheno_acylsugars.tsv",header = T,stringsAsFactors = F,check.names = F)
acylsugars = acylsugars[,4:ncol(acylsugars)]
acylsugars[acylsugars > 0] <- 1

# make classes of occurences
classes.occurences.acylsugars = make_classes_of_occurences(acylsugars)

#######
# Plots
#######
# bar plot of volatile occurences
p.bars.volatiles = ggplot(data = classes.occurences.terpenes,aes(x = class.occurence,y=percentage,label=percentage)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 0)) +
  ylim(c(0,40)) +
  labs(x = "Number of times a volatile compound is occuring in the N = 19 tomato genotype collection",y = "Percentage of occurence (%, 100% corresponds to N=19 genotypes)") +
  geom_label(size=3) +
  theme(axis.text = element_text(color = "black"))

p.bars.as = ggplot(data = classes.occurences.acylsugars,aes(x = class.occurence,y=percentage,label=percentage)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 0)) +
  ylim(c(0,60)) +
  labs(x = "Number of times a volatile compound is occuring in the N = 19 tomato genotype collection",y = "Percentage of occurence (%, 100% corresponds to N=19 genotypes)") +
  geom_label(size=3) + 
  theme(axis.text = element_text(color = "black"))


#################
# Save plots
#################
ggsave(filename = "FigureSX_occurences/FigureS2A.volatiles.barplot.pdf",plot = p.bars.volatiles,width = 10,height = 5)
ggsave(filename = "FigureSX_occurences/FigureS2B.acylsugars.barplot.pdf",plot = p.bars.as,width = 10,height = 5)

##################################
# Make heatmaps of sparse matrices
##################################
terpenes <- read.delim("Figure6_heatmaps/pheno_terpenoids.tsv",header = T,stringsAsFactors = F,check.names = F)
row.names(terpenes) = terpenes$sample
terpenes = terpenes[,4:ncol(terpenes)] # keep only compound abundances
terpenes[terpenes > 0] <- 1

png("FigureSX_occurences/FigS2C.volatiles.heatmap.png",height = 900,width=800)
superheat(terpenes,bottom.label = "none")
dev.off()

acylsugars <- read.delim("Figure6_heatmaps/pheno_acylsugars.tsv",header = T,stringsAsFactors = F,check.names = F)
row.names(acylsugars)=acylsugars$sample
acylsugars = acylsugars[,4:ncol(acylsugars)]
acylsugars[acylsugars > 0] <- 1

png("FigureSX_occurences/FigS2D.acylsugars.heatmap.png",height = 900,width=800)
superheat(acylsugars,bottom.label = "none")
dev.off()

#################
# session info
##################
writeLines(capture.output(sessionInfo()), "FigureS2_volatiles_19acc_vs_F2/sessionInfo.txt")