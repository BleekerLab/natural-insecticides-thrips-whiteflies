################
# load libraries
################
library(tidyverse)
require(gridExtra) # for side to side plotting


##################
# load dataframes
##################
volatiles.19acc = read.delim("Figure5_heatmaps/pheno_terpenoids.tsv",header=T,stringsAsFactors = F,check.names = F)
volatiles.F2 = read.delim("Figure8_volatiles_F2/20180625_F2metabolite_table.tsv",header = T,stringsAsFactors = F,check.names = F)

###############################################
# Count number of times a volatile was detected
###############################################
volatiles.19acc.subset = volatiles.19acc[,-c(2,3)] # remove phenotypes
counts = colSums(volatiles.19acc.subset != 0) # count the number of times a certain volatile is detected 
counts = as.data.frame(counts)
counts.19accessions = data.frame(name=row.names(counts),occurence=counts$counts)
counts.19accessions$percentage = counts.19accessions$occurence / 19 * 100
counts.19accessions = counts.19accessions[2:nrow(counts.19accessions),] # remove first line (sample line)


counts = colSums(volatiles.F2 != 0)
counts=as.data.frame(counts)
counts.F2 = data.frame(name=row.names(counts),occurence=counts$counts)
counts.F2$percentage = counts.F2$occurence / 25 * 100
counts.F2 = counts.F2[2:nrow(counts.F2),] # remove first line (line line)

# order dataframes
counts.19accessions = counts.19accessions[order(counts.19accessions$percentage,decreasing = T),]
counts.F2 = counts.F2[order(counts.F2$percentage,decreasing = T),]

# reorder factors
counts.19accessions$name = factor(counts.19accessions$name,levels = counts.19accessions$name)
counts.F2$name = factor(counts.F2$name,levels = counts.F2$name)

###########
# plots
##########
p1 = ggplot(data = counts.19accessions,aes(x = name,y=occurence)) +
  geom_bar(stat = "identity",width = 0.7) +
  theme(axis.text.x = element_text(angle = 40,hjust = 1,size = 6)) +
  labs(x = "Volatile compound",y="Number of compound occurences") +
  scale_y_continuous(limits=c(0,20)) +
  ggtitle("Occurences of volatiles in the selected 19 genotypes")
p1

p2 = ggplot(data = counts.F2,aes(x = name,y=occurence)) +
  geom_bar(stat = "identity",width = 0.7) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 8)) +
  labs(x = "Volatile compound",y="Number of compound occurences") +
  scale_y_continuous(limits=c(0,30)) +
  ggtitle("Occurences of volatiles in the 22 F2 lines, their parents and their F1 (n=25)")
p2

pdf("FigureS2_volatiles_19acc_vs_F2/occurences.pdf",width = 7,height = 10)
grid.arrange(p1,p2)
dev.off()
