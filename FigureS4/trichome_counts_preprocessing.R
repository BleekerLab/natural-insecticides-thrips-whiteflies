library("tidyverse")
library("ggpubr")
library("gridExtra")
library("agricolae")

###########################
# Data import and wrangling
###########################
# import data
# make it tidy
# remove date (not necessary)
df = read.delim("Figure5_densities/trichome_counts.tsv",header = T,stringsAsFactors = F)
df = gather(df,key = "type",value = "counts",-genotype,-accession,-plant,-leaf.side,-leaf.disc,-name,-date)
df$date = NULL

#######################
# custom theme for plot 
#######################
mytheme <- function(){
  theme_grey(base_size = 12) +
    theme(
      text = element_text(colour="black"),
      plot.title = element_text(size=18)
    )
}

#################################################
# Can some factors be dropped if not significant?
#################################################


#######
# Plant
#######
# for each genotype individual, several plants were assayed
# are individual plants different?
df$plant = factor(df$plant,levels = unique(df$plant))
res = aov(data=df,formula = counts ~ plant)
res = TukeyHSD(res) # no significant difference: there is no difference between plant individuals. 
df$plant = NULL     # the plant factor can be thus removed

###########
# Leaf disc
###########
# for each genotype individual, two leaf discs were cut from one leaflet.
# Are the counts from the two leaf discs different?
leaf.disc.1 = df %>% filter(leaf.disc == 1) %>% select(counts)  
leaf.disc.2 = df %>% filter(leaf.disc == 2) %>% select(counts) 
shortest <- min(length(leaf.disc.1$counts), length(leaf.disc.2$counts))
leaf.disc.1 <- tail(leaf.disc.1$counts, shortest)
leaf.disc.2 <- tail(leaf.disc.2$counts, shortest)
t.test(leaf.disc.1,leaf.disc.2) # p-value > 0.05 so we can treat leaf discs as equivalent and average
df = df %>% group_by(genotype,accession,leaf.side,name,type) %>% summarise(counts = mean(counts)) # average leaf discs
df = ungroup(df)


###########
# Leaf side
###########
df$type = factor(df$type,levels = unique(df$type))

# is there a difference between leaf side depending on the type of trichomes? 
types.of.trichomes = as.vector(unique(df$type))
models <- sapply(types.of.trichomes, function(x) {
  lm(counts ~ leaf.side, data=df, type==x)
}, simplify=FALSE)
anova.tables.type = sapply(models, anova, simplify=FALSE)       # One-way ANOVA for each trichome type    
anova.tables.type = plyr::ldply(anova.tables.type)              # converts list into one dataframe
colnames(anova.tables.type)[6] = "pval"
type.aov = filter(anova.tables.type,Df == 1 & pval < 0.05)     # keeps only leaf side factor (degree freedom == 1)


##############################
# type 4 counts = f(scientist)
##############################
# Are scientists counting differently type 6 trichomes?
hsd.type14 = df %>% 
  filter(type == "typeIandIV") %>% 
  filter(counts > 0) %>%
  aov(formula = counts ~ name) %>% 
  TukeyHSD() 
hsd.type14 = as.data.frame(hsd.type14$name)
hsd.type14$name = row.names(hsd.type14)
colnames(hsd.type14)[4]="padj"
hsd.type14[which(hsd.type14$padj < 0.05),] # only significant difference (padj = 0.04420046 < 0.05) is between Michelle & Marc. But is marginally significant so we keep all data from all scientists.

p.type4_counts_per_scientist = df %>% 
  filter(type == "typeIandIV") %>% 
  ggplot(.,aes(x=name,y=counts,fill=name)) + 
  geom_boxplot() + 
  facet_wrap(~ genotype) +
  theme(axis.text.x = element_text(angle=90)) + 
  ggtitle("type 4 trichome counts per scientist & genotype")


##############################
# type 6 counts = f(scientist)
##############################
# Are scientists counting differently type 6 trichomes?
hsd.type6 = df %>% 
  filter(type == "typeVI") %>% 
  filter(counts > 0) %>%
  aov(formula = counts ~ name) %>% 
  TukeyHSD() 
hsd.type6 = as.data.frame(hsd.type6$name)
hsd.type6$name = row.names(hsd.type6)
colnames(hsd.type6)[4]="padj"
dim(hsd.type6[which(hsd.type6$padj < 0.05),]) # no signif diff between scientists for type 6 (we keep all measurements)

p.type6_counts_per_scientist_and_genotype = df %>% 
  filter(type == "typeVI") %>% 
  ggplot(.,aes(x=name,y=counts,fill=name)) + 
  geom_boxplot() + 
  facet_wrap(~ genotype) +
  theme(axis.text.x = element_text(angle=90)) + 
  ggtitle("type 6 trichome counts per scientist & genotype")

p.type6_counts_per_scientist = df %>% 
  filter(type == "typeVI") %>% 
  ggplot(.,aes(x=name,y=counts,fill=name)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  theme(axis.text.x = element_text(angle=90)) + 
  ggtitle("type 6 trichome counts per scientist & genotype")

p.ng_counts_per_scientist = df %>% 
  filter(type == "non.glandular") %>% 
  ggplot(.,aes(x=name,y=counts,fill=name)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  theme(axis.text.x = element_text(angle=90)) + 
  ggtitle("non-glandular trichome counts per scientist & genotype")


###############
# Desired plots
###############
# add density column
leaf.disc.radius = 2 # mm2
leaf.disc.area = pi * leaf.disc.radius^2
df = mutate(df,density = counts / leaf.disc.area)

# add species and corresponding color for species
species.info = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)
df = left_join(df,species.info[,c("species","color","accession")],by="accession")

# what trichomes?
types_to_keep = c("non.glandular","typeIandIV","typeVI")
df = filter(.data = df,type %in% types_to_keep) # only keep the trichomes of interest 

# write table
write.table(x = df,file = "Figure5_densities/trichome.counts.processed.tsv",quote = F,sep = "\t",row.names = F)


###########
# Figure 5A
############
p.leafside = df %>%
  ggplot(aes(x = leaf.side,y=density,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  ggtitle("Trichome density (trichomes/mm2) on each leaf side (all trichome types)") +
  mytheme() +
  stat_compare_means(method = "t.test",label.x = 1.5)


p.leafside.per.type = df %>%
  ggplot(aes(x = leaf.side,y=density,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  ggtitle("Trichome density (trichomes/mm2) on each leaf side per trichome type") +
  facet_wrap(~ type,scales = "free") +
  mytheme() +
  stat_compare_means(method = "t.test",label.x.npc = "center",label.y.npc = "top")


p.leafside.per.genotype = df %>%
  ggplot(aes(x = leaf.side,y=density,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  facet_wrap(~ genotype,scales = "free") +
  mytheme() +
  ggtitle("Trichome density (trichomes/mm2) per genotype and per leaf side")  +
  stat_compare_means(method="t.test",label.x.npc = "center",label.y.npc =1)

#ggsave(filename = "Figure5_densities/plots/leafside.pdf",plot = p.leafside,width = 7,height = 5)
ggsave(filename = "Figure5_densities/plots/Figure5A.pdf",plot = p.leafside.per.type,width = 10,height = 7)
#ggsave(filename = "Figure5_densities/plots/leafside.per.genotype.pdf",plot = p.leafside.per.genotype,width = 20,height = 10)

#ggsave(filename = "Figure5_densities/plots/leafside.svg",plot = p.leafside,width = 7,height = 5)
ggsave(filename = "Figure5_densities/plots/Figure5A..svg",plot = p.leafside.per.type,width = 7,height = 5)
#ggsave(filename = "Figure5_densities/plots/leafside.per.genotype.svg",plot = p.leafside.per.genotype,width = 7,height = 5)


##########################
# ANOVA and post-hoc tests
##########################
# helper function
add_row_names <- function(df){
  data.table::setDT(df,keep.rownames = "genotype")[]
}

### adaxial (upper panel)
adaxial = df %>% filter(leaf.side == "adaxial")
types.of.trichomes = as.vector(unique(adaxial$type))
models.adaxial <- sapply(types.of.trichomes, function(x) {
  lm(density ~ genotype, data=filter(adaxial,type == x))
}, simplify=FALSE)
models.adaxial = sapply(models.adaxial, aov, simplify=FALSE)       # One-way ANOVA for each trichome type    
models.adaxial = lapply(X = models.adaxial,FUN = function(x){HSD.test(x,trt = "genotype",alpha = 0.05,group = TRUE)$groups})

models.adaxial$non.glandular$genotype = row.names(models.adaxial$non.glandular)
models.adaxial$typeIandIV$genotype = row.names(models.adaxial$typeIandIV)
models.adaxial$typeVI$genotype = row.names(models.adaxial$typeVI)

models.adaxial$non.glandular$type = "non.glandular"
models.adaxial$typeIandIV$type = "typeIandIV"
models.adaxial$typeVI$type = "typeVI"


### abaxial (lower panel)
abaxial = df %>% filter(leaf.side == "abaxial")
models.abaxial <- sapply(types.of.trichomes, function(x) {
  lm(density ~ genotype, data=filter(abaxial,type == x))
}, simplify=FALSE)
models.abaxial = sapply(models.abaxial, aov, simplify=FALSE)       # One-way ANOVA for each trichome type    
models.abaxial = lapply(X = models.abaxial,FUN = function(x){HSD.test(x,trt = "genotype",alpha = 0.05,group = TRUE)$groups})

models.abaxial$non.glandular$genotype = row.names(models.abaxial$non.glandular)
models.abaxial$typeIandIV$genotype = row.names(models.abaxial$typeIandIV)
models.abaxial$typeVI$genotype = row.names(models.abaxial$typeVI)

models.abaxial$non.glandular$type = "non.glandular"
models.abaxial$typeIandIV$type = "typeIandIV"
models.abaxial$typeVI$type = "typeVI"

# combine HSD group results into two dataframes
adaxial = with(data = models.adaxial,rbind.data.frame(non.glandular,typeIandIV,typeVI))
abaxial = with(data = models.abaxial,rbind.data.frame(non.glandular,typeIandIV,typeVI))

# merge into one unique dataframe containing both leaf side info
adaxial$leaf.side = "adaxial"
abaxial$leaf.side = "abaxial"


adaxial
test = left_join(adaxial,species.info,by="genotype")


#adaxial_and_abaxial = rbind.data.frame(adaxial,abaxial)                              

# write to file
write.table(x = adaxial_and_abaxial,file = "Figure5_densities/hsd.groups.tsv",sep = "\t",row.names = F,quote = F)

# add info to main dataframe
#df = left_join(df,adaxial_and_abaxial,by = c(
#  "genotype" = "genotype",
# "type" = "type",
#  "leaf.side" = "leaf.side"
#))



############
# Figure 5B
############
p.adaxial = df %>% 
  filter(leaf.side == "adaxial") %>%
  ggplot(.,aes(x=genotype,y=density,fill=species)) + 
  geom_boxplot() + 
  geom_point() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5, hjust=1))  +
  facet_wrap(~ type,scales = "free") +
  labs(y = "Adaxial trichome density (trichomes/mm2)")
p.adaxial + geom_text(data = adaxial,label="groups")

p.abaxial = df %>% 
  filter(leaf.side == "abaxial") %>%
  ggplot(.,aes(x=genotype,y=density,fill=species)) + 
  geom_boxplot() + 
  geom_point() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5, hjust=1))  +
  facet_wrap(~ type,scales = "free") +
  labs(y = "Abaxial trichome density (trichomes/mm2)")

svg("Figure5_densities/plots/Figure5B.svg",width = 10,height = 7)
grid.arrange(p.adaxial,p.abaxial,nrow=2)
dev.off()

pdf("Figure5_densities/plots/Figure5B.pdf",width = 10,height = 7)
grid.arrange(p.adaxial,p.abaxial,nrow=2)
dev.off()

ggsave(filename = "Figure5_densities/plots/Figure5B.adaxial.pdf",plot = p.adaxial,width = 10,height = 5)
ggsave(filename = "Figure5_densities/plots/Figure5B.abaxial.pdf",plot = p.abaxial,width = 10,height = 5)

#################
# session info
##################
writeLines(capture.output(sessionInfo()), "Figure5_densities/sessionInfo.txt")
