library("tidyverse")
library("ggpubr") # publication-ready plots


###########################
# Data import and wrangling
###########################
# import data
df = read.delim("Figure5_densities/trichome_counts.tsv",header = T,stringsAsFactors = F)
species.info = read.delim("genotype2species.txt",header = T,stringsAsFactors = F)

# make it tidy
df = gather(df,key = "type",value = "counts",-genotype,-plant,-leaf.side,-leaf.disc,-name,-date)

# remove date (not necessary)
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
df = df %>% group_by(genotype,leaf.side,name,type) %>% summarise(counts = mean(counts)) # average leaf discs
df = ungroup(df)


###########
# Leaf side
###########
df$type = factor(df$type,levels = unique(df$type))

# plots
p.leafside = df %>%
  ggplot(aes(x = leaf.side,y=counts,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  ggtitle("Trichome counts on each leaf side (all trichome types)") +
  mytheme() +
  stat_compare_means(method = "t.test",label.x = 1.5)


p.leafside.per.type = df %>%
  ggplot(aes(x = leaf.side,y=counts,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  ggtitle("Trichome counts on each leaf side per trichome type") +
  facet_wrap(~ type,scales = "free") +
  mytheme() +
  stat_compare_means(method = "t.test",label.x.npc = "center",label.y.npc = "top")


p.leafside.per.genotype = df %>%
  ggplot(aes(x = leaf.side,y=counts,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  facet_wrap(~ genotype,scales = "free") +
  mytheme() +
  ggtitle("Trichome counts per genotype and per leaf side")  +
  stat_compare_means(method="t.test",label.x.npc = "center",label.y.npc =1)


p.adaxial = df %>%
  filter(leaf.side == "adaxial") %>%
  ggplot(aes(x = type,y=counts,fill=type)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.2) +
  mytheme() +
  ggtitle("Trichome counts on the adaxial leaf side")

ggsave(filename = "Figure5_densities/plots/leafside.png",plot = p.leafside,width = 7,height = 5)
ggsave(filename = "Figure5_densities/plots/leafside.per.type.png",plot = p.leafside.per.type,width = 7,height = 5)
ggsave(filename = "Figure5_densities/plots/leafside.per.genotype.png",plot = p.leafside.per.genotype,width = 10,height = 6)

ggsave(filename = "Figure5_densities/plots/leafside.svg",plot = p.leafside,width = 7,height = 5)
ggsave(filename = "Figure5_densities/plots/leafside.per.type.svg",plot = p.leafside.per.type,width = 7,height = 5)
ggsave(filename = "Figure5_densities/plots/leafside.per.genotype.svg",plot = p.leafside.per.genotype,width = 7,height = 5)


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


###############
# Desired plots
###############
types_to_keep = c("non.glandular","typeIandIV","typeVI","typeVII")



#df = df[which(df$type %in% types_to_keep),] %>%
#  ggplot(.,aes(x=genotype,y=counts,fill=genotype)) + 
# geom_violin() + 
#  geom_jitter(width = 0.2) +
#  theme(axis.text.x = element_text(angle=90)) 
  



