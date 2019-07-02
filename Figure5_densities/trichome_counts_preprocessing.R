library("tidyverse")

###########################
# Data import and wrangling
###########################
# import data
df = read.delim("Figure5_densities/trichome_counts.tsv",header = T,stringsAsFactors = F)

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
      plot.title = element_text(size=20)
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

# plots
p.leafside = df %>%
  ggplot(aes(x = leaf.side,y=counts,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.1) +
  ggtitle("Trichome counts on each leaf side (all trichome types)") +
  mytheme()

p.leafside.per.type = df %>%
  ggplot(aes(x = leaf.side,y=counts,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.1) +
  ggtitle("Trichome counts on each leaf side per trichome type") +
  facet_wrap(~ type,scales = "free") +
  mytheme()

p.leafside.per.genotype = df %>%
  ggplot(aes(x = leaf.side,y=counts,fill=leaf.side)) + 
  geom_violin() + 
  geom_jitter(width = 0.2) +
  facet_wrap(~ genotype,scales = "free") +
  mytheme()

p.leafside
p.leafside.per.type
p.leafside.per.genotype

# is there a difference between leaf side depending on the type of trichomes? 
types.of.trichomes = unique(df$type)
models <- sapply(types.of.trichomes, function(type) {
  lm(counts ~ leaf.side, data=df, type==type)
}, simplify=FALSE)
ANOVA.tables <- sapply(models, anova, simplify=FALSE) # they are all significant: leaf sides 


# are leaf sides different?
leaf.side.1 = df %>% filter(leaf.side == "abaxial") %>% select(counts)  
leaf.side.2 = df %>% filter(leaf.side == "adaxial") %>% select(counts) 
t.test(leaf.side.1$counts,leaf.side.2$counts,paired = T,var.equal = F) # highly significant so they are different



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
types_to_keep = c("typeIandIV","typeVI","non.glandular")




##### adaxial side 
p1 = df %>% 
  filter(side == "adaxial") %>%
  ggplot(.,aes(genotype,counts,fill=genotype)) +
  geom_boxplot() +
  facet_wrap(~ type,scales = "free") +
  theme(axis.text.x = element_text(angle=90)) 


nMeasures = filter(df,side =="adaxial") %>% with(data = .,table(genotype,type)) %>% as.data.frame()
# add variables
#df = mutate(df,percent.glandular = round(Glandular / All * 100,digits = 1))
#df = mutate(df,disc.area = pi*2^2) 

anova.adaxial = df %>% filter(side == "adaxial") %>% aov(data = .,formula = counts ~ type + genotype)
