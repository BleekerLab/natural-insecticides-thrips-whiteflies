library("tidyverse")

###########################
# Data import and wrangling
###########################
# import data
df = read.delim("Figure5_densities/trichome_counts.tsv",header = T,stringsAsFactors = F)

# make it tidy
df = gather(df,key = "type",value = "counts",-genotype,-plant,-side,-leaf.disc,-name,-date)

# remove date (not necessary)
df$date = NULL

# for each genotype individual, two leaf discs were cut from one leaflet.
# Are the counts from the two leaf discs different?
leaf.disc.1 = df %>% filter(leaf.disc == 1) %>% select(counts)  
leaf.disc.2 = df %>% filter(leaf.disc == 2) %>% select(counts) 
shortest <- min(length(leaf.disc.1$counts), length(leaf.disc.2$counts))
leaf.disc.1 <- tail(leaf.disc.1$counts, shortest)
leaf.disc.2 <- tail(leaf.disc.2$counts, shortest)
t.test(leaf.disc.1,leaf.disc.2) # gives a p-value of 0.1156 so we can treat leaf discs as equivalent and average
df = df %>% group_by(genotype,plant,side,name,type) %>% summarise(counts = mean(counts)) # average leaf discs
df = ungroup(df)

# for each genotype individual, several plants were assayed
df$plant = factor(df$plant,levels = unique(df$plant))
res = aov(data=df,formula = counts ~ plant)
res = TukeyHSD(res) # no significant difference
df$plant = NULL

# Are scientists counting differently?
# res = TukeyHSD(aov(data = df,formula = counts ~ Name))
# res = as.data.frame(res$Name);res$name = row.names(res);colnames(res)[4]="padj"
# res = res[which(res$padj < 0.05),] # it seems that Michelle is different but only relatively to Sanne and Ruy. Let's remove her values.
# df = ungroup(df) # to allow column removal
# df = df %>% filter(Name != "Michelle")
# df$Name = NULL # no need to keep that column

p = df %>% 
  filter(type == "TypeVI") %>% 
  ggplot(.,aes(x=Name,y=counts,fill=Name)) + 
  geom_boxplot() + 
  facet_wrap(~ Accession) 
# this plots shows:
#   - that LA1578 is underestimated by Marc compared to Michelle and Sanne's results
#   - that LA1718 is overestimated by Michelle compared to Marc, Ruy and Sanne

# are leaf sides different?
leaf.side.1 = df %>% filter(side == "abaxial") %>% select(counts)  
leaf.side.2 = df %>% filter(side == "adaxial") %>% select(counts) 
t.test(leaf.side.1$counts,leaf.side.2$counts,paired = T,var.equal = F) # highly significant so they are different


###############
# Desired plots
###############

##### adaxial side
p1 = df %>% 
  filter(side == "adaxial") %>%
  ggplot(.,aes(genotype,counts,fill=genotype)) +
  geom_boxplot() +
  facet_wrap(~ type,scales = "free") +
  theme(axis.text.x = element_text(angle=90)) 
p1

nMeasures = filter(df,side =="adaxial") %>% with(data = .,table(genotype,type)) %>% as.data.frame()
# add variables
#df = mutate(df,percent.glandular = round(Glandular / All * 100,digits = 1))
#df = mutate(df,disc.area = pi*2^2) 

anova.adaxial = df %>% filter(side == "adaxial") %>% aov(data = .,formula = counts ~ type + Accession)
