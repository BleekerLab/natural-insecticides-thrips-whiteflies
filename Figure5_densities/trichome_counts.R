library("tidyverse")

###########################
# Data import and wrangling
###########################
# import data
df = read.delim("~/Documents/workspace/trichome_counts.tsv",header = T,stringsAsFactors = F)

# make it tidy
df = gather(df,key = "type",value = "counts",-Accession,-Plant,-Side,-Leaf_disc,-Name,-Date)

# remove date (not necessary)
df$Date = NULL

# for each genotype individual, two leaf discs were cut from one leaflet.
# Are the counts from the two leaf discs different?
leaf.disc.1 = df %>% filter(Leaf_disc == 1) %>% select(counts)  
leaf.disc.2 = df %>% filter(Leaf_disc == 2) %>% select(counts) 
shortest <- min(length(leaf.disc.1$counts), length(leaf.disc.2$counts))
leaf.disc.1 <- tail(leaf.disc.1$counts, shortest)
leaf.disc.2 <- tail(leaf.disc.2$counts, shortest)
t.test(leaf.disc.1,leaf.disc.2) # gives a p-value of 0.1156 so we can treat leaf discs as equivalent and average
df = df %>% group_by(Accession,Plant,Side,Name,type) %>% summarise(counts = mean(counts)) # average leaf discs

# Are scientists different?
p1 = df %>% filter(type == "All") %>% ggplot(.,aes(x=Name,y=counts,fill=Name)) + geom_boxplot() + facet_wrap(~ Accession)
res = TukeyHSD(aov(data = df,formula = counts ~ Name))
res = as.data.frame(res$Name);res$name = row.names(res);colnames(res)[4]="padj"
res = res[which(res$padj < 0.05),] # it seems that Michelle is different but only relatively to Sanne and Ruy. Let's keep her values.
df = ungroup(df) # to allow column removal
df$Name = NULL # no need to keep that column

# are leaf sides different?
leaf.side.1 = df %>% filter(Side == "abaxial") %>% select(counts)  
leaf.side.2 = df %>% filter(Side == "adaxial") %>% select(counts) 
t.test(leaf.side.1$counts,leaf.side.2$counts,paired = T,var.equal = F) # highly significant so they are different


###############
# Desired plots
###############
# plot 2 = which accession have the most trichomes ()
p2 = ggplot(df,aes(Accession,counts)) +
  geom_boxplot() +
  geom_point(aes(colour=Side)) + 
  facet_wrap(~ Side)
 
# add variables
#df = mutate(df,percent.glandular = round(Glandular / All * 100,digits = 1))
#df = mutate(df,disc.area = pi*2^2) 



