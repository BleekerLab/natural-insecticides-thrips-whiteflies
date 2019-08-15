# library
suppressMessages(library(tidyverse))

############ Data import and wrangling #########
# import whitefly no-choice data
df = read.delim("Figure1_wf/whitefly_no-choice_19_accessions.tsv",header = T,stringsAsFactors = F)

# remove unecessary variables
df$alive = NULL
df$dead = NULL
df$total = NULL

# average clip-cages results
df = df %>% group_by(accession,plant) %>% summarise(average = mean(percentage,na.rm = T))

# import accession to species
accession2species = read.delim("genotype2species.txt",header = T,stringsAsFactors=F)
df = dplyr::left_join(x = df,y = accession2species)

# extract ordering by increasing median survival values
newOrder = group_by(df,accession) %>%
  summarise(median = median(average)) %>%
  arrange(.,median) %>%
  select(accession)
newOrder = newOrder$accession


########## Plot ################
# reorder factor levels for accession in melted dataframe
df$accession = factor(df$accession,levels = newOrder)

# make the plot
survival <- ggplot(data = df,aes(x = accession,y = average,fill=species)) +
  geom_boxplot() +
  stat_summary(fun.y="mean",geom="point",shape=23,size=2,fill="white") +
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1)) + 
  ggtitle("Whitefly survival after five days") +
  labs(x="Tomato genotype",y="Whitefly survival (%)") +
  theme(plot.title = element_text(size=18),axis.text = element_text(colour = "black"))

# print the plot in the final document
print(survival)

# save plot into a file
ggsave(filename = file.path("Figure1_wf/","wf_survival_19_accessions.png"),plot = survival,width=10,height=5,dpi = 400)
ggsave(filename = file.path("Figure1_wf/","wf_survival_19_accessions.svg"),plot = survival,width=10,height=5)

################################################