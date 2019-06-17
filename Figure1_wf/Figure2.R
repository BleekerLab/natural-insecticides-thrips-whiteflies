# library
library(tidyverse)
library(lme4)

############ Data import and wrangling #########
# import whitefly no-choice data
df = read.delim("Figure2_wf/whitefly_no-choice_19_accessions.txt",header = T,stringsAsFactors = F)

# remove unecessary variables
df$alive = NULL
df$dead = NULL
df$total = NULL

# average clip-cages results
df = df %>% group_by(accession,plant) %>% summarise(average = mean(percentage,na.rm = T))

# import accession to species
accession2species = read.delim("Figure2_wf/accession2species.txt",header = T,stringsAsFactors=F)
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
ggsave(filename = file.path("Figure2_wf/","wf_survival_19_accessions.png"),plot = survival,width=10,height=5,dpi = 400)
ggsave(filename = file.path("Figure2_wf/","wf_survival_19_accessions.svg"),plot = survival,width=10,height=5)

############ Generalized Linear Model ##########
# remove missing values
# convert cage to factor
# Calculate a probability to die for each observation
df = read.delim("Figure2_wf/whitefly_no-choice_19_accessions.txt",header = T,stringsAsFactors = F)
df = na.omit(df)

# Remove LA0716 since it has only zero values (cannot estimate the coefficient)
df = dplyr::filter(.data = df,accession != "LA0716")

# Fit models
# dead / total - dead
# -1 means no intercept
fit1 = glm(cbind(dead,total-dead) ~ -1 + accession, data=df, family = binomial(link = logit)) # model with accession
fit2 = glmer(cbind(dead,total-dead) ~ 1 + accession + (1|accession/plant),data = df,family=binomial(link = logit)) # random fixed plant effects
fit3 = glmer(cbind(dead,total-dead) ~ 1 + accession + (1|accession/plant/cage),data = df,family=binomial(link = logit)) # random fixed plant/cage effects

# Comparing models
anova(fit1,fit2,test="LRT")
anova(fit2,fit3,test="LRT")

# extract coefficients
coefs4 = as.data.frame(summary(fit4)$coefficients)
colnames(coefs4)=c("coeff","stderr","zval","pval")
coefs4$accession = row.names(coefs4$accession)
coefs4$accession = sub("^accession",replacement = "",x = row.names(coefs))
coefs4$accession[1] = "Intercept"

# add species
coefs4 = dplyr::left_join(coefs4,accession2species)

# which conditions have the significant pvalues?
coefs4$text = NA
coefs4[which(coefs4$pval > 0.05),]$text <- "ns"
coefs4[which(coefs4$pval < 0.05),]$text <- "***"

