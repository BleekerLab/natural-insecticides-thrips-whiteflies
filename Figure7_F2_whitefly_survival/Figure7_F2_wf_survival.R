# load libraries
library(reshape2)
library(ggplot2)
suppressMessages(library(dplyr))
library(RColorBrewer)
library(lme4)
library(tibble)

# load dataset
df = read.delim("Figure7_F2_whitefly_survival/20170912_wf_bioassay_F2.tsv",header = T,stringsAsFactors = F)

######## Plot of whitefly survival per plant line #######
# percentage per plant line
df = df %>% 
  mutate(percentage = alive / (alive + dead) * 100) 

# calculate median survival per plant (to order boxplots)
plant_ranking = df %>%
  group_by(line) %>%
  summarise(median=median(percentage)) %>% arrange(median) %>% select(line) 
plant_ranking = plant_ranking$line

# create a gradient of colors
# colorRampPalette returns a function that takes an integer argument (the required number of colors) and returns a character vector of colors
color_func <- colorRampPalette(colors = c("darkorange1","mediumpurple"))
colorPalette <- color_func(length(plant_ranking))

# plot
df$line = factor(df$line,levels = plant_ranking)

g <- ggplot(data = df,aes(x=line,y=percentage,fill=line)) + 
  geom_boxplot() +
  geom_point() +
  geom_jitter(width = 0,height = 1) +
  guides(fill=FALSE) +
  scale_fill_manual(values = colorPalette) +
  scale_y_continuous(limits = c(0,100)) +
  labs(x="Plant line id",y="Whitefly survival (%)") +
  theme(axis.text = element_text(color="black"),
        axis.text.x = element_text(angle = 90))

print(g)

# save plot
ggsave(filename = "Figure7_F2_whitefly_survival/Figure7A.png",plot = g,width=10,height=5,dpi = 400)
ggsave(filename = "Figure7_F2_whitefly_survival/Figure7A.svg",plot = g,width=10,height=5)

######### Logistic regressions (Generalized Linear Models and Generalized Linear Mixed Models) ##########
# remove missing values
# convert cage to factor
# Calculate a probability to die for each observation
df = na.omit(df)
df$cage = as.factor(df$cage)
df$line = as.factor(df$line)

# add total column
df = df %>% 
  mutate(total=dead+alive)

########### Fit models #########
# -1 means no intercept
# Generalized linear model for binomial distribution (dead / total - dead)
# We first remove lines that have only 0s because they impair coefficient estimation 
# Lines 294, 324, 373, 387, 396 have no surviving whiteflies
zero_lines = df %>% group_by(line) %>% summarise(alive = sum(alive)) %>% filter(alive ==0) %>% select(line) %>% unlist(line)
df.minus.lines.with.only.zeros = df[! df$line %in% zero_lines,]

# We then fit a Generalized Linear Model (line is the only parameter with an influence. Cage number does not)
fit = glm(cbind(dead,total-dead)~-1 + line,data=df.minus.lines.with.only.zeros,family = binomial(link=logit))

# extract coefficients from fitted model
stats.model = as.data.frame(coefficients(summary(fit)))
lines = gsub("line","",x = row.names(stats.model)) # extract line numbers
stats.model = tibble::add_column(stats.model,lines,.before = "Estimate")
write.table(x = stats.model,file = "Figure7_F2_whitefly_survival/stats.model.tsv",quote = F,row.names = F,sep = "\t")

######### Extracts F2 lines and make a table with "toxic" and "non-toxic" effects  
colnames(stats.model) = c("line","coef","stderr","zvalue","pvalue")

# toxic lines = significant coefs & positively contribute (coef>0) to death probability 
toxic = stats.model %>% 
  filter(pvalue < 0.05 & coef > 0) 

# adding the lines with the "supertoxic zero survivors" phenotype
# these lines were not included in the GLM regression model but are toxic (no surviving whiteflies!)
supertoxic = data.frame(line = as.character(zero_lines),coef="NA",stderr="NA",zvalue="NA",pvalue="NA",stringsAsFactors = F)

# combine toxic with the "supertoxic zero survivors" lines
supertoxic_with_toxic = rbind(supertoxic,toxic,stringsAsFactors = F) 
write.table(supertoxic_with_toxic,file = "Figure7_F2_whitefly_survival/toxic_lines.tsv",sep = "\t",quote = F,row.names = F)

# non-toxic lines = all other lines with a coefficient < 0 
nontoxic = stats.model[which(! stats.model$line %in% toxic$line),]
nontoxic = nontoxic %>% filter(coef < 0)
write.table(nontoxic,file = "Figure7_F2_whitefly_survival/nontoxic_lines.tsv",sep = "\t",quote = F,row.names = F)

### Session info
writeLines(capture.output(sessionInfo()), "Figure7_F2_whitefly_survival/sessionInfo.txt")
