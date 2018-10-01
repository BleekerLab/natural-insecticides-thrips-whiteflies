# load libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(lme4)

# load .Rprofile file
source("./.custom_R_profile")

# where to save results?
resdir = file.path(getwd(),"Figure7_F2_whitefly_survival/")

# load dataset
df = read.delim("Figure7_F2_whitefly_survival/20170912_wf_bioassay_F2.txt",header = T,stringsAsFactors = F)

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
ggsave(filename = file.path(resdir,"Figure7A.png"),plot = g,width=10,height=5,dpi = 400)
ggsave(filename = file.path(resdir,"Figure7A.svg"),plot = g,width=10,height=5)

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

# Fit models
# dead / total - dead
fit1 = glm(cbind(dead,total-dead) ~ 1,data = df,family = binomial(link = logit)) # intercept only model
fit2 = glm(cbind(dead,total-dead) ~ 1 + line, data=df, family = binomial(link = logit)) # model with accession
fit3 = glmer(cbind(dead,total-dead) ~ 1 + line + (1|cage),data = df,family=binomial(link = logit)) # random plant effects

# AIC criteria and goodness of fit of the models
fits_no_mixed = list(fit1,fit2)
fit_mixed = list(fit3)

AICs = cbind.data.frame(
  models = sapply(c(fits_no_mixed,fits_mixed),function(x) {deparse(formula(x))}),
  AIC = sapply(c(fits_no_mixed,fits_mixed),function(x) AIC(x)),
  dAIC = as.vector(bbmle::AICtab(fit1,fit2,fit3)$dAIC)
)
AICs$models = c("model1","model2","model3")

plotAIC = ggplot(AICs,aes(models,AIC)) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  ggtitle("AIC for each model")
print(plotAIC)

