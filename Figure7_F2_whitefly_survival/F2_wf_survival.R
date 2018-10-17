# load libraries
library(reshape2)
library(ggplot2)
suppressMessages(library(dplyr))
library(RColorBrewer)
library(lme4)

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

########### Fit models #########
# -1 means no intercept
# Generalized linear model for binomial distribution (dead / total - dead)
# We first remove lines that have only 0s because they impair coefficient estimation 
# Lines 294, 324, 373, 387, 396 have no surviving whiteflies
zero_lines = df %>% group_by(line) %>% summarise(alive = sum(alive)) %>% filter(alive ==0) %>% select(line) %>% unlist(line)
df.minus.lines.with.only.zeros = df[! df$line %in% zero_lines,]

# We then fit a Generalized Linear Model (line is the only parameter with an influence. Cage number does not)
fit2 = glm(cbind(dead,total-dead)~-1 + line,data=df.minus.lines.with.only.zeros,family = binomial(link=logit))


