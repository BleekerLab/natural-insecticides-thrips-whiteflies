# load libraries
library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

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
  geom_jitter(width = 0,height = 0) +
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
