library("tidyverse")


#################
# Acylsugars data
#################
as = read.delim("Figure6_heatmaps/pheno_acylsugars.tsv",header = T,stringsAsFactors = F)
as = select(.data = as,-c(wf,thrips))

# sum over rows (sum of acylsugar peak areas per genotype)
sampleRowSums = rowSums(as[,2:ncol(as)])
as.summed = data.frame(genotype = as$sample,as.summed=sampleRowSums)

# compute log10 of summed peak area
as.summed = mutate(as.summed,log.as.summed = log10(as.summed))

################
# Trichomes data
################
trichomes = read.delim("FigureS4_spearman_trichome_vs_acylsugars/trichome.counts.processed.tsv",header = T,stringsAsFactors = F)

trichomes = trichomes %>% 
  filter(type == "typeIandIV") %>%  # keep only trichomes relevant for acylsugar synthesis
  dplyr::group_by(genotype,species,color,type) %>% 
  dplyr::summarise(density = mean(density)) 

# transform 0 into 1 for log scale plotting
trichomes$density[trichomes$density == 0] <- 1

trichomes = dplyr::mutate(.data = trichomes,log.density = log10(density))

###################################
# Spearman correlation (all values)
###################################
# join the two datasets
df = dplyr::inner_join(trichomes,as.summed,by="genotype")

# calculate R squared
r2 = with(data = df,
          round(
            cor(density,log.as.summed)^2,
            digits = 2)
)

p = ggplot(df,aes(x = density,log.as.summed)) + 
  geom_point(aes(color=species),size=3) + 
  geom_smooth(method = "lm",se = T) 
p = p + annotate("text", x = 8, y = 8, label = paste("R2:",r2,sep = " "))
p

#############################################
# Spearman correlation (only non-zero values)
#############################################
df.nonzeros = filter(df,density > 1) # filter zero values (that were set to 1 before)

# calculate R squared
r2.nonzeros = with(data = df,
          round(
            cor(density,log.as.summed)^2,
            digits = 2)
)

p.nonzeros = ggplot(df.nonzeros,aes(x = density,log.as.summed)) + 
  geom_point(aes(color=species),size=3) + 
  geom_smooth(method = "lm",se = T) 
p.nonzeros = p.nonzeros + annotate("text", x = 8, y = 8, label = paste("R2:",r2.nonzeros,sep = " "))
p.nonzeros


############
# Save plots
############
ggsave(filename = "FigureS4_spearman_trichome_vs_acylsugars/FigureS4A.all.pdf",
       plot = p,
       width = 7,
       height = 5)


ggsave(filename = "FigureS4_spearman_trichome_vs_acylsugars/FigureS4B.nonzeros.pdf",
       plot = p.nonzeros,
       width = 7,
       height = 5)

#################
# session info
##################
writeLines(capture.output(sessionInfo()), "FigureS4_spearman_trichome_vs_acylsugars/sessionInfo.txt")




