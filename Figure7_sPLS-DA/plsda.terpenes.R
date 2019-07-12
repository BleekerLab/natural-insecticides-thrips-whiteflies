library(mixOmics)
library(dplyr)



#############
# Load data
#############

# load data
# remove unncessary columns
df = read.delim("Figure7_sPLS-DA/pheno_terpenoids.tsv",header = T,stringsAsFactors = F,check.names = F)
row.names(df)=df$sample

##############
# PCA analysis
##############
X = df %>% select(-c(sample,wf,thrips))
pca.analysis = pca(X,ncomp = 10,center = T,scale = T)
plot(pca.analysis) 

# extracts the number of components necessary to explain more than 50% of the variance 
nComp = as.integer(which(pca.analysis$cum.var > 0.5)[1]) 

plotIndiv(object = pca.analysis,comp = c(1,2),group = df$wf,ind.names = T,legend = T,title = "PCA analysis (terpenes) - PC1 & PC2")
plotIndiv(object = pca.analysis,comp = c(1,3),group = df$wf,ind.names = T,legend = T,title = "PCA analysis (terpenes) - PC1 & PC3")
plotIndiv(object = pca.analysis,comp = c(2,3),group = df$wf,ind.names = T,legend = T,title = "PCA analysis (terpenes) - PC2 & PC3")



###################################
# PLS-DA on whitefly toxicity class
###################################
Y = df$wf

plsda.analysis <- plsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later
plotIndiv(plsda.analysis , comp = c(1,2),
          group = Y, 
          ind.names = FALSE, 
          ellipse = TRUE,
          legend = TRUE,
          title = 'PLSDA on terpenes & whitefly categories')



# performance evaluation
set.seed(2543) # for reproducibility
perf.plsda <- perf(plsda.analysis, dist = "max.dist",validation = "Mfold", folds = 5,progressBar = TRUE, auc = TRUE, nrepeat = 10) 

# error rates
plot(perf.plsda,measure = "BER", col = color.mixo(5), sd = TRUE, legend.position = "horizontal")


#########
# sPLS-DA
#########
# grid of possible keepX values that will be tested for each component
list.keepX <- seq(from = 10, to = ncol(X), by = 5)

tune.results <- tune.splsda(X,
                            Y, 
                            ncomp = 3,
                            validation = 'loo',
                            folds = 5, 
                            progressBar = F,
                            dist = 'max.dist', 
                            measure = "BER",
                            test.keepX = list.keepX,
                            nrepeat = 1)


#################
# session info
##################
writeLines(capture.output(sessionInfo()), "Figure7_sPLS-DA/sessionInfo.txt")