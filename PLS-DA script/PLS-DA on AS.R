library(mixOmics)
library(mdatools)
library(tidyverse)

df = read.delim(file = "Random_Forest/phenotypes_vs_acylsugars.tsv")
wf.pheno <- df$wf
matrix <- df %>% select(-thrips, -wf ) %>% column_to_rownames(var = "sample")

# Check amount of components  
plsda.results <- splsda(matrix, wf.pheno, ncomp = 10)
set.seed(30)
plsda.performance <- perf(plsda.results, validation = "Mfold", folds = 3, 
     progressBar = FALSE, nrepeat = 50)
plot(plsda.performance, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")


list.keepX <- c(5:10,  seq(20, 100, 10))
tune.plsda<- tune.splsda(matrix, wf.pheno, ncomp = 10, 
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)

error <- tune.plsda$error.rate
ncomp <- tune.plsda$choice.ncomp$ncomp

select.keepX <- tune.plsda$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.plsda, col = color.jet(ncomp))

MyResult.splsda.final <- splsda(matrix, wf.pheno, ncomp = 2)

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")
selectVar(MyResult.splsda.final, comp=1)$value

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')

plotVar(plsda.results, cutoff = 0.65)
