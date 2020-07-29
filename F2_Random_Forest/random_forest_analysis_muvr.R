##############
# Libraries
##############
library("MUVR")
library("tidyverse")
library("doParallel")

#########
# params
#########
# fixed params
n_cores = 3
n_outer = 7
n_inner = 4

# grid search params
start_ratio <- 0.6
end_ratio <- 0.9
step_ratio <- 0.1 

start_nreps <- 5
end_nreps <- 50
step_nreps <- 5

# cluster generation
cl = makeCluster(n_cores)   
registerDoParallel(cl)


##########
# Input data
##########
metabolites = read.delim("F2_Random_Forest/Figure8_volatiles_F2/20180625_F2metabolite_table.tsv",
                         header = T,
                         stringsAsFactors = F,
                         check.names = F)


whitefly <- read.delim("F2_Random_Forest/Figure7_F2_whitefly_survival/whitefly_table_for_random_forest.tsv",
                            header = T,
                            stringsAsFactors = F)

df <- inner_join(metabolites, whitefly, by = "line")

###############################################
# MUVR Random Forest search for best parameters
###############################################
IDs <- df$line
X = df %>% dplyr::select(- line, - percentage)  
Y = df$percentage

# hyperparameter grid search --> same as above but with increased mtry values
hyper_grid <- expand.grid(
  ratio     = 
    seq(from = start_ratio, 
        to   = end_ratio, 
        by   = step_ratio),
  rep       = 
    seq(from = start_nreps, 
        to   = end_nreps, 
        by   = step_nreps),
  q2     = 0 # will be used to collect the Q2 metric for each RF model
)

for (i in 1:nrow(hyper_grid)){
  print(paste0("testing combination ", i, " out of ", nrow(hyper_grid)))
  
  rf_model <- MUVR(X = X, 
                   Y = Y,
                   ID = IDs,
                   nRep = hyper_grid$rep[i],
                   nOuter = n_outer,
                   nInner = n_inner,
                   varRatio = hyper_grid$ratio[i],
                   scale = FALSE, 
                   DA = FALSE, 
                   fitness = "RMSEP", 
                   method = "RF", 
                   parallel = TRUE)
  
  hyper_grid$q2[i] <- rf_model$fitMetric$Q2[1] # min model
}

# plot evolution of Q2 against parameters
hyper_grid %>% 
  mutate(ratio = as.factor(ratio), 
         rep = as.factor(rep)) %>% 
  rownames_to_column("combination") %>% 
  mutate(combination = as.numeric(combination)) %>% 
  ggplot(., aes(combination, y = q2, color = ratio)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0,1))

# extract the best parameters
best_params <- hyper_grid[which.max(hyper_grid$q2),]

###############################################
# MUVR Random Forest search *with* best parameters
###############################################

rf_model <- MUVR(X = X, 
                 Y = Y,
                 ID = IDs,
                 nRep = best_params$rep,
                 nOuter = n_outer,
                 nInner = n_inner,
                 varRatio = best_params$ratio,
                 scale = FALSE, 
                 DA = FALSE, 
                 fitness = "RMSEP", 
                 method = "RF", 
                 parallel = TRUE)

plotVAL(rf_model)
plotMV(rf_model, model = "min")

stopCluster(cl)


