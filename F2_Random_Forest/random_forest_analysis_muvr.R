##############
# Libraries
##############
library("MUVR")
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("doParallel"))

#########
# params
#########
# fixed params
n_cores = 3
n_outer = 7
n_inner = 4
n_permutations = 10
model_choosen = "min"             # Choose between "min", "mid" or "max" 

# grid search params
start_ratio <- 0.7
end_ratio <- 0.7
step_ratio <- 0.1 

start_nreps <- 5
end_nreps <- 10
step_nreps <- 5

# cluster generation
cl = makeCluster(n_cores)   
registerDoParallel(cl)


##########
# Input data
##########
metabolites = read.delim("F2_Random_Forest/volatiles/20180625_F2metabolite_table.tsv",
                         header = T,
                         stringsAsFactors = F,
                         check.names = F)


whitefly <- read.delim("F2_Random_Forest/whitefly/whitefly_table_for_random_forest.tsv",
                            header = T,
                            stringsAsFactors = F) %>% 
  group_by(line) %>% 
  dplyr::summarise(percentage = median(percentage))

df <- inner_join(metabolites, whitefly, by = "line")

###############################################
# MUVR Random Forest search for best parameters
###############################################
IDs <- df$line
X = df %>% dplyr::select(- line, - percentage)  
Y = df$percentage

# hyperparameter grid search --> same as above but with increased mtry values
hyper_grid <- expand.grid(
  ratio     = seq(from = start_ratio, 
                  to   = end_ratio, 
                  by   = step_ratio),
  rep       = seq(from = start_nreps, 
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
optimization_plot <- hyper_grid %>% 
  mutate(ratio = as.factor(ratio), 
         rep = as.factor(rep)) %>% 
  rownames_to_column("combination") %>% 
  mutate(combination = as.numeric(combination)) %>% 
  ggplot(., aes(combination, y = q2, color = rep)) + 
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

par(mar=c(4,4,1,1))
plotVAL(rf_model)
plotMV(rf_model, model = "min")

par(mar=c(3,10,1,1))
plotVIP(rf_model)

#######################
# Permutation analysis
#######################

### Part 1: permutations ###
# compute permuted models 
# collect Q2 from permutations
perm_fit = numeric(n_permutations)

# computed permuted models
# collect permuted feature importances
features_permuted_pvalues_matrix <- matrix(0, 
                                           nrow = ncol(X), # One row = one feature
                                           ncol = n_permutations)

rownames(features_permuted_pvalues_matrix) = colnames(X)
colnames(features_permuted_pvalues_matrix) = paste0("permutation",seq(1:n_permutations))

# permutations
cat("\nStarting permutations")
for (p in 1:n_permutations) {
  cat('\nPermutation',p,'of', n_permutations)
  YPerm = sample(Y)
  
  perm = MUVR(X         = X, 
              Y         = YPerm,
              ID =        IDs,
              nRep      = best_params$rep,
              nOuter    = n_outer,
              varRatio  = best_params$ratio,
              scale     = FALSE, 
              DA        = FALSE, 
              fitness   = "RMSEP", 
              method    = "RF", 
              parallel  = TRUE)
  # for model
  perm_fit[p] = perm$fitMetric$Q2
  
  # for each variable
  features_permuted_pvalues_matrix[,p] = as.vector(perm$VIP[,"min"])
}

### Part 2: plot of RF model (actual fit vs permuted fits) ###
cat("\nCreating permutation plot for Random Forest model")
# actual (original RF mode)
actual_fit <- rf_model$fitMetric$Q2[1] #(1 = "min model")

# Parametric (Student’s) permutation test significance
pvalue <- pPerm(actual = actual_fit, 
                h0 = perm_fit,
                side = "greater",
                type = "t")

plotPerm(actual = actual_fit, 
         xlab = "Q2 metric",
         h0 = perm_fit) 


perm_fit_df = data.frame(permutation = seq(1:length(perm_fit)), 
                         q2 = perm_fit) 

model_permutation_plot <- ggplot(perm_fit_df, aes(x = q2)) + 
  geom_histogram(bins = 10) + 
  geom_vline(xintercept = actual_fit, col = "blue") + 
  labs(x = "Q2 metric", y = "Frequency") +
  ggtitle(paste("Distribution of Q2 p-values based on \n",
                n_permutations,
                "permutations of the Y variable \n p-value = ",
                format(pvalue,digits = 3, decimal.mark = "."),sep = " "
                ))


### Part 3: extract significant p-values for each variable ###
features_permuted_pvalues_matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("feature") %>% 
  pivot_longer(cols = - feature) %>% 
  split(.$feature) %>% 
  

mutate(original_vip = rf_model$VIP[,model_choosen]) %>% 


# plot
features_permuted_pvalues_df %>% 
  pivot_longer(cols = - feature, names_to = "permutation", values_to = "var_imp") %>% 
  ggplot(., aes(x = var_imp)) +
  geom_histogram() +
  facet_wrap(~ feature)
features_permuted_pvalues_df



#########
# Save
########
save(df,
     X,
     Y,
     rf_model,
     best_params,
     optimization_plot,
     model_permutation_plot,
     file = "rf_analysis.RData",
     compress = "gzip",
     compression_level = 6)

stopCluster(cl)