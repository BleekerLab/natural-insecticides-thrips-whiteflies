# library
library(dplyr)

############## read input data #########
# volatiles
metabolites = read.delim("Figure8_volatiles_F2/volatiles_sorted.tsv",header=T,stringsAsFactors = F,check.names = F)

# whitefly classification
toxic = read.delim("Figure7_F2_whitefly_survival/toxic_lines.tsv",header = T,stringsAsFactors = F)
nontoxic = read.delim("Figure7_F2_whitefly_survival/nontoxic_lines.tsv",header = T,stringsAsFactors = F)

######## Data wrangling ######
# Add a class name to the F2 whitefly survival analyses
toxic$class = "R"
nontoxic$class = "S"

# combine these two tables
toxic= toxic %>% select(line,class) # remove coef, pvalue etc.
nontoxic= nontoxic %>% select(line,class) # remove coef, pvalue etc.
all_lines_with_pheno_class = rbind.data.frame(toxic,nontoxic)

########## intersection between volatiles and phenotype tables 
df = inner_join(metabolites,all_lines_with_pheno_class,by="line")

######## Extract input table for Random Forest analysis

# feature_vector
feature_vector = as.data.frame(df$class)
write.table(feature_vector,file = "Table3_RF_F2_wf_volatiles_candidates/feature_vector.txt",sep = "\t",quote = F,row.names = F,col.names = F)

# feature matrix 
feature_matrix = select(df,-c(class))
write.table(feature_matrix,file = "Table3_RF_F2_wf_volatiles_candidates/feature_matrix.txt",sep = "\t",quote = F,row.names = F)

# terpenoids
write.table(select(df,-c(class)),file = "Table3_RF_F2_wf_volatiles_candidates/volatiles.txt",sep = "\t",quote = F,row.names = F)

# whiteflies
write.table(select(df,line,class),file = "Table3_RF_F2_wf_volatiles_candidates/whiteflies.txt",sep = "\t",quote = F,row.names = F)

### Session info
writeLines(capture.output(sessionInfo()), "Table3_RF_F2_wf_volatiles_candidates/sessionInfo.txt")