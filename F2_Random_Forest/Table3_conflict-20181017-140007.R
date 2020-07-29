

############## read input data #########
# volatiles
metabolites = read.delim("Figure8_volatiles_F2/volatiles_sorted.tsv",header=T,stringsAsFactors = F,check.names = F)

# whitefly classification
toxic = read.delim("Figure7_F2_whitefly_survival/toxic_lines.tsv",header = T,stringsAsFactors = F)
nontoxic = read.delim("Figure7_F2_whitefly_survival/nontoxic_lines.tsv",header = T,stringsAsFactors = F)

lines_with_zero_surviving_wf = c(294,324,373,387,396) # these lines were not included in the WF phenotype analysis but are toxic (no surviving whiteflies!)

######## Data wrangling ######
# Add a class name to the F2 whitefly survival analyses
toxic$class = "R"
nontoxic$class = "S"

# combine these two tables
lines_with_a_pheno_class = rbind.data.frame(toxic,nontoxic)

# remove unnessary tables

########## intersection between volatiles and phenotype tables 
df = inner_join(metabolites,lines_with_a_pheno_class,by="line")


######## Extract input table for Random Forest analysis




lines_with_zero_surviving_wf = c(294,324,373,387,396)
