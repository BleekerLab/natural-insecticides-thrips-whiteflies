"""
Module containing helper functions for the Random Forest analysis
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
import random
import seaborn as sns


###################################
# define main random forest routine 
###################################

def single_random_forest_run(X,y,rs,disp=False,nb_of_splits = 6,nb_of_trees=1000,njobs=2):
    """
    This function takes a feature matrix (X) and a label array (y) and fits 
    a number of random forest classifiers based on the specified number of splits and trees.
    One can fix the random state /rs/ to consistently generate the same models between runs. 

    Arguments
    y: should be a list.
    X: should be a dataframe (index contains rows id, columns contains variable values)
    """
  
    # create dataframes to accomodate:
      # computed values for each variable importance (variableImportance)
      # the actual and predicted y values
      # the normalised accuracy score
    variableImportance = pd.DataFrame(np.zeros([X.shape[1],nb_of_splits]))    
    dfy = pd.DataFrame(y,columns=['tox'])
    yhat = pd.DataFrame(['']*len(y),columns=['predicted'],dtype=np.str)

    # Stratified K-fold cross validator: provides train/test indices to split data in train/test sets
    kfold = StratifiedKFold(n_splits=nb_of_splits,random_state=rs)
    
    splitnr = 0
    
	# stratify info
    stratify_info = pd.DataFrame(data=y,columns=["y"])
    stratify_info = stratify_info["y"].map({'non-toxic': 0, 'toxic': 1})

    # Splits the X and y arrays n times (as specified by number of splits)
    # Initialise a Random Forest classifier with the specified number of trees and using class weights. 
    # Random state is fixed so that results are reproducible from run to run.
    for train_index, test_index in kfold.split(X,y):
        nt = 1/sum(stratify_info[train_index]==0) # calculates a weight for the non-toxic class (corrects for class imbalance)
        tx = 1/sum(stratify_info[train_index]==1) # calculates a weight for the toxic class (corrects for class imbalance)   
        rf = RandomForestClassifier(n_estimators=nb_of_trees,
                                    class_weight={"toxic":tx,"non-toxic":nt},
                                    random_state=rs,
                                    n_jobs=njobs
                                    )

        # train the Random Forest model 
        rf = rf.fit(X.iloc[train_index,:],dfy.iloc[train_index].values.ravel())

        # predict class values and adds it to the yhat dataframe
        y_pred = rf.predict(X.iloc[test_index,:]) # this yields ["toxic","non-toxic","non-toxic","toxic"] (same length as test_index)
        for j in range(len(y_pred)):    
            yhat.iloc[test_index[j]] = y_pred[j]  

        # extracts feature importance for each run
        variableImportance.iloc[:,splitnr] = rf.feature_importances_
        
        splitnr+=1

    # For each dataset split, calculate the accuracy score 
    # predict the accuracy score: takes a vector of y_true and compares to y_pred
    accuracyScore = round(accuracy_score(yhat["predicted"],y),ndigits=2)

    if disp:
        print(classification_report(y, yhat))
        pdata = variableImportance.copy()
        pdata.rename(columns={i:'split{0}'.format(i) for i in range(6)},inplace=True)
        pdata['x'] = range(variableImportance.shape[0])
        pdata.head()
        sns.lineplot(data=pdata)
    
    return [variableImportance,yhat,accuracyScore] # two dataframes and one float number


############################################
# Extracts variable importance and predicted 
############################################

def extract_feature_importance_avg_and_sd_from_multiple_random_forest_runs(
        X,
        y,
        nb_of_splits=6,
        nb_of_trees=1000,
        nb_of_runs=5,
        njobs=2
        ):
    """
    This function runs the single_random_forest_run multiple times to assess the effect 
    of the different splits on each feature importance.
    It returns a dataframe with the average and standard deviation for each feature importance.

    Arguments:
    X: A Pandas dataframe containing the feature values (index contains rows id, columns contains variable values). 
    y: a list containing the labels 
    nb_of_splits: the number of splits used for stratified K-folds cross-validator.
    nb_of_trees: number of trees in each single RandomForest model.
    nb_of_runs: the number of times the single random forest will be run (with the specified number of splits)

    Returned value: a dataframe containing each feature importance for each single run 
    and the average and standard deviation.
    The dataframe has nb_of_runs columns (+ 2 columns for average and std)
    The dataframe has X.shape[1] rows (one row per feature).


    """
    # creating column names for the result dataframe
    colNames = ["run" + str(i) for i in range(nb_of_runs)] 

    # initialize two dataframe that will accomodate each run feature importance averages and standard deviations.
    feature_importance_averages = pd.DataFrame(
        index=X.columns.tolist(),
        columns=colNames)
    feature_importance_sd = pd.DataFrame(
        index=X.columns.tolist(),
        columns=colNames)

    # For each run:
    #  1) A single random forest is computed. 
    #  2) Each feature has nb_of_splits importances. 
    #  3) Average and standard deviation are computed and added in their corresponding final dataframes. 
    for i in range(nb_of_runs):
        single_run = single_random_forest_run(
            X,y,rs=i,disp=False,nb_of_splits = nb_of_splits,nb_of_trees=nb_of_trees,njobs=njobs)[0]
        
        feature_importance_averages.loc[:,"run" + str(i)] = single_run.mean(axis=1).tolist()    
        feature_importance_sd.loc[:,"run" + str(i)] = single_run.std(axis=1).tolist()     

    return [feature_importance_averages,feature_importance_sd]


#######################
# Permutation functions
#######################

def extract_feature_importances_from_random_forests_on_permuted_y(X,y,nperm=100,randomNumber=random.SystemRandom(),nb_of_splits=6,nb_of_trees=1000,njobs=2):
    """
    This function runs a single random forest analysis for the specified number of permutations. 
    It does the shuffling of the y label for the specified number of 
    It returns

    Arguments:
    X: A Pandas dataframe containing the feature values (index contains rows id, columns contains variable values). 
    y: a list containing the labels.
    n_perm: number of permutations.
    randomNumber: this has to be a number used to suffle  
    nb_of_splits: the number of splits used for stratified K-folds cross-validator.
    nb_of_trees:  number of trees in each single RandomForest model.


    """

    # functions used to create a generator of permutations
    def perm_generator(seq):
        """
        Returns a generator that can be used to iterate over possible permutations.
        It takes a Python list as an input
        """
        seen = set()
        length = len(seq)
        while True:
            perm = tuple(random.sample(seq, length)) # produces a random draw of same length as initial list
            if perm not in seen:
                seen.add(perm)
                yield perm

    def get_perms(seq, N=nperm):
        """
        This function encapsulate the previous one (perm_generator)
        It returns a list with all possible permutations given a list + a number of permutations
        """
        rand_perms = perm_generator(seq)
        return [next(rand_perms) for _ in range(N)]

    ##################################
    # Random Forest runs on permuted y
    ##################################

    # Make a dataframe to accomodate the permuted feature importance
    colNames = ["perm" + str(perm) for perm in range(nperm)]
    permuted_feature_importance_df = pd.DataFrame(columns = colNames, index=X.columns.tolist())

    # generate nperm y permutations and run a single Random Forest on it
    y_permutations = list(get_perms(y,N=nperm))

    for i in range(nperm):
        # get one of the possible permutations
        y_permuted = list(y_permutations[i])

        # run a single Random Forest analysis (arbitrarily fixed the random state)
        feature_importance_from_permuted_y = single_random_forest_run(X,y_permuted,rs=1234,nb_of_splits=nb_of_splits,nb_of_trees=nb_of_trees,njobs=njobs)[0]
        # The feature importances from the different split are averaged. 
        # This average feature importance is added to the dataframe as column nperm

        permuted_feature_importance_df["perm" + str(i)] = feature_importance_from_permuted_y.mean(axis=1).tolist()

    return permuted_feature_importance_df


###############################################################################
# Get a dataframe of features with a significant p-value (based on permutations)
###############################################################################


def iperc(feature_importance_original,feature_importance_permuted):
    """
    Function to determine pvalues of average feature importance values based on results from permuted feature importance values.

    Arguments:
    feature_importance_original: a 1D numpy array (Pandas series) containing the averaged feature importance values (original dataset)
    feature_importance_permuted: a Pandas dataframe containing each feature importances from each permutations.

    Returns a Pandas dataframe with variables as index and one column containing the p-values.
    """
    


    df = pd.DataFrame(index=feature_importance_permuted.index,columns=['p-value'])
    
    for i in range(len(feature_importance_original)):            
        pn = sum(feature_importance_permuted.iloc[i,:] < feature_importance_original.iloc[i]) / feature_importance_permuted.shape[0]
        pp = sum(feature_importance_permuted.iloc[i,:] >= feature_importance_original.iloc[i]) / feature_importance_permuted.shape[0]
    
        df.iloc[i] = min(pn,pp)
            
    return df


 
def get_significant_features(X,original_feature_importances,permuted_feature_importances,pval=0.05):
    """

    Arguments:
    X: A Pandas dataframe containing the feature values (index contains rows id, columns contains variable values).
    original_feature_importances: the output of the extract_feature_importance_avg_and_sd_from_multiple_random_forest_runs function.
    permuted_feature_importances: the output of the extract_feature_importances_from_random_forests_on_permuted_y function.
    pval: a p-value threshold to filter the features based on the computed p-values.

    Returns a pandas dataframe with only the significant features (pvalue < pval) and 
    with additional metrics useful for filtering (average, standard deviation and pooled standard deviation)
    """
    
    # average feature importance from original dataset
    mean_varimportance = original_feature_importances[0].mean(axis=1)

    # how many times this original feature importance was lower or higher than a set of random feature importances?
    df = iperc(
    	feature_importance_original = original_feature_importances[0].mean(axis=1),
    	feature_importance_permuted = permuted_feature_importances
    	)

    # add the average en std deviations to the dataframe
    df["average"] = original_feature_importances[0].mean(axis=1).tolist()
    pooled_std = np.sqrt((original_feature_importances[1].mean(axis=1)**2))
    df["sd"] = pooled_std.tolist()
    df["rsd"] = pooled_std/mean_varimportance

    return df


####################
# Plotting functions
####################

def plot_candidate_heatmap(df,class_col = "wf", significant="signif.index.values",log=True):
    """
    Plots a heatmap from a dataframe containing a column
    with the phenotypic class and the other columns containing the measured variables.

    Arguments: 
    df: a Pandas dataframe with the first column containing the pheno class and the rest being measured variables.
    class_col: the name of the category column from which a new index will be made.
    significant: the index values from the dataframe returned by the get_significant_features function.
    log: True or False. If True, 

    """
    # sort dataframe by the column class to separate categories
    df_sorted = df.sort_values(by=[class_col])

    # make a new column to host the new index + reindex
    df_sorted["new_index"] = df_sorted.index.values + "_" + df_sorted[class_col]
    df_sorted = df_sorted.set_index(df_sorted["new_index"])
    
    # remove column with class info
    # remove new_index column
    df_sorted = df_sorted.drop(class_col,axis=1)
    df_sorted = df_sorted.drop("new_index",axis=1)

    # log correction
    if log ==  True:
        # replace 0 by 1
        df_sorted_replaced = df_sorted.replace(0,1) # log2(1) = 0
        # apply log2 function to every value
        log_function = lambda x: np.log2(x)
        df_sorted_replaced_logged = df_sorted_replaced.apply(log_function)
    else: 
        pass # do nothing

    # keep only significant variables
    df_final = df_sorted_replaced_logged.loc[:,significant]

    # final plot
    return sns.heatmap(df_final,cmap="Blues")

def plot_candidate_sum_by_class(df,class_col="wf",significant="signif.index.values",ax=None):
    """
    Plots the summed abundance per class of the selected candidates.

    df: df: a Pandas dataframe with the first column containing the class info and the rest being measured variables.
    class_col: the name of the column from which the sums will be computed (sum over class A // sum over class B)
    significant: the index values from the dataframe returned by the get_significant_features function. 
    """
    # how many unique class levels?
    classes = list(set(df.loc[:,class_col].tolist()))

    # filters out the non-significant variables
    df_filt = df.loc[:,significant]

    # add the class column as the last column
    df_filt[class_col] = df[class_col]

 	# sum over each category (ignoring the last column that contains the class info)
    if len(classes) == 2:
    	classA= classes[0]
    	classB = classes[1]
    	classA_sum = df_filt.iloc[:,:-1][df_filt[class_col] == classA].sum() 
    	classB_sum = df_filt.iloc[:,:-1][df_filt[class_col] == classB].sum()
    else:
    	print("no more than two levels for a category implemented at this time")

    # prepare a dataframe for plotting
    sumByClass = pd.DataFrame({classA:classA_sum, classB:classB_sum})
    sumByClass = sumByClass.reset_index()
    sumByClass.columns.values[0]="candidate"
    sumByClass = sumByClass.melt(id_vars="candidate",var_name="class",value_name="value")

    fig, ax = plt.subplots(figsize=(12,10))
    ax.set_xscale("log",basex=2)
    return sns.barplot(x="value",y="candidate",hue="class",data=sumByClass,ax=ax)

