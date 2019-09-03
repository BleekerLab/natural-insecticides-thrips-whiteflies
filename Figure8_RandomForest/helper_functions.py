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


###################################
# define main random forest routine 
###################################

def single_random_forest_run(X,y,rs,disp=False,nb_of_splits = 6,nb_of_trees=1000):
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
    yhat = pd.DataFrame(['']*len(y),columns=['predict_tox'],dtype=np.str)
    accuracyScores = pd.DataFrame(['']*nb_of_splits,columns=["normalised_accuracy"],dtype=np.str)

    
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
                                    random_state=rs)

        rf = rf.fit(X.iloc[train_index,:],dfy.iloc[train_index].values.ravel())

        # predict class values and adds it to the yhat dataframe
        y_pred = rf.predict(X.iloc[test_index,:])
        for j in range(len(y_pred)):    
            yhat.iloc[test_index[j]] = y_pred[j]    
        
        # extracts feature importance for each run
        variableImportance.iloc[:,splitnr] = rf.feature_importances_

        # extracts accuracy score for each run
        accuracyScores.iloc[splitnr] = accuracy_score(y_pred, dfy.iloc[test_index,:])

        splitnr+=1

        
    if disp:
        print(classification_report(y, yhat))
        pdata = variableImportance.copy()
        pdata.rename(columns={i:'split{0}'.format(i) for i in range(6)},inplace=True)
        pdata['x'] = range(variableImportance.shape[0])
        pdata.head()
        sns.lineplot(data=pdata)
    
    return [variableImportance,yhat,accuracyScores] # three dataframes


 ############################
 # Extracts variable importance and predicted 
 ############################