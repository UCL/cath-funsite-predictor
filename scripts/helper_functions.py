
import shap
import pandas as pd
import numpy as np
from collections import OrderedDict

import matplotlib.pyplot as plt
import warnings
import seaborn as sns
sns.despine()

from sklearn import preprocessing
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import binarize

from sklearn.pipeline import Pipeline

from sklearn import metrics
from sklearn.metrics import log_loss
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import matthews_corrcoef

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import SelectFromModel
from sklearn.externals import joblib

from sklearn import svm
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier

import xgboost as xgb
from xgboost import XGBClassifier
from xgboost import plot_importance
from xgboost import plot_tree

def groupshufflesplit_trainingdata(X, y, groups, testsize, randomstate):
    
    train_inds, test_inds = next(GroupShuffleSplit(test_size=testsize, random_state=randomstate).split(X, groups=groups))
    X_train, X_test, y_train, y_test = X.iloc[train_inds], X.iloc[test_inds], y.iloc[train_inds], y.iloc[test_inds]
    
    return(X_train, X_test, y_train, y_test, train_inds)

def preprocess_dataframes(df, annotation_col, ratio):
    
    '''If ratio is 100, use all negative samples'''
    df[annotation_col] = df.apply(lambda x: x[annotation_col], axis=1)
    print('#Prediction annotation:', annotation_col)

    # Get the list of domains that have site annotations
    site_data = df[(df[annotation_col] == 1)]
    positive_sample_num = site_data.shape[0]
    print('#Postive samples:',positive_sample_num)
    
    # make a subset dataframe containing only NOSITE (negative) residues
    nonsite_data = df[(df[annotation_col] == 0)]
    negative_sample_num = nonsite_data.shape[0]
    print('#Negative samples:',negative_sample_num)
    
    print ('Use these samples for training the model:')
    use_sample_num = positive_sample_num * ratio
    print ('- Used negative samples:',use_sample_num)
    if (ratio == 100):
        use_sample_num=negative_sample_num
    print ('- Used negative samples:',use_sample_num)
    if (use_sample_num > negative_sample_num):
        use_sample_num = negative_sample_num
        
    total_samples = use_sample_num + positive_sample_num
    print ('- Used negative samples:',use_sample_num)
    print ('- Total samples:',total_samples)
    
    dom_groups_df = nonsite_data.groupby(['dom_group']).size()
    dom_group_num = dom_groups_df.shape[0]
    print ('- No. of groups of samples:',dom_group_num)
    sample_size = round(use_sample_num/dom_group_num)
    print ('- Min. sample size in no_sites:', sample_size)
    if (sample_size > 40):
        sample_size=40
    print ('- Min. sample size in no_sites:', sample_size)
    
    nonsite_data = nonsite_data.groupby(['dom_group']).filter(lambda x: len(x) > sample_size)
    nonsite_data_randomsubset=nonsite_data.groupby('dom_group').apply(lambda x: x.sample(n=sample_size, random_state=10)).reset_index(drop=True)
    
    # COMBINE selected csa and non-csa data for the desired dataset ratio
    frames = [site_data, nonsite_data_randomsubset]
    concatenated_feature_data = pd.concat(frames)
    dataset_sample_num = concatenated_feature_data.shape[0]
    feature_data_ML = concatenated_feature_data.set_index('index').sample(n=dataset_sample_num, random_state=10)
    feature_data_ML.index.name = None
    
    return(feature_data_ML)
  
#get validation dataset in list

def plot_pr_curve(X_test, y_test, model, labelname, plotfilename):
    
    plt.rcParams['figure.figsize'] = (6.0, 6.0)
    plt.rcParams['axes.facecolor']='white'
    precision, recall, thresholds = precision_recall_curve(y_test, model.predict_proba(X_test)[:, 1])
    plt.plot(recall, precision, label=labelname)
    
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.legend(loc="best")
    #plt.title("Catalytic Site Prediction")
    plt.savefig(plotfilename, dpi=300, bbox_inches='tight')
    plt.show()
    
def get_validation_list(file):
    validation_domlist=[]
    infile = open(file,'r')
    
    for line in infile:
        list=line.split('\t')
        domain=list[0]
        #print(domain)
        validation_domlist.append(domain)

    infile.close()
    
    #print(validation_domlist)
    #len(validation_domlist)
    
    return(validation_domlist)

# plot feature importance

def plot_xgboost_inbuilt_feature_importance(model, plot_title, image_name):
    plt.rcParams['figure.figsize'] = (15.0, 15.0)
    plt.rcParams.update({'font.size': 20})
    
    plot_importance(model, max_num_features=50, importance_type='weight', show_values=False, height = 0.4)
    
    plt.title(plot_title)
    plt.savefig(image_name, dpi=300, bbox_inches='tight')
    plt.show()
    
def plot_xgboost_tree(model, plot_filename):
    plot_tree(model)
    fig = plt.gcf()
    plt.rcParams['figure.figsize'] = (5.0, 5.0)
    fig.set_size_inches(50, 50)
    fig.savefig(plot_filename)
    
def reduce_features_with_no_effect(modelname, model_features):
    
    features=model_features
    model_features=[]
    
    DD = dict(zip(features, modelname.feature_importances_))
    for item, score in DD.items():
        #print(item, score)
        if(score > 0):
            model_features.append(item)

    return model_features

def plot_shap_importance(model, X, plot_filename):
    
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)
    shap.summary_plot(shap_values, X, plot_type="bar", show=False)
    plt.rcParams['figure.figsize'] = (5.0, 5.0)
    plt.tight_layout()
    plt.savefig(plot_filename)

    
# plot feature characteristics

def plot_feature_characteristics_violin(feature_name, annotation, df, feature_label, annotation_label):
    plt.rcParams['figure.figsize'] = (5.0, 5.0)
    
    fig=sns.violinplot(x=annotation, y=feature_name, data=df)
    
    fig.set(xlabel=annotation_label, ylabel=feature_label)

    plt.show()
    
    
def perform_grid_search(classifier, params_test, score_type, X_train, y_train, features, cv_num):
    #score_type = f1_score, roc_auc
    
    gsearch = GridSearchCV(estimator = classifier, param_grid = params_test, scoring=score_type,iid=False, cv=cv_num)
    gsearch.fit(X_train[features], y_train)
    
    return(gsearch.grid_scores_, gsearch.best_params_, gsearch.best_score_)

def preprocess_dataframes_sites(site_data, nonsite_data, ratio):
    '''If ratio is 100, use all negative samples'''
    
    # Get the list of domains that have site annotations
    positive_sample_num = site_data.shape[0]
    print('#Postive samples:',positive_sample_num)
    
    # make a subset dataframe containing only NOSITE (negative) residues
    negative_sample_num = nonsite_data.shape[0]
    print('#Negative samples:',negative_sample_num)
    
    print ('Use these samples for training the model:')
    use_sample_num = positive_sample_num * ratio
    if (ratio == 100):
        use_sample_num=negative_sample_num
        
    if (use_sample_num > negative_sample_num):
        use_sample_num = negative_sample_num
    total_samples = use_sample_num + positive_sample_num
    print ('- Used negative samples:',use_sample_num)
    print ('- Total samples:',total_samples)
    
    dom_groups_df = nonsite_data.groupby(['dom_group']).size()
    dom_group_num = dom_groups_df.shape[0]
    print ('- No. of groups of samples:',dom_group_num)
    sample_size = round(use_sample_num/dom_group_num)
    #if (sample_size > 40):
        #sample_size=40
    print ('- Min. sample size in no_sites:', sample_size)
    
    nonsite_data = nonsite_data.groupby(['dom_group']).filter(lambda x: len(x) > sample_size)
    nonsite_data_randomsubset=nonsite_data.groupby('dom_group').apply(lambda x: x.sample(n=sample_size, random_state=10)).reset_index(drop=True)
    
    # COMBINE selected csa and non-csa data for the desired dataset ratio
    frames = [site_data, nonsite_data_randomsubset]
    concatenated_feature_data = pd.concat(frames)
    dataset_sample_num = concatenated_feature_data.shape[0]
    feature_data_ML = concatenated_feature_data.set_index('index').sample(n=dataset_sample_num, random_state=10)
    feature_data_ML.index.name = None
    
    return(feature_data_ML)

def preprocess_dataframes_sites_all(site_data, nonsite_data):  
    # COMBINE selected csa and non-csa data
    frames = [site_data, nonsite_data]
    concatenated_feature_data = pd.concat(frames)
    dataset_sample_num = concatenated_feature_data.shape[0]
    feature_data_ML = concatenated_feature_data.set_index('index')
    #feature_data_ML.index.name = None

def get_pr_curve_pred_anno_dataset(X, y, groups, pred_anno_column_name):
    
    y_real=[]
    y_proba=[]
    
    i = 0
    
    for train, test in group_kfold.split(X, y, groups):
        
        print("FOLD:", i+1)
        print("TRAIN:", train, "TEST:", test)

        X_train = X.iloc[train]
        y_train = y.iloc[train]
        X_test = X.iloc[test]
        y_test = y.iloc[test]
        
        pred_proba = X_test[pred_anno_column_name]

        # Compute PR curve and area the curve
        precision, recall, thresholds = precision_recall_curve(y_test, pred_proba)

        lab = 'Fold %d AUC=%.4f' % (i+1, auc(recall, precision))
        print(lab)
        
        #plt.plot(recall, precision, lw=1, alpha=0.3, label=lab)
        
        y_real.append(y_test)
        y_proba.append(pred_proba)

        i += 1
        
    y_real = np.concatenate(y_real)
    y_proba = np.concatenate(y_proba)

    precision_overall, recall_overall, _ = precision_recall_curve(y_real, y_proba)

    pr_auc_overall = auc(recall, precision)
    
    return (precision_overall, recall_overall, pr_auc_overall)

def get_pr_curve_model_pred_dataset(X, y, groups, classifier):
    
    y_real=[]
    y_proba=[]
    
    i = 0
    
    for train, test in group_kfold.split(X, y, groups):
        
        print("FOLD:", i+1)
        print("TRAIN:", train, "TEST:", test)

        X_train = X.iloc[train]
        y_train = y.iloc[train]
        X_test = X.iloc[test]
        y_test = y.iloc[test]

        classifier.fit(X_train, y_train)
        
        pred_proba = classifier.predict_proba(X_test)

        # Compute PR curve and area the curve
        precision, recall, thresholds = precision_recall_curve(y_test, pred_proba[:, 1])

        lab = 'Fold %d AUC=%.4f' % (i+1, auc(recall, precision))
        print(lab)
        
        #plt.plot(recall, precision, lw=1, alpha=0.3, label=lab)
        
        y_real.append(y_test)
        y_proba.append(pred_proba[:,1])

        i += 1
        
    y_real = np.concatenate(y_real)
    y_proba = np.concatenate(y_proba)

    precision_overall, recall_overall, _ = precision_recall_curve(y_real, y_proba)

    pr_auc_overall = auc(recall, precision)
    
    return (precision_overall, recall_overall, pr_auc_overall)