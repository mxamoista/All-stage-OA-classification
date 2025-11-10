import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_curve, roc_auc_score,auc
from sklearn.model_selection import train_test_split,KFold,StratifiedKFold
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier   
from sklearn.metrics.pairwise import chi2_kernel
from sklearn.model_selection import train_test_split
import warnings

def feature_selection_rf(data, label, feature_num):
    """
    use random forest to select feature
    :param data:  data [dataframe]
    :param label: label column
    :param feature_num: the number of selected feature
    :return: selected data [dataframe]
    """
    rfc = RandomForestClassifier(n_estimators=1500)
    feature_importance = rfc.fit(data, label).feature_importances_
    sort_ind = np.argsort(feature_importance, )[::-1]
    reduce_col = [data.columns[sort_ind[i]] for i in range(feature_num)]
    sel_data = data.loc[:, reduce_col]
    return sel_data

def get_cv_score_svm(data,label,c):
    n_splits = 5
    folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=2022)
    f1_test, f1_train = [],[]
    pre_test, pre_train = [],[]
    rec_test, rec_train = [],[]
    auc_test = []

    # data,label = data.iloc[:,:-1], data.iloc[:,-1]
    for train_ind, test_ind in folds.split(data, label):
        x_train, y_train = data.iloc[train_ind,:], label.iloc[train_ind]
        x_test, y_test = data.iloc[test_ind,:], label.iloc[test_ind]
        if gamma != None:
            k_train = chi2_kernel(x_train,gamma=gamma)
            k_test = chi2_kernel(x_test,x_train,gamma=gamma)
        else:
            k_train = chi2_kernel(x_train)
            k_test = chi2_kernel(x_test,x_train)
        
        svm = SVC(C=c,kernel='precomputed',cache_size=200,probability=True,class_weight='balanced').fit(k_train.T,y_train)
        y_pred = svm.predict(k_test)
        y_pred_prob = svm.predict_proba(k_test)
        y_train_pred = svm.predict(k_train)
        y_train_pred_prob = svm.predict_proba(k_train)[:,1]


        f1_test.append(f1_score(y_test, y_pred, average='micro'))
        f1_train.append(f1_score(y_train, y_train_pred, average='micro'))
        pre_test.append(precision_score(y_test, y_pred, average='micro'))
        pre_train.append(precision_score(y_train, y_train_pred, average='micro'))
        rec_test.append(recall_score(y_test, y_pred, average='micro'))
        rec_train.append(recall_score(y_train, y_train_pred, average='micro'))     
        auc_test.append(roc_auc_score(y_test, y_pred_prob,multi_class='ovr'))

    f = sum(f1_test) / n_splits
    p = sum(pre_test) / n_splits
    r = sum(rec_test) / n_splits
    auc = sum(auc_test) / n_splits
    return round(p,3),round(r,3),round(f,3),round(auc,3),c

def data_normalization(data):
    scaler = preprocessing.MinMaxScaler()
    norm_data = scaler.fit_transform(data)
    #norm_data = preprocessing.normalize(norm_data, norm='l2')
    norm_data = pd.DataFrame(norm_data,columns = data.columns)
    return norm_data

# main
gamma = None #if set gamma, use kernel function
c = 3 #SVM's hyperparameter
feature_num = 300
iter_num = 50
best_col = {} #chosen features and their chosen times
aucs = []

# Load data
matrix = pd.read_csv('./data/combat_edata.txt', sep='\t', index_col=0)
annotation_col = pd.read_csv('./data/annotation_col.txt', sep='\t', index_col=0)
# Combine the two dataframes
combi_data = pd.concat([matrix.T, annotation_col], axis=1)
# Map subtype to numerical values for clustering
subtype_mapping = {'TC1': 0, 'TC2': 1, 'TC3': 2, 'TC4': 3, 'TC5': 4, 'TC6': 5}
combi_data['subtype'] = combi_data['Subtype'].map(subtype_mapping)
# Separate features and target
X = combi_data.drop(['Subtype', 'subtype'], axis=1)
y = combi_data['subtype']
# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
# Save the train/test sets
X_train.to_csv('X_train.txt', sep='\t')
X_test.to_csv('X_test.txt', sep='\t')
y_train.to_csv('y_train.txt', sep='\t')
y_test.to_csv('y_test.txt', sep='\t')
X_train = data_normalization(X_train)

# select features by RF and get auc score by SVM
for _ in range(iter_num):
    threshold = 0.9
    st1_data, st1_label = X_train, y_train
    st1_data = feature_selection_rf(st1_data, st1_label, feature_num=feature_num) #select features
    # st1_data_label = pd.concat([st1_data,st1_label],axis=1)
    output = get_cv_score_svm(st1_data, st1_label,c) #get svm score
    print("{} attempt: {}".format(_+1,output))
    # if aucscore>threshold then save
    if output[3]>threshold:
        for name in st1_data.columns.tolist():
            if name in best_col.keys():
                best_col[name]+=1
            else: 
                best_col[name] = 1
with open('./feature_selection/RF_select.txt', 'w') as f:
    for feature,times in best_col.items():
        f.write(str(feature)+' '+str(times)+'\n')

