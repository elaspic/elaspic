import numpy as np
import pylab as pl
from sklearn import ensemble
from sklearn import datasets
from sklearn.utils import shuffle
from scipy.stats.stats import pearsonr
#from sklearn.metrics import mean_squared_error
from sklearn.cross_validation import KFold
from matplotlib.collections import LineCollection
from sklearn.linear_model import LinearRegression
from sklearn.utils import check_random_state
import sys




print "Processing interface mutations"
intr_train_X=np.loadtxt(open("SKEMPINR.X.dat","rb"),delimiter="\t",skiprows=0);
intr_train_Y=np.loadtxt(open("SKEMPINR.Y.dat","rb"),delimiter="\t",skiprows=0);
intr_disease_X=np.loadtxt(open("disease_interface.X.dat","rb"),delimiter="\t",skiprows=0);
maxIter=2000;maxDepth=5;alpha=0.8;learnRate=0.05;minSplitSize=2;minLeafSize=4;
params1 = {'subsample':alpha,'n_estimators':maxIter, 'max_depth': maxDepth, 'learning_rate':learnRate,'min_samples_split': minSplitSize,'min_samples_leaf': minLeafSize, 'loss': 'ls'}
intr_clf = ensemble.GradientBoostingRegressor(**params1)
intr_clf.fit(intr_train_X,intr_train_Y)
intr_disease_Y=intr_clf.predict(intr_disease_X);
np.savetxt("disease_interface.Y.dat", intr_disease_Y)


print "Processing core mutations"
core_train_X=np.loadtxt(open("PRTRMNR.X.dat","rb"),delimiter="\t",skiprows=0);
core_train_Y=np.loadtxt(open("PRTRMNR.Y.dat","rb"),delimiter="\t",skiprows=0);
core_disease_X=np.loadtxt(open("disease_core.X.dat","rb"),delimiter="\t",skiprows=0);
maxIter=1500;maxDepth=5;alpha=0.8;learnRate=0.05;minSplitSize=2;minLeafSize=4;
params1 = {'subsample':alpha,'n_estimators':maxIter, 'max_depth': maxDepth, 'learning_rate':learnRate,'min_samples_split': minSplitSize,'min_samples_leaf': minLeafSize, 'loss': 'ls'}
core_clf = ensemble.GradientBoostingRegressor(**params1)
core_clf.fit(core_train_X,core_train_Y)
core_disease_Y=core_clf.predict(core_disease_X);
np.savetxt("disease_core.Y.dat", core_disease_Y)



