import numpy as np 
from matplotlib import pyplot as plt
import pandas as pd
import os
import seaborn as sns 
from sklearn.cluster import DBSCAN 
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
from sklearn import datasets
import networkx as nx #type:ignore 

def dodbscan(file, eps, min_samples):
    """
    Perform DBSCAN clustering on a dataset and return the number of clusters.
    
    Parameters:
    - file_path: str, path to the CSV file containing the dataset.
    - eps: float, the maximum distance between two samples for one to be considered as in the neighborhood of the other.
    - min_samples: int, the number of samples (or total weight) in a neighborhood for a point to be considered as a core point.

    Returns:
    - n_clusters: int, the number of clusters found by DBSCAN (excluding noise).
    """
    # Load the dataset
    X = file[['ra', 'dec']].values

    # Fit the DBSCAN model
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    
    return n_clusters


cwd = os.getcwd() 
os.chdir(cwd+"\\CSV files")
gaia = file = pd.read_csv( 'table_irsa_catalog_search_results_gaiaforngc.csv' )
'''gaia = file = pd.read_csv( 'table_Gaia-gaia_dr3_source-Box900.csv' )'''
notnullfilter= (file['pmra'].notnull() & 
                file['pmdec'].notnull() & 
                ((file['pmra_error']**2 + file['pmdec_error']**2) < 1 ) )
file = file[notnullfilter] 

clusternumbereps = []
for i in range(10):
    clusternumbereps.append(dodbscan(file,0.00386*(i+1),100))
maxclusternumber = np.max(clusternumbereps)
clusternumberminpts = []
for i in range(10):
    clusternumberminpts.append(dodbscan(file,0.00386 * maxclusternumber,10*(i+1)))
print(clusternumberminpts)

X = file[['ra','dec']].values

db = DBSCAN(eps=0.00386*7, min_samples= 100).fit(X) #np.where(np.max(clusternumberminpts))
labels = db.labels_

# Plot the results
unique_labels = set(labels)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = (labels == k)

    xy = X[class_member_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=6)

plt.title('DBSCAN clustering')
plt.show()