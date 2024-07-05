import numpy as np 
from matplotlib import pyplot as plt
import pandas as pd
import os
from scipy.spatial import distance,KDTree
import numpy as np

def grid(x_min, x_max, y_min, y_max, num_points):
    x = np.linspace(x_min, x_max, num_points)
    y = np.linspace(y_min, y_max, num_points)
    xx, yy = np.meshgrid(x, y)
    xxflat = xx.ravel()
    yyflat = yy.ravel()
    grid_points = []
    for i in range(len(xxflat)):
        grid_points.append([xxflat[i],yyflat[i]])
    return(grid_points)

def densityfunction(array1,array2,precision,k):
    max1 = np.max(array1) 
    min1 = np.min(array1) 
    max2 = np.max(array2) 
    min2 = np.min(array2) 
    
    grid_array = grid(min1,max1,min2,max2,precision)
    np.reshape(grid_array,(2*precision,precision))
    #print(grid_array)
    data = []
    for i in range(len(array1)):
        data.append((array1[i],array2[i]))
    #print(data)
    tree = KDTree(data)

    # 5th Nearest Neighbor Search
    density = []
    for point in grid_array:
        distances, indices = tree.query(point, k)
        density.append(k/(3600 * np.pi * distances[-1]**2)) # density is in stars/arcsec^2
    demsity = np.reshape(density,(precision,precision))
    renormalised1 = []
    renormalised2 = []
    for i in range(len(array1)):
        renormalised1.append(precision * (array1[i] - min1 ) /(max1-min1))
        renormalised2.append(precision * (array2[i] - min2 ) /(max2-min2))
    
    return demsity,renormalised1,renormalised2
os.chdir("C:\\Users\\Pratham\\Desktop\\IIA Internship\\CSV files")
'''file = pd.read_csv( 'table_irsa_catalog_search_results_gaiaforngc.csv' )'''
file = pd.read_csv( 'table_Gaia-gaia_dr3_source-Box900.csv' )
yso = pd.read_csv('Final file.csv')
notnullfilter= (file['pmra'].notnull() & 
                file['pmdec'].notnull() & 
                file['parallax'].notnull() &
                ((file['pmra_error']**2 + file['pmdec_error']**2) < 1) &
                (file['pm'] < 50))
file = file[notnullfilter] 
coordinates_ra_array = np.asarray(file['ra'])
coordinates_dec_array = np.asarray(file['dec'])
yso_ra = np.asarray(yso['ra'])
yso_dec = np.asarray(yso['dec'])
density,a1,a2 = densityfunction(coordinates_ra_array,coordinates_dec_array,250,20)
plt.imshow(density,cmap='viridis')
plt.colorbar()
plt.show()

pmdecarray = np.asarray(file['pmdec'])
pmraarray = np.asarray(file['pmra'])
yso_pmra = np.asarray(yso['pmra'])
yso_pmdec = np.asarray(yso['pmdec'])
pmdensity,aa,ab = densityfunction(pmraarray,pmdecarray,500,80)
plt.imshow(pmdensity,cmap = 'afmhot')
#plt.scatter(aa,ab,color = 'white',marker = 'o',s = [0.2 for i in range(len(aa))])
plt.colorbar()
plt.show()
print(np.min(pmraarray)+((np.max(pmraarray)-np.min(pmraarray))*242.5/500))
print(np.min(pmdecarray)+((np.max(pmdecarray)-np.min(pmdecarray))*275/500))
print(np.min(pmraarray)+((np.max(pmraarray)-np.min(pmraarray))*242.5/500))
print(np.min(pmdecarray)+((np.max(pmdecarray)-np.min(pmdecarray))*282.5/500))
print(yso_pmra)
print(yso_pmdec)

