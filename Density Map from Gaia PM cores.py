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

cwd = os.getcwd() 
os.chdir(cwd+"\\CSV files")
file = pd.read_csv( 'table_Gaia-gaia_dr3_source-Box900.csv' )
yso = pd.read_csv('Final file.csv')
notnullfilter= (file['pmra'].notnull() & 
                file['pmdec'].notnull() & 
                file['parallax'].notnull() &
                ((file['pmra_error']**2 + file['pmdec_error']**2) < 1) &
                (file['pm'] < 50) )

source1 = [-5.0617,-2.11039]
source2 = [-3.98761,-2.11039]

filter1 = (((file['pmra']-source1[0])**2 + (file['pmdec']-source1[1])**2)**(0.5) < 1)
filter2 = (((file['pmra']-source2[0])**2 + (file['pmdec']-source2[1])**2)**(0.5) < 1)
file1 = file[notnullfilter & filter1] 
file2 = file[notnullfilter & filter2]
ra1 = np.asarray(file1['ra'])
dec1 = np.asarray(file1['dec'])
ra2 = np.asarray(file2['ra'])
dec2 = np.asarray(file2['dec'])
ysora = np.asarray(yso['ra'])
ysodec = np.asarray(yso['dec'])
plt.scatter(ra1,dec1,s = [5 for n in range(len(dec1))],marker = '*',color = 'red')
plt.scatter(ra2,dec2,s = [5 for n in range(len(dec2))],marker = '+',color = 'green')
plt.scatter(ysora,ysodec,s = [15 for n in range(len(ysodec))],marker = '.',color = 'black')
plt.gca().invert_xaxis()
plt.show()
density,a1,a2 = densityfunction(ra1,dec1,250,10)
density2,a1,a2 = densityfunction(ra2,dec2,250,10)
plt.imshow(density,cmap='afmhot')
plt.colorbar()
plt.show()
plt.imshow(density2,cmap ='afmhot')
plt.colorbar()
plt.show()