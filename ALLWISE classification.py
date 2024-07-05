import numpy as np
import pandas as pd
import scipy
import os
import matplotlib.pyplot as plt

cwd = os.getcwd() 
os.chdir(cwd+"\\CSV files")
file = pd.read_csv( 'table_irsa_catalog_search_results.csv' )
#print(file)
file['w1-w2'] = file['w1mpro']-file['w2mpro']
file['w2-w3'] = file['w2mpro']-file['w3mpro']

'''
#IRAS 18391+0805 source
rasource = (18*15)+ (41*15/60) + (34.8*15/3600)
decsource = (8) + (8/60) + (22/3600)
# 1 star in J-H, H-K,just below CTTS line,1 star just next to dwarf line, 1 Class I source and 1 Class II source in H-K, w1-w2

#dobashi1641 source
rasource = (18*15) + (41*15/60) + (30.4*15/3600)
decsource = (8) + (8/60) + (54/3600)
# 1 star in J-H, H-K,just below CTTS line,2 star just next to dwarf line, 1 Class I source and 1 Class II source in H-K, w1-w2

#leda 62265
rasource = (18*15) + (41*15/60) + (35.25*15/3600)
decsource = (8) + (8/60) + (18.5/3600)
# 1 star in J-H, H-K,just below CTTS line,  1 Class I source and 1 Class II source in H-K, w1-w2

#bd+073807 source
rasource = (18*15) + (41*15/60) + (47.7206*15/3600)
decsource = (8) + (7/60) + (1.115/3600)
#6 stars in J-H, H-K, 1 Class I source in H-K, w1-w2

#IRAS 18393+0759 source
rasource = (18*15) + (41*15/60) + (44.3*15/3600)
decsource = (8) + (2/60) + (16/3600)
# 2 star in J-H, H-K, 1 Class I source in H-K, w1-w2

#2Mass J18414362+0811365 source
rasource = (18*15) + (41*15/60) + (43.62*15/3600)
decsource = (8) + (11/60) + (36.5/3600) 
# 1 star in J-H, H-K, 1 Class I source in H-K, w1-w2

#BD+073811
rasource = (18*15) + (42*15/60) + (6.47*15/3600)
decsource = (8) + (4/60) + (1/3600) 
# No source anywhere
'''

#IRAS 18391+0805 source
rasource = (18*15)+ (41*15/60) + (34.8*15/3600)
decsource = (8) + (8/60) + (22/3600)

loc_filter = (np.abs(file['ra']-rasource) <=(1/60))&(np.abs(file['dec']-decsource) <= (1/60))

Decontamination_Filter = ((file["w1rchi2"]<(file["w1snr"]-3)/7) &
                          (file["w2rchi2"] < 0.1 * file["w2snr"]-0.3) & 
                          ((file['w3snr']>=5) & ((0.45<file['w3rchi2']) & ((file['w3rchi2']<1.15) | (file["w3rchi2"] < 0.125 * file["w3snr"]-1)))) & 
                          (file['w4rchi2'] < 0.2* file['w4snr']-2))

Real_Data = file[Decontamination_Filter]
print(Real_Data)

ax = Real_Data.plot(kind = 'scatter', y = 'w1-w2', x = 'w2-w3')
'''xpoints = [2,2]
ypoints = [10,-10]
xpoints2= [10,-10]
ypoints2=[-2,6.4]
ypoints3=[3.7,-5.5]
ypoints4 = [0.25,0.25]
ypoints5 = [8.75,-10.25]
ypoints6 = [-12.9,17.1]
plt.plot(xpoints,ypoints,linestyle ='dashed')
plt.plot(xpoints2,ypoints2,linestyle = 'dashed')
plt.plot(xpoints2,ypoints3,linestyle = 'dashed')
plt.plot(xpoints2,ypoints4,linestyle = 'dashed')
plt.plot(xpoints2,ypoints5,linestyle = 'dashed')
plt.plot(xpoints2,ypoints6,linestyle = 'dashed')'''
y = z = np.linspace(-1, 5, 200)
x = 0*y+2
plt.plot(x, y, color = 'red')
x = 0*y+4.5
plt.plot(x, y,color = 'green')
x = (-y+2.2)/0.42
plt.plot(x, y, color = 'red')
y = 0.46*(x)-0.9
plt.plot(x, y, color = 'red')
y = 0.25+0*z
plt.plot(z,y, color = 'red')
y = 0.9 * x -0.25
plt.plot(x, y,color = 'green')
y = -1.5 * x + 2.1
plt.plot(x, y, color = 'red')
plt.xlim(-1,5)
plt.ylim(-1,3)
plt.show()

