import numpy as np 
from matplotlib import pyplot as plt
import pandas as pd
import os

os.chdir("C:\\Users\\Pratham\\Desktop\\IIA Internship\\CSV files")
#gaia = file = pd.read_csv( 'table_irsa_catalog_search_results_gaiaforngc.csv' )
gaia = file = pd.read_csv( 'table_Gaia-gaia_dr3_source-Box900.csv' )
notnullfilter= (file['pmra'].notnull() & 
                file['pmdec'].notnull() & 
                ((file['pmra_error']**2 + file['pmdec_error']**2) < 1 ) )


source1 = [-5.0617,-2.11039]
source2 = [-3.98761,-2.11039]

filter = (((file['pmra']-source1[0])**2 + (file['pmdec']-source1[1])**2)**(0.5) < 1) | (((file['pmra']-source2[0])**2 + (file['pmdec']-source2[1])**2)**(0.5) < 1)
file = file[notnullfilter & filter] 
ra = np.asarray(file['ra'])
dec = np.asarray(file['dec'])
plt.scatter(ra,dec,s = [1 for n in range(len(dec))],marker = '*')
plt.show()