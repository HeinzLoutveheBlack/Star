import numpy as np
import pandas as pd # type: ignore
import scipy # type: ignore
import os
from matplotlib import pyplot as plt # type: ignore


def closestmatch(array1,array2,array3,array4,array5):
    # The variables are arrays in order of : Source RA, Source DEC, Destination RA, Destination DEC, Destination PMRA, Destination PMDEC.
    # The Destination is mainly Gaia csv file

    dist = []
    for i in range(len(array1)):
        ra_source = array1[i]
        dec_source = array2[i]
        for j in range(len(array3)):
            ra_destinaition = array3[j]
            dec_destination = array4[j]
            dist.append((((ra_source-ra_destinaition)**2)+(dec_source-dec_destination)**2))

    dist_array = np.asarray(dist)
    dist_array_reshaped = dist_array.reshape(len(array1),len(array3))
    indices = np.argmin(dist_array_reshaped, axis = 1)
    candidate_distance = []
    for i in range(len(indices)):
        append_index = indices[i]
        candidate_distance.append(array5[append_index])

    return candidate_distance
    # The return order is declination parameters followed by right ascension parameters

os.chdir("C:\\Users\\Pratham\\Desktop\\IIA Internship\\CSV files")
file = pd.read_csv( 'Final File.csv' )
gaia = pd.read_csv( 'table_Gaia-gaia_dr3_source-Box900.csv' )

gaiafilter = (gaia['pmra'].notnull() & gaia['pmdec'].notnull() & gaia['parallax'].notnull())
gaia = gaia[gaiafilter]
gaiara = np.asarray(gaia['ra'])
gaiadec = np.asarray(gaia['dec'])
filera = np.asarray(file['ra'])
filedec = np.asarray(file['dec'])
gaiadist = np.asarray(gaia['distance_gspphot'])

distarray = closestmatch(filera,filedec,gaiara,gaiadec,gaiadist)
print(file['designation'])
print(distarray)


