import numpy as np
import pandas as pd # type: ignore
import scipy # type: ignore
import os
from matplotlib import pyplot as plt # type: ignore
import math
from mpl_toolkits import mplot3d
import networkx as nx #type:ignore 
from scipy.spatial.distance import pdist, squareform

def propermotion(array1,array2,array3,array4,array5,array6):
    # The variables are arrays in order of : Source RA, Source DEC, Destination RA, Destination DEC, Destination PMRA, Destination PMDEC.
    # The Destination is mainly Gaia csv file

    dist = []
    for i in range(len(array1)):
        ra_source = array1[i]
        dec_source = array2[i]
        for j in range(len(array3)):
            ra_destination = array3[j]
            dec_destination = array4[j]
            dist.append((((ra_source-ra_destination)**2)+(dec_source-dec_destination)**2))
            # This is to find the distance between one point of the source array to the points of the distance array

    dist_array = np.asarray(dist)
    dist_array_reshaped = dist_array.reshape(len(array1),len(array3)) #reshaping the array in the shape of source arrays, so that for each point 'i', we can find the minimum along each row.
    indices = np.argmin(dist_array_reshaped, axis = 1)
    candidate_pmra = []
    candidate_pmdec = []
    for i in range(len(indices)):
        append_index = indices[i]
        candidate_pmra.append(array5[append_index])
        candidate_pmdec.append(array6[append_index])
    
    return candidate_pmdec, candidate_pmra
    # The return order is declination parameters followed by right ascension parameters

def closestmatch(array1,array2,array3,array4):
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
    candidate_ra = []
    candidate_dec = []
    for i in range(len(indices)):
        append_index = indices[i]
        candidate_ra.append(array3[append_index])
        candidate_dec.append(array4[append_index])
    
    return candidate_dec,candidate_ra
    # The return order is declination parameters followed by right ascension parameters

def coordinatetodistance(array1,array2,array3,array4,array5,array6,cap):
# The function appends the distance of the corresponding closest point in gaia to a given point in any other csv file to the return array, if the distance is positive, not nan or lesser than a cap value
    dist = []
    for i in range(len(array1)):
        ra_source = array1[i]
        dec_source = array2[i]
        for j in range(len(array3)):
            ra_destination = array3[j]
            dec_destination = array4[j]
            dist.append((((ra_source-ra_destination)**2)+(dec_source-dec_destination)**2))

    dist_array = np.asarray(dist)
    dist_array_reshaped = dist_array.reshape(len(array1),len(array3))
    indices = np.argmin(dist_array_reshaped, axis = 1)
    distance =[]
    for i in range(len(indices)):
        append_index = indices[i]
        if np.isnan(array5[append_index]):
            if np.isnan(array6[append_index]):
                distance.append(0)
            elif array6[append_index]>cap:
                distance.append(0)
            elif array6[append_index]<0:
                distance.append(0)
            else: 
                distance.append(array6[append_index])
        else:
            distance.append(array5[append_index])
    return distance
    # The return order is declination parameters followed by right ascension parameters

def kruskal_mst(graph):
    mst = nx.Graph()
    mst.add_nodes_from(graph.nodes)
    
    edges = sorted(graph.edges(data=True), key=lambda x: x[2]['weight'])
    
    disjoint_set = {node: node for node in graph.nodes}
    
    def find(node):
        if disjoint_set[node] != node:
            disjoint_set[node] = find(disjoint_set[node])
        return disjoint_set[node]
    
    def union(node1, node2):
        root1 = find(node1)
        root2 = find(node2)
        disjoint_set[root1] = root2
    
    for u, v, data in edges:
        if find(u) != find(v):
            mst.add_edge(u, v, weight=data['weight'])
            union(u, v)
    
    return mst


cwd = os.getcwd() 
os.chdir(cwd+"\\CSV files")
'''wise = pd.read_csv( 'table_WISE-allwise_p3as_psd-Cone1200_ngc.csv' )
gaia = pd.read_csv( 'table_irsa_catalog_search_results_gaiaforngc.csv' )'''

wise = file = pd.read_csv( 'table_irsa_catalog_search_results.csv' )
gaia = pd.read_csv( 'table_Gaia-gaia_dr3_source-Box900.csv' )

#akari = pd.read_csv('table_AKARI-akari_irc-Box900.csv')
#print(wise)

wise['w1-w2'] = wise['w1mpro']-wise['w2mpro']

'''
wise['H-K'] = wise['h_m_2mass']-wise['k_m_2mass']
wise['J-H'] = wise['j_m_2mass']-wise['k_m_2mass']
'''

wise['H-Ks (2Mass)'] = wise['h_m_2mass']-wise['k_m_2mass']
wise['J-H (2Mass)'] = wise['j_m_2mass']-wise['h_m_2mass']

wise['H-K'] = (wise['H-Ks (2Mass)']-0.028)/0.996 # 2Mass Color to Bessel and Brett colors
wise['J-H'] = (wise['J-H (2Mass)']+0.045)/0.980 # From Carpenter et al 2001


ldn_637_location = [280.25,8.10]
#IRAS18391+0805 source
source_iras18391_0805 = [(18*15)+ (41*15/60) + (34.8*15/3600),(8) + (8/60) + (22/3600) ] #= 280.395 #= 8.13944 
# 1 star in J-H, H-K,just below CTTS line,1 star just next to dwarf line, 1 Class I source and 1 Class II source in H-K, w1-w2

#dobashi1641 source
source_dobashi1641 = [(18*15) + (41*15/60) + (30.4*15/3600),(8) + (8/60) + (54/3600)] #= 380.376  #= 8.1483
# 1 star in J-H, H-K,just below CTTS line,2 star just next to dwarf line, 1 Class I source and 1 Class II source in H-K, w1-w2

#leda_62265
source_leda_62265 = [(18*15) + (41*15/60) + (35.25*15/3600),(8) + (8/60) + (18.5/3600)] #= 280.396875  #= 8.138472
# 1 star in J-H, H-K,just below CTTS line,  1 Class I source and 1 Class II source in H-K, w1-w2

#bd+073807 source
source_bd_073807 = [(18*15) + (41*15/60) + (47.7206*15/3600),(8) + (7/60) + (1.115/3600)] #= 280.44883583 #= 8.116976
#6 stars in J-H, H-K, 1 Class I source in H-K, w1-w2

#IRAS_18393+0759 source
source_IRAS_18393_0759 = [(18*15) + (41*15/60) + (44.3*15/3600),(8) + (2/60) + (16/3600)] #= 280.434583  #= 8.03777
# 2 star in J-H, H-K, 1 Class I source in H-K, w1-w2

#2Mass_J18414362+0811365 source
source_2Mass_J18414362_0811365 = [(18*15) + (41*15/60) + (43.62*15/3600), (8) + (11/60) + (36.5/3600)] #= 280.43175 #= 8.193472
# 1 star in J-H, H-K, 1 Class I source in H-K, w1-w2

#BD+073811
source_BD_073811 = [(18*15) + (42*15/60) + (6.47*15/3600),(8) + (4/60) + (1/3600)] #= 280.526958#= 8.066944
# No source anywhere

source_g03909_0585 = [(18*15) + (41*15/60) + (41.177*15/3600),(8) + (8/60) + (10.1/3600)]

#Tycho Sources
source_J18414110_0805039 =[280.42124056,8.08453139,760.3026]
source_J18412261_0810424 =[280.34434944,8.17850500,976.6436]
source_J18415588_0805369 =[280.48288333,8.09358333,0] ##
source_J18411116_0811122 =[280.29656750,8.18691750,0] ####
source_J18410587_0810413 =[280.27446944,8.17827028,285.1608]
source_J18410123_0809452 =[280.25520056,8.16264278,0] ## #Tycho-1025-2122-1
source_J18413373_0758356 =[280.39051500,7.97659278,758.5649] 

sourcearray =[source_IRAS_18393_0759, source_2Mass_J18414362_0811365, source_bd_073807,source_BD_073811,source_dobashi1641, source_iras18391_0805, source_leda_62265, source_g03909_0585]
tychosourcearray= [source_J18414110_0805039,source_J18412261_0810424,source_J18415588_0805369,source_J18411116_0811122,source_J18410587_0810413,source_J18410123_0809452,source_J18413373_0758356]
'''
loc_filter = (np.abs(wise['ra']-rasource) <=(2/60)) & (np.abs(wise['dec']-decsource) <= (2/60))
'''

not_null_filter = ((wise['h_m_2mass'].notnull()) & 
                  (wise['k_m_2mass'].notnull()) &
                  (wise['j_m_2mass'].notnull()) )

Decontamination_Filter = (not_null_filter &
                         (wise["w1rchi2"] < (wise["w1snr"]-3) / 7 ) & 
                         (wise["w2rchi2"] < 0.1 * wise["w2snr"] - 0.3 ))

wise_2mass_filter = ((wise['H-K'] > 0) &
                     (wise['H-K'] < 3.4375 * (wise['w1-w2']) - 0.85) &
                     (wise['H-K'] > -1.76 * (wise['w1-w2']) + 0.9))

Decontamination_Filter_2 = ((file["w1rchi2"]<(file["w1snr"]-3)/7) &
                          (file["w2rchi2"] < 0.1 * file["w2snr"]-0.3) & 
                          ((file['w3snr']>=5) & ((0.45<file['w3rchi2']) & ((file['w3rchi2']<1.15) | (file["w3rchi2"] < 0.125 * file["w3snr"]-1)))) & 
                          (file['w4rchi2'] < 0.2* file['w4snr']-2))

Real_Data_2 = file[Decontamination_Filter_2]
Real_Data = wise[Decontamination_Filter & wise_2mass_filter] 
print(Real_Data)
Real_Data.plot(kind = 'scatter', y = 'H-K', x = 'w1-w2',color = 'red')
plt.scatter(0.578,(0.568-0.028)/0.996,color = 'black')
#plt.scatter(0.152,(-0.059-0.028)/0.996,color = 'black') #Tycho-1025-2122-1
#plt.scatter(0.129,(-0.071-0.028)/0.996,color= 'black') #Tycho-1025-2578-1
xpoints= [10,-10]
ypoints = [0,0] # y= 0
ypoints2=[-16.7,18.5] # y = -1.76 x + 0.9
ypoints3=[33.525,-35.225] #y = 3.4375 x - 0.85
ypoints4 = [-15.05,20.15] # y = -1.76 x + 2.55
plt.plot(xpoints,ypoints,linestyle ='dashed',color = 'red')
plt.plot(xpoints,ypoints2,linestyle = 'dashed',color = 'green')
plt.plot(xpoints,ypoints3,linestyle = 'dashed',color = 'blue')
plt.plot(xpoints,ypoints4,linestyle = 'dashed',color = 'orange')
plt.xlim((-1,3))
plt.ylim((-1,4))
plt.show()

yso_filter = (((wise['J-H'] +0.043)/(1.076) > 0.58 * (wise['H-K']-0.028)/(1.026) + 0.52) &
              (wise['J-H'] > 1.6928 * wise['H-K'] - 0.5928) & 
              (wise['J-H'] < 1.6928 * wise['H-K'] - 0.033664) )

wise_yso = wise[not_null_filter & yso_filter]
wise_yso.plot(kind = 'scatter', y = 'J-H', x = 'H-K',color = 'red')
plt.scatter((0.568+0.045)/0.980,(1.04-0.028)/0.996,color='Black')
#plt.scatter((0.152+0.045)/0.980,(0.407-0.028)/0.996,color='Black') #Tycho-1025-2122-1 # A6 type dwarf, very close to being a G0 type Giant
#plt.scatter((0.129+0.045)/0.980,(0.564-0.028)/0.996,color='Black') #Tycho-1025-2578-1 # G4 type dwarf
# Dwarf sequence
plt.plot([-0.035,0.00,0.005,0.015,0.025,0.03,0.035,0.04,0.045,0.05,0.052,0.055,0.06,0.075,0.09,0.105,0.11,0.13,0.165,0.20,0.21,0.25,0.275,0.32,0.37],[-0.05,0.00,0.02,0.06,0.09,0.13,0.165,0.23,0.285,0.305,0.32,0.33,0.37,0.45,0.50,0.58,0.61,0.66,0.695,0.68,0.665,0.62,0.60,0.62,0.66], linestyle = 'dashdot',color = 'blue',ms = 20)
# Giant sequence
plt.plot([0.065,0.08,0.085,0.085,0.095,0.10,0.115,0.14,0.15,0.165,0.19,0.205,0.215,0.235],[0.37,0.47,0.50,0.50,0.54,0.58,0.63,0.68,0.73,0.79,0.83,0.85,0.87,0.90],linestyle= 'dashdot',color = 'orange', ms= 20)
plt.plot([0.235,2.5],[0.9,4.745825],linestyle = 'dashed', color = 'black', ms = 20)
plt.plot([0.37,2.5],[0.66,4.277425],linestyle = 'dashed', color = 'black', ms = 20)
plt.plot([1,(0.4-0.028)/1.026],[1.1,(0.752+0.043)/1.076],linestyle = 'dotted', color = 'green', ms = 20) #CTTS Line, colors changed to 2mass
'''plt.plot([1,0.4],[1.1,0.752],linestyle = 'dotted', color = 'red', ms = 20) #CTTS Line'''
plt.plot([1,2.5],[1.1,3.647495],linestyle = 'dashed', color = 'black', ms = 20)
plt.show()

Combined_Yso_array = [wise_yso,Real_Data,Real_Data_2]
Combined_Yso_DF = pd.concat(Combined_Yso_array)
print(Combined_Yso_DF)
ra_array = np.asarray(Combined_Yso_DF['ra'])
#print(len(ra_array))
dec_array = np.asarray(Combined_Yso_DF['dec'])
gaia_ra_array = np.asarray(gaia['ra'])
#print(len(gaia_ra_array))
gaia_dec_array = np.asarray(gaia['dec'])
'''akari_ra_array = np.asarray(akari['ra'])
akari_dec_array = np.asarray(akari['dec'])'''
gaia_pmra_array = np.asarray(gaia['pmra'])
gaia_pmdec_array = np.asarray(gaia['pmdec'])
gaia_distance = np.asarray(gaia['distance_gspphot'])
gaia_distance_parallax = np.asarray(1000/gaia['parallax'])

'''e,f = closestmatch(akari_ra_array,akari_dec_array,gaia_ra_array,gaia_dec_array)
check_dec,check_ra = closestmatch(akari_ra_array,akari_dec_array,ra_array,dec_array)'''
gaia_ysodec,gaia_ysora = closestmatch(ra_array,dec_array,gaia_ra_array,gaia_dec_array)
wisepmdec,wisepmra = propermotion(ra_array,dec_array,gaia_ra_array,gaia_dec_array,gaia_pmra_array,gaia_pmdec_array)
#akaripmdec, akaripmra = propermotion(check_ra,check_dec,gaia_ra_array,gaia_dec_array,gaia_pmra_array,gaia_pmdec_array)
ysodistance = coordinatetodistance(gaia_ysora,gaia_ysodec,gaia_ysora,gaia_ysodec,gaia_distance,gaia_distance_parallax,2000)
#print(ysodistance)

gaia.plot(kind = 'scatter', y ='pmra', x = 'pmdec', marker="*", s=1)
plt.scatter(wisepmdec,wisepmra,color = 'red', s= [2 for n in range(len(wisepmdec))])
plt.scatter(-5.359, -3.415, color = 'Black',s = [1])
source = [-5.359, -3.415]
size_of_box = (2)
plt.plot([source[0]+(size_of_box),source[0]+(size_of_box)],[source[1]-(size_of_box),source[1]+(size_of_box)],linestyle = 'dotted',color = 'red')
plt.plot([source[0]-(size_of_box),source[0]-(size_of_box)],[source[1]-(size_of_box),source[1]+(size_of_box)],linestyle = 'dotted',color = 'red')
plt.plot([source[0]+(size_of_box),source[0]-(size_of_box)],[source[1]-(size_of_box),source[1]-(size_of_box)],linestyle = 'dotted',color = 'red')
plt.plot([source[0]+(size_of_box),source[0]-(size_of_box)],[source[1]+(size_of_box),source[1]+(size_of_box)],linestyle = 'dotted',color= 'red')
plt.show()
'''
plt.scatter(akaripmdec,akaripmra, color = 'black', s = [5 for n in range(len(akaripmdec))])
plt.show()
plt.scatter(e,f,color = 'black', s= [5 for n in range(len(f))])'''


plt.scatter(ra_array, dec_array, color = 'red', s= [20 for n in range(len(ra_array))])
'''plt.title('RA vs DEC for YSOs with sources of interest')
plt.scatter(ldn_637_location[0], ldn_637_location[1], marker = 's', color = 'black', s= [5])
plt.scatter(source_2Mass_J18414362_0811365[0], source_2Mass_J18414362_0811365[1],marker = 'v', color = 'black', s= [5])
plt.scatter(source_bd_073807[0], source_bd_073807[1], marker = '*', color = 'blue', s= [5])
plt.scatter(source_BD_073811[0], source_BD_073811[1], marker = 'P', color = 'cyan', s= [5])
plt.scatter(source_dobashi1641[0], source_dobashi1641[1], marker = 'H', color = 'green', s= [5])
plt.scatter(source_iras18391_0805[0], source_iras18391_0805[1],marker = '+', color = 'pink', s= [5])
plt.scatter(source_IRAS_18393_0759[0], source_IRAS_18393_0759[1], marker = 'D', color = 'orange', s= [5])
plt.scatter(source_leda_62265[0], source_leda_62265[1], marker = 'X', color = 'purple', s= [5])
plt.scatter(source_g03909_0585[0], source_g03909_0585[1], marker = '.', color = 'black', s= [5])
plt.legend(['YSOs','LDN 637','2MASS J18414362+0811365','BD+073807','BD+073811','DOBASHI 1641','IRAS 18391+0805','IRAS 18393+0759','LEDA 62265','Gsource'],bbox_to_anchor=(1,1), ncol =1 )'''
plt.gca().invert_xaxis()
plt.show()

fig = plt.figure()
ax = plt.axes(projection ='3d')
ax.scatter(ra_array, dec_array, ysodistance, color = 'red', s= [20 for n in range(len(ra_array))])
ax.set_title('RA vs DEC for YSOs with sources of interest')
'''ax.scatter(ldn_637_location[0], ldn_637_location[1],0, marker = 's', color = 'black', s= [5])
ax.scatter(source_2Mass_J18414362_0811365[0], source_2Mass_J18414362_0811365[1],0,marker = 'v', color = 'black', s= [5])
ax.scatter(source_bd_073807[0], source_bd_073807[1],527.8901, marker = '*', color = 'blue', s= [5])
ax.scatter(source_BD_073811[0], source_BD_073811[1],550.2355, marker = 'P', color = 'cyan', s= [5])
ax.scatter(source_dobashi1641[0], source_dobashi1641[1],0, marker = 'H', color = 'green', s= [5])
ax.scatter(source_iras18391_0805[0], source_iras18391_0805[1],0, marker = '+', color = 'pink', s= [5])
ax.scatter(source_IRAS_18393_0759[0], source_IRAS_18393_0759[1],0, marker = 'D', color = 'orange', s= [5])
ax.scatter(source_leda_62265[0], source_leda_62265[1],0, marker = 'X', color = 'purple', s= [5])
ax.scatter(source_g03909_0585[0], source_g03909_0585[1],0, marker = '.', color = 'black', s= [5])

for source in tychosourcearray:
    ax.scatter(source[0],source[1],source[2],color = 'green')
'''

'''size_of_window = (4/60)
size_of_box = size_of_window/2
for source in sourcearray:
    plt.plot([source[0]+(size_of_box),source[0]+(size_of_box)],[source[1]-(size_of_box),source[1]+(size_of_box)],linestyle = 'dotted',color = 'red')
    plt.plot([source[0]-(size_of_box),source[0]-(size_of_box)],[source[1]-(size_of_box),source[1]+(size_of_box)],linestyle = 'dotted',color = 'red')
    plt.plot([source[0]+(size_of_box),source[0]-(size_of_box)],[source[1]-(size_of_box),source[1]-(size_of_box)],linestyle = 'dotted',color = 'red')
    plt.plot([source[0]+(size_of_box),source[0]-(size_of_box)],[source[1]+(size_of_box),source[1]+(size_of_box)],linestyle = 'dotted',color= 'red')
'''
plt.legend(['YSOs','LDN 637','2MASS J18414362+0811365','BD+073807','BD+073811','DOBASHI 1641','IRAS 18391+0805','IRAS 18393+0759','LEDA 62265','Gsource'],bbox_to_anchor=(1,1), ncol =1 )
plt.gca().invert_xaxis()
plt.show()

'''plt.hist(ysodistance,(0,10000))
plt.show()'''

#Real_Data.to_csv('RealData.csv')
abc = Combined_Yso_DF[['ra','dec']].values
'''G = nx.Graph()
for source in abc:     
    for dest in abc:
        G.add_edges_from([(source,dest,{"weight": 1})])
T = nx.minimum_spanning_tree(G)
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, node_color="lightblue", node_size=5)
nx.draw_networkx_edges(T, pos, edge_color="green", width=2)
plt.axis("off")
plt.show()'''

points = abc
# Plot the scatter plot
for source in abc:
    plt.scatter(source[0],source[1],s=[1])
plt.title('Scatter Plot Points')
plt.show()

# Compute the distance matrix
dist_matrix = squareform(pdist(points, 'euclidean'))

# Create a graph from the distance matrix
graph = nx.Graph()
for i in range(len(points)):
    for j in range(i + 1, len(points)):
        graph.add_edge(i, j, weight=dist_matrix[i, j])

mst = kruskal_mst(graph)

# Plot the scatter plot points
plt.scatter(points[:, 0], points[:, 1], color='blue')

distances =[]
# Plot the MST edges
for u, v, data in mst.edges(data=True):
    plt.plot([points[u, 0], points[v, 0]], [points[u, 1], points[v, 1]], color='red')
    distances.append(data['weight'])
plt.title('Minimal Spanning Tree on Scatter Plot')
plt.show()
distances.sort()
print(np.mean(distances))
plt.hist(distances)
plt.show()

#Combined_Yso_DF.to_csv('Final File.csv')