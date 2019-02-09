from .utils import Atom, Residue, ActiveSite
import rmsd
import numpy as np
from itertools import permutations
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    site_a_coords = []
    #unpack the ActiveSite class to extract the atom coords of interest
    for i in range(0, len(site_a.residues)): #loop through all residues
        for j in range(0, len(site_a.residues[i].atoms)): #I chose to loop through all atoms instead of the backbone atoms because all atoms improved the shilhouette mean score
            a = site_a.residues[i].atoms[j].coords
            site_a_coords.append(a)
    site_b_coords = []
    #unpack the site_b and save the coords
    for i in range(0, len(site_b.residues)):
        for j in range(0, len(site_b.residues[i].atoms)): 
            b = site_b.residues[i].atoms[j].coords
            site_b_coords.append(b)
    #compute the Root-mean-square deviation
    #low RMSD values (close to 0) are similar structures, compared to higher RMSD values (>>1) are dissimilar
    #RMSD values are never negative
    #Generally, the RMSD calculation is the square root of the standard deviation
    #in atomic structures, the RMSD calculation is the square root of the distance between two atoms with x,y,z coords
    #the package implemented here performs rotations to find the lowest RMSD
    similarity = rmsd.rmsd(site_a_coords, site_b_coords)

    return similarity

def cluster_by_partitioning(active_sites, num_of_clusters):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.
    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    n = num_of_clusters
    rmsd_matrix = np.eye(len(active_sites)) #initialize an empty matrix 
    for i in range(0, len(active_sites)):
        for j in range(0,len(active_sites)):
            rmsd_matrix[i,j] = compute_similarity(active_sites[i], active_sites[j]) #perform RMSD between all PDBs and create matrix
   #K-means imported form sklearn
   #K means randomly selects K clusters and divides the elements into voronoi cells 
   #the next step is to update the cluster centers to be the mean of the centroids of the assigned clusters
   #this is repeated until the labels no longer change. I used the default number of iterations of 10
    kmeans = KMeans(n_clusters=n, random_state=10).fit(rmsd_matrix)
    
    return [rmsd_matrix, kmeans.labels_]

    

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    single_linkage_matrix = np.eye(len(active_sites)) #initialize an empty matrix
    for i in range(0, len(active_sites)):
        for j in range(0, len(active_sites)):
            single_linkage_matrix[i,j] = compute_similarity(active_sites[i], active_sites[j]) #perform RMSD between all PDBs and create matrix
    #Single linkage agglomerative clustering works by assigning every element to its own cluster and finding the nearest neighbor to bottom-up create clusters.
    #The single linkage method uses the closest neighbor of the clusters as the distance metric compared the the farthest or the average member.
    clustering = AgglomerativeClustering(linkage='single').fit(single_linkage_matrix)
    
    return [single_linkage_matrix, clustering.labels_]

def silhouette(similarity_matrix, predicted_labels):
    #silhouette coefficient is between -1 and 1, with higher values indicate a member being very similar to the cluster its in. This function
    #outputs the mean silhouette of the all elements. The silhouette is generally a score of how well a given element is assigned to its cluster and how different it is from the next nearest cluster. 
    MSC = silhouette_score(similarity_matrix, predicted_labels)
    return MSC

def pca(similarity_matrix, predicted_labels):
    n_components = 2
    pca = PCA(n_components = n_components)
    X_pca = pca.fit_transform(similarity_matrix)
    #colors = ["navy", "red"]
    return X_pca

def cluster_comparision(labels):
    lab1 = labels[0] #labels for one clustering method
    lab2 = labels[1] #labels for another clustering method
    counter = 0
    for i in range(0, len(lab1)):
        if lab1[i] == lab2[i]:
            counter += 1
    return counter/len(lab1) #this returns a shared labels percentage.

def unpack(sites): #this just returns the active sites of inputs
    res_list = []
    for site in sites:
        for i in range(0, len(site.residues)):
            res = site.residues
            res_list.append(res)
    return res_list

