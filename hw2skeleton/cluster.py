from .utils import Atom, Residue, ActiveSite
import rmsd
import numpy as np
from itertools import permutations
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    site_a_coords = []
    for i in range(0, len(site_a.residues)):
        for j in range(0, 2): #only look at the first 3 atoms which correspond to the backbone
            a = site_a.residues[i].atoms[j].coords
            site_a_coords.append(a)
    site_b_coords = []
    for i in range(0, len(site_b.residues)):
        for j in range(0, 2): 
            b = site_b.residues[i].atoms[j].coords
            site_b_coords.append(b)

    similarity = rmsd.rmsd(site_a_coords, site_b_coords)

    return similarity

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.
    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    rmsd_matrix = np.eye(len(active_sites)) 
    for i in range(0, len(active_sites)):
        for j in range(0, len(active_sites)):
            rmsd_matrix[i,j] = compute_similarity(active_sites[i], active_sites[j])
    kmeans = KMeans(n_clusters=5, random_state=0).fit(rmsd_matrix)
    
    return list(zip(active_sites, kmeans.labels_))

    

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    single_linkage_matrix = np.eye(len(active_sites)) 
    for i in range(0, len(active_sites)):
        for j in range(0, len(active_sites)):
            single_linkage_matrix[i,j] = compute_similarity(active_sites[i], active_sites[j])
    clustering = AgglomerativeClustering(linkage='single').fit(single_linkage_matrix)
    
    return list(zip(active_sites, clustering.labels_))
