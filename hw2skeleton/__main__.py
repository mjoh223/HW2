import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import compute_similarity, cluster_by_partitioning, cluster_hierarchically, silhouette, pca
from sklearn.decomposition import PCA
import numpy as np
# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])
#silhouette score for multiple clusters
msc_list = np.zeros(len(range(2,10)))
for ii, i in enumerate(range(2, 10)):
    kmeans = cluster_by_partitioning(active_sites, i)
    similarity_matrix, labels  = kmeans
    msc = silhouette(similarity_matrix, labels)
    msc_list[ii] = msc
np.savetxt("msc.csv", msc_list)

single_link = cluster_hierarchically(active_sites)
sim_matrix, labels_single = single_link

#plot pca
n_components = 3
pca = PCA(n_components = n_components)
X_pca = pca.fit_transform(similarity_matrix) 
X_pca = np.array(X_pca)
labels = np.array(labels).astype(int)
stack = np.column_stack((X_pca, labels, labels_single))
#print(stack)
#np.savetxt("pca_clustering.csv", stack, delimiter=",")
#sns.scatterplot(x = X_pca[:,0], y = X_pca[:,1])



# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites, 2)
    #write_clustering(sys.argv[3], clustering)
    #print(clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites)
    #write_mult_clusterings(sys.argv[3], clusterings)
    #print(clusterings)

