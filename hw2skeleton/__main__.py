import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import compute_similarity, cluster_by_partitioning, cluster_hierarchically, silhouette, pca

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])
#silhouette score for multiple clusters
for i in range(2, 2):
    kmeans = cluster_by_partitioning(active_sites, i)
    similarity_matrix, labels  = kmeans
    msc = silhouette(similarity_matrix, labels)
    format_output = [i, msc]
    print("for {} clusters, the mean silhouette score is {} ".format(*format_output))

#plot pca
kmeans = cluster_by_partitioning(active_sites, 2)
similarity_matrix, labels  = kmeans
print(pca(similarity_matrix, labels))



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

