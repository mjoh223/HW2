from .utils import Atom, Residue, ActiveSite

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    site_a.Residue(type, 1)
    site_b.Residue(type, 2)
    similarity = 0.0

    # Fill in your code here!
    return [site_a, site_b]
    #return similarity


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!
    list_of_pairs = []
    for i in range( 0, len(active_sites) ):
        for j in range( i, len(active_sites) ):
            pair = [ active_sites[i], active_sites[j] ]
            site_a = ActiveSite( str(pair[0]) )
            site_b = ActiveSite( str(pair[1]) )
            
    return [site_a.name, Residue]


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
