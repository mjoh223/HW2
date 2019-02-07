from hw2skeleton import cluster
from hw2skeleton import io
import os
import random
def test_similarity():
    #set seed here
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "276.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion
    assert cluster.compute_similarity(activesite_a, activesite_b) == 0.0

def test_partition_clustering():
    random.seed(40)
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    label=cluster.cluster_by_partitioning(active_sites, 3)
    #assert cluster.cluster_by_partitioning(active_sites) == []
    assert all(label == [0, 2, 1])

def test_hierarchical_clustering():
    random.seed(40)
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    label=cluster.cluster_hierarchically(active_sites)
    assert all(label == [0, 0, 1])
