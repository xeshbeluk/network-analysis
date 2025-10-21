
!pip install python-igraph
import igraph
import pandas as pd
import numpy as np
import collections
import pprint
import operator
from typing import List, Set

!curl https://csx46.s3-us-west-2.amazonaws.com/PathwayCommons9.All.hgnc.sif.gz --output PathwayCommons9.All.hgnc.sif.gz
!gunzip -f PathwayCommons9.All.hgnc.sif.gz

sif_data = pd.read_csv("PathwayCommons9.All.hgnc.sif",
                       sep="\t", names=["species1", "interaction_type", "species2"])

interaction_types_ppi = set(["interacts-with",
                             "in-complex-with"])
interac_ppi = sif_data[sif_data.interaction_type.isin(interaction_types_ppi)].copy()

boolean_vec = interac_ppi['species1'] > interac_ppi['species2']
interac_ppi.loc[boolean_vec, ['species1', 'species2']] = interac_ppi.loc[boolean_vec, ['species2', 'species1']].values

interac_ppi_unique = interac_ppi[["species1", "species2"]].drop_duplicates()


ppi_igraph = igraph.Graph.TupleList(interac_ppi_unique.values.tolist(), directed=False)
igraph.summary(ppi_igraph)

# call the `clusters` method on the `ppi_igraph` object, and assign the
# resulting `VertexClustering` object to have object name `ppi_components`
ppi_components = ppi_igraph.clusters()

# call the `sizes` method on the `ppi_components` object, and assign the
# resulting list object to have the name `ppi_component_sizes`.
ppi_component_sizes = ppi_components.sizes()
print("component sizes, in order of component ID: ")
print(ppi_component_sizes)

# find its maximum value using the `max` built-in function
print("maximum component size: ")
max(ppi_component_sizes)

component_ids = ppi_components.membership
component_ids[0:5]

ctr = collections.Counter(component_ids)
print("from our Counter based counting: ")
print(dict(ctr))
print("from calling VertexClustering.sizes: ")
print({cid: csize for cid, csize in enumerate(ppi_component_sizes)})

list(np.argsort(-np.array(ppi_component_sizes)))

[id for id, _ in sorted(ctr.items(), reverse=True, key=operator.itemgetter(1))]

three_component_inds = [vertex_index for vertex_index, component_id in
                        enumerate(component_ids) if component_id == 3]
ppi_igraph.vs(three_component_inds)["name"]

# define a helper function to make an "edge key"
# (as a string) from the min(n,m) and max(n,m), separated by a hyphen
def make_edge_key(n: int, m: int) -> str:
  return str(min(n,m)) + '-' + str(max(n,m))

def is_eulerian_path(graph: List[Set], path: List[int]) -> bool:
    # create a dictionary `edge_counts` to hold the count of each edge "key"
    edge_counts = dict()

    # iterate over enumerate(graph), assigning n to index and neighbor_set to the value
    for n, neighbor_set in enumerate(graph):
        # if n is in `neighbor_set`, this is not a simple graph; raise ValueError
        if n in neighbor_set:
          raise ValueError("this graph contains a loop, and therefore is not simple")
        # iterate over neighbors m of n:
        for m in neighbor_set:
            # make an "edge key" (as a string) using the `make_edge_key` function
            edge_key = make_edge_key(n, m)
            # for that edge key, set the value in `edge_counts` to 0 (zero)
            edge_counts[edge_key] = 0
    # iterate over all indexes of entries in the path, except the last:
    for i in range(0, len(path) - 1):
        # get the starting vertex for this edge
        n = path[i]
        # get the ending vertex for this edge
        m = path[i + 1]
        # construct the edge key, as above
        edge_key = make_edge_key(n, m)
        # increment the counter for this edge key
        edge_counts[edge_key] += 1
    # every value in the counter dictionary should be exactly 1
    return set(list(edge_counts.values())) == {1}

print(is_eulerian_path([{1,2},{0},{0}], [2, 0]))
print(is_eulerian_path([{1,2},{0},{0}], [2, 0, 1]))