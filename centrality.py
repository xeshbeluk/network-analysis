!pip install python-igraph
import igraph, timeit, matplotlib.pyplot
import numpy as np
import pandas as pd

!curl https://csx46.s3-us-west-2.amazonaws.com/PathwayCommons9.All.hgnc.sif.gz --output PathwayCommons9.All.hgnc.sif.gz
!gunzip -f PathwayCommons9.All.hgnc.sif.gz
sif_data = pd.read_csv("PathwayCommons9.All.hgnc.sif",
                       sep="\t", names=["species1","interaction_type","species2"])

interac_grn = sif_data[sif_data.interaction_type == "controls-expression-of"]
interac_grn_unique = interac_grn[["species1","species2"]].drop_duplicates()

grn_igraph = igraph.Graph.TupleList(interac_grn_unique.values.tolist(), directed=False)
grn_igraph.summary()

# get the number of vertices _N_
N = len(grn_igraph.vs)

# allocate a vector to contain the vertex closeness centralities; initialize to zeroes
# (so if a vertex is a singleton we don't have to update its closeness centrality)
closeness_centralities = np.zeros(N)

# initialize a counter to zero
ctr = 0

# start the timer
start_time = timeit.default_timer()

# for each `my_vertex` in `grn_igraph.vs`
for my_vertex in grn_igraph.vs:

    # compute the geodesic distance to every other vertex, from my_vertex, using
    # the `igraph.Graph.distances` instance method and passing `my_vertex`
    # as the `source` named parameter; put the result in a numpy.array (dtype float)
    my_dists = np.array(grn_igraph.distances(source=my_vertex))

    # filter the numpy array to include only entries that are nonzero and finite,
    # using `> 0 & numpy.isfinite(...)`
    my_dists = my_dists[np.isfinite(my_dists) & (my_dists > 0.0)]

    # if there are any distance values that survived the filtering, take their
    # element-wise reciprocals, then compute the sum, then divide by N-1
    # (following Eq. 7.30 in Newman)
    if len(my_dists) > 0:
        closeness_centralities[ctr] = np.sum(1.0/my_dists)/(N - 1.0)

    # increment ctr
    ctr += 1

# compute and print the elapsed time
ci_elapsed =  timeit.default_timer() - start_time
print(ci_elapsed)

matplotlib.pyplot.hist(closeness_centralities)
matplotlib.pyplot.xlabel("C")
matplotlib.pyplot.ylabel("Freq")
matplotlib.pyplot.show()

ax = matplotlib.pyplot.gca()
ax.scatter(grn_igraph.degree(), closeness_centralities)
ax.set_xscale("log")
matplotlib.pyplot.xlabel("degree")
matplotlib.pyplot.ylabel("closeness")
matplotlib.pyplot.show()

max_cc = np.max(closeness_centralities)
print(f"Maximum closeness centrality value is: {max_cc:.3f}")
print("Maximum closeness centrality protein has gene name: " +
      grn_igraph.vs[np.argmax(closeness_centralities)]["name"])

grn_igraph.vs[np.argsort(closeness_centralities)[::-1][0:9].tolist()]["name"]

cc_df = pd.DataFrame(list(zip(grn_igraph.vs["name"],
                              closeness_centralities.tolist())),
                     columns=["protein","CC"])
cc_df = cc_df.set_index("protein")
cc_df.sort_values("CC", ascending=False).head(n=12)

import random

!curl https://csx46.s3-us-west-2.amazonaws.com/neph_gene_network.txt > neph_gene_network.txt

!head neph_gene_network.txt

edge_list_neph = pd.read_csv("neph_gene_network.txt", sep="\t", names=["regulator","target"])

neph_graph = igraph.Graph.TupleList(edge_list_neph[["target","regulator"]].values.tolist(), directed=True)
neph_graph.summary()

pageranks = np.array(neph_graph.pagerank())

max_pr = np.max(pageranks)
print(f"Maximum Pagerank centrality value is: {max_pr:.2f}")
print(" the highest Pagerank centrality is: " + neph_graph.vs[np.argmax(pageranks)]["name"])

ax = matplotlib.pyplot.gca()
ax.set_xscale("log")
ax.set_yscale("log")
degrees = np.array(neph_graph.indegree())
inds_keep = np.where(degrees > 0)
ax.scatter(degrees[inds_keep],
           pageranks[inds_keep])
matplotlib.pyplot.xlabel("degree")
matplotlib.pyplot.ylabel("pagerank")
matplotlib.pyplot.show()

def pagerank(g):
    # N is the number of vertices
    N = len(g.vs)

    # alpha is the damping parameter, 0.85
    alpha = 0.85

    # beta = (1-alpha)/N
    beta = (1.0 - alpha) / N

    # compute the out-degree of each vertex, using igraph.Graph.degree with mode=igraph.OUT
    degree_values = g.degree(mode = 'out')

    # Construct a floating-point adjacency matrix M=A*D^(-1) in the Newman orientation,
    M = np.matrix(g.get_adjacency().data).transpose().astype(float)

    # and then we will do the "multiply by D^(-1)" step column by column:
    for j in range(0, N):

        # get the out degree of the vertex as a float
        degree_value = float(degree_values[j])

        # if degree is nonzero, normalize the column of M
        # otherwise, set the column to zero
        if degree_value > 0:
            M[:,j] /= degree_value # FILL IN HERE
        else:
            M[:,j] = 0# FILL IN HERE

    # Now implement Newman Eq. 7.19, by computing the difference of I (the
    # identity matrix) and alpha*M, and then inverting it, and then finally
    # multiplying by the beta vector and taking the vector transpose.
    pr = np.linalg.inv(np.diag([1.]*N) - alpha * M) * np.matrix([beta]*N).transpose()

    # Normalize the pagegrank centrality vector by dividing by the sum, so the
    # vector has unit L1 norm.
    pr = pr/np.sum(pr)

    # return the transpose vector, using the numpy.matrix.transpose method
    return(pr.transpose())

g = igraph.Graph.Barabasi(n=5, m=2)
print(pagerank(g).tolist()[0])
print(g.pagerank())