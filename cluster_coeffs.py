# -*- coding: utf-8 -*-
!apt-get install libcairo2-dev libjpeg-dev libgif-dev
!pip install pycairo
!pip install cython
!pip install python-igraph pympler
!pip install --only-binary=:all: bintrees
import cairo, igraph, pandas, numpy, timeit, pympler.asizeof, bintrees, matplotlib

!curl https://csx46.s3-us-west-2.amazonaws.com/PathwayCommons9.All.hgnc.sif.gz --output PathwayCommons9.All.hgnc.sif.gz
!gunzip -f PathwayCommons9.All.hgnc.sif.gz

sif_data = pandas.read_csv("PathwayCommons9.All.hgnc.sif",
                           sep="\t", names=["species1", "interaction_type", "species2"])

interaction_types_ppi = {"interacts-with", "in-complex-with"}
interac_ppi = sif_data[sif_data.interaction_type.isin(interaction_types_ppi)]

boolean_vec = interac_ppi['species1'] > interac_ppi['species2']
interac_ppi.loc[boolean_vec, ['species1', 'species2']] = interac_ppi.loc[boolean_vec, ['species2', 'species1']].values
interac_ppi.head()

interac_ppi_unique = interac_ppi[["species1","species2"]].drop_duplicates()

ppi_igraph = igraph.Graph.TupleList(interac_ppi_unique.values.tolist(), directed=False)
igraph.summary(ppi_igraph)

ppi_adj_list = ppi_igraph.get_adjlist()

def get_bst_forest(theadjlist):
    g_adj_list = theadjlist
    n = len(g_adj_list)
    theforest = []
    for i in range(0, n):
        itree = bintrees.AVLTree()
        for j in g_adj_list[i]:
          itree.insert(j,1)
        theforest.append(itree)
            # insert at j,1
        # append itree to theforest
    return theforest

def find_bst_forest(bst_forest, i, j):
    return j in bst_forest[i]
    # is j in bst_forest[i]?

ppi_adj_forest = get_bst_forest(ppi_adj_list)

N = len(ppi_adj_list)
civals = numpy.zeros(100)
civals[:] = numpy.NaN
start_time = timeit.default_timer()
for n in range(0, 100): # loop over all vertices
    neighbors = ppi_adj_list[n]
    nneighbors = len(neighbors)
    if nneighbors > 1:
        # initialize nctr to 0
        nctr = 0
        for i in range(0, nneighbors):# 0 to nneighbors
            for j in range(i + 1, nneighbors):# 1+1 to nneighbors
                if neighbors[j] in ppi_adj_forest[neighbors[i]]:
                    nctr += 1 # check for j neighbor of i
                    # increment nctr
        civals[n] = nctr/(nneighbors*(nneighbors-1)/2)

ci_elapsed = timeit.default_timer() - start_time
# print elapsed time, to two decimal points
print("%0.2f s" % ci_elapsed)

start_time = timeit.default_timer()
civals_igraph = ppi_igraph.transitivity_local_undirected(vertices=list(range(0,100)))
# remember to restrict to the first 100 vertices!!
ci_elapsed = timeit.default_timer() - start_time
# print elapsed time, to two decimal points
print("%0.2f s" % ci_elapsed)

# plot civals vs. civals_igraph, using matplotlib
matplotlib.pyplot.plot(civals, civals_igraph)
matplotlib.pyplot.xlabel("Ci (my code)")
matplotlib.pyplot.ylabel("Ci (igraph)")
matplotlib.pyplot.show()

civals_igraph = numpy.array(ppi_igraph.transitivity_local_undirected())
deg_igraph = ppi_igraph.degree()
deg_npa = numpy.array(deg_igraph)
deg_binids = numpy.rint(deg_npa/50)
binkvals = 50*numpy.array(range(0,25))
civals_avg = [numpy.mean(civals_igraph[deg_binids == j]) for j in range(0, 25)]# for each j in 0 to 24, compute mean of vertex Ci vals for all vertices i for which the vertex's degree is in degree bin j

# on log-log scale, using matplotlib.pyplot, plot civals_avg vs. binkvals
matplotlib.pyplot.loglog(
     binkvals,
     civals_avg)
matplotlib.pyplot.ylabel("<Ci>")
matplotlib.pyplot.xlabel("k")
matplotlib.pyplot.show()

civals = numpy.zeros(len(ppi_adj_list))
civals[:] = numpy.NaN

# build an adjacency "list of hashtables" (list of length N = number of vertices, containing a python `set` of neighbors)
ppi_adj_hash = []
for i in range(0, len(ppi_adj_list)):
    newhash = set()
    for j in ppi_adj_list[i]:
        newhash.add(j)
    ppi_adj_hash.append(newhash)
        # add j to newhash
    # append newhash to ppi_adj_hash

start_time = timeit.default_timer()
for n in range(0, len(ppi_adj_list)):
    neighbors = ppi_adj_hash[n]
    nneighbors = len(neighbors)
    if nneighbors > 1:
        nctr = 0
        for i in neighbors:
            for j in neighbors:
                if (j > i) and (j in ppi_adj_hash[i]):
                    nctr += 1
        civals[n] = nctr/(nneighbors*(nneighbors-1)/2)

ci_elapsed = timeit.default_timer() - start_time
print("%0.2f" % ci_elapsed)

pympler.asizeof.asizeof(ppi_adj_hash)/1000000

pympler.asizeof.asizeof(ppi_adj_forest)/1000000