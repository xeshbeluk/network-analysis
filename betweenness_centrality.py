!apt-get install libcairo2-dev libjpeg-dev libgif-dev
!pip install pycairo
import cairo
!pip install python-igraph
import collections, random, igraph, pprint

from typing import List, Tuple

random.seed(1337)

g = igraph.Graph.Barabasi(10)
igraph.drawing.plot(g, bbox=[0, 0, 300, 300], vertex_label=list(range(0, 10)))

g_adjlist = g.get_adjlist()
pprint.pprint(dict(enumerate(g_adjlist)))

def single_vertex_bfs_with_paths_and_weights(s: int, g: List[List[int]]) -> \
Tuple[List[int], List[List[int]], List[int], List[int]]:
    # Let's review the function arguments:
    # `s` is the index of the vertex that we are at (type `int`)
    # `g` is the graph in adjacency list format (type `List[List[int]]`)

    # get the integer number of vertices in the graph (length of the list `g`)
    # and assign to `N`
    N = len(g)

    # define a `list` of length `N` containing the integer value -1,
    # [-1, -1, ... -1]; use the * operator; assign it to variable `dists`
    dists = [-1]*N

    # we know the distance (0) from vertex "s" to "s" (as the ending vertex);
    # so, set that entry in the `dists` list to value zero
    dists[s] = 0

    # define a `list` of length `N` containing the integer value 0,
    # [0, 0, ..., 0]; use the * operator; assign it to variable `weights`
    weights = N*[0]

    # we'll say that there is a single path from s to itself, so set the `s`
    # entry of `weights` to integer value 1
    weights[s] = 1

    # define a `list` of N empty lists, to contain the next-hop backpath edges to
    # get back to `s` from each vertex; assign it to the variable `paths`
    paths = [ [] for i in range(0,N)]


    # create an empty `collections.deque` object and assign it to the variable
    # `work`
    work = collections.deque()

    # add the vertex `s` to the queue `work`, using the `collections.deque.append`
    # method

    work.append(s)

    # create an empty `list` to contain the vertex IDs in decreasing order of
    # the vertices' geodesic distances from `s`; assign to the variable `orders`
    orders = []

    # start a `while` loop that will run until `work` has zero length (`len()`)
    while len(work) > 0:
        # grab the next vertex in the `work` queue using the `popleft` method,
        # and assign it to the variable `u`
        u = work.popleft()

        # since `u` is an integer vertex ID, use it to index into the adjacency
        # list `g` to get the list of vertex IDs of neighbors of `u`; assign that
        # list to variable `u_neighbors`
        u_neighbors = g[u]

        # iterate over `u_neighbors`, each time assigning the element to `v`
        for v in u_neighbors:
            # branch on whether list `dists` at position `v` is negative
            # (meaning that we have not visited vertex `v` before)
            if dists[v] < 0:  # we haven't visited this vertex before
                # compute dists[u]+1 and set dists[v] to this value
                dists[v] = dists[u] + 1

                # set weights[v] to the value weights[u]
                weights[v] = weights[u]

                # add `v` to the queue `work` using the `append` method
                work.append(v)

                # add item `v` to the *end* of the `orders` list (build up orders
                # in the reverse of the orientation we want, then *later* we will
                # reverse it before returning it in the final tuple return value)
                orders.append(v)

                # construct a list with `u` appended to `paths[v]`, and
                # store that list in `paths[v]`
                paths[v] = paths[v] + [u]
            else:             # we have visited this vertex before
                # check if there is already a path of length `dists[u] + 1`
                if dists[v] == dists[u] + 1:
                    # construct a list with `u` appended to `paths[v]`, and
                    # store that list in `paths[v]`
                    paths[v] = paths[v] + [u]

                    # increment `weights[v]` by `weights[u]`
                    weights[v] += weights[u]

    # return a tuple of `dists`, `paths`, `weights`, `orders[::-1]`;
    # the [::-1] slice syntax reverses a list
    return (dists, paths, weights, orders[::-1])

single_vertex_bfs_with_paths_and_weights(0, g_adjlist)

tg = igraph.Graph.TupleList([[0, 1], [1, 2],[2, 3], [2, 4],[3, 5],[4, 6],[4, 7],[6, 8], [5, 9],[7, 8],[8, 10], [9, 10]])
igraph.drawing.plot(tg, bbox=[0, 0, 150, 150], vertex_label=range(0, 11))

tg_adjlist = tg.get_adjlist()
print(tg_adjlist)

single_vertex_bfs_with_paths_and_weights(0, tg_adjlist)

def all_vertices_betweenness_centrality(g: List[List[int]]) -> List[float]:
    # compute the number of vertices as the `len` of list `g` and assign to `N`
    N = len(g)

    # define a list [0., 0., ..., 0.] of length `N`, and assign to `final_scores`
    final_scores = [0.]*N

    # for `s` in 0, ..., N-1: (this is pseudocode)
    for s in range(0,N):
        # call `single_vertex_bfs_with_paths_and_weights` with `s`, `g`;
        # assign the return tuple to (dists, paths, weights, orders)
        (dists, paths, weights, orders) = single_vertex_bfs_with_paths_and_weights(s,g)

        # define a list [1., 1., ..., 1.] of length `N`, and assign to `scores`
        scores = [1.]*N

        # iterate over `orders`, assigning each element to `v`:
        for v in orders:

            # get the list element at `paths[v]` and assign to `neighbors`
            neighbors = paths[v]

            # iterate over each item in `neighbors`, and assign to `neighbor`:
            for neighbor in neighbors:
                # compute the ratio `weights[neighbor]` to `weights[v]` and
                # assign it to the variable `ratio`
                ratio = weights[neighbor] / weights[v]

                # increment `scores[neighbor]` by `scores[v]*ratio`
                scores[neighbor] += scores[v]*ratio

        # iterate over `range(0, N)` and assign each to variable `u`:
        for u in range(0, N):
            # increment `final_scores[u]` by `scores[u]`
            final_scores[u] += scores[u]

    # iterate over `range(0, N)` and assign each element to variable `s`:
    for s in range(0, N):
        # compute (final_scores[s] - 2*N + 1)/2.0 and assign to `final_scores[s]`
        final_scores[s] = (final_scores[s] - 2*N + 1)/2.0

    # return `final_scores`
    return final_scores

all_vertices_betweenness_centrality(g_adjlist )

g.betweenness(directed=False)