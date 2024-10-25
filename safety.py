import graph
import logging
from queue import Queue

logger    = logging.getLogger(__name__)

def find_all_bridges(adj_list, s, t) -> list:

    #find arbitrary s-t path
    s_aux = s
    p = [s_aux]
    while s_aux!=t:
        x = adj_list[s_aux].pop() #remove edge in O(1)
        p.append(x)               #keep track of the path
        s_aux = x
    #add reversed path to G
    for i in range(len(p)-1):
        adj_list[p[i+1]].append(p[i])

    n            = len(adj_list)
    i            = 1
    bridges      = []
    component    = [0] * n
    q            = Queue(maxsize = n+1)
    component[s] = 1
    first_node   = 0
    q.put(s)

    while component[t]==0: #do while :(

        if i!=1:
            #find first node u of P with component[u]=0. all in all we pay |P| time for this
            while component[p[first_node]] != 0:
                first_node += 1

            y = p[first_node-1]
            z = p[first_node  ]

            bridges.append((y,z))
            q.put(z)
            component[z] = i

        while not q.empty():
            u = q.get()
            for v in adj_list[u]:
                if component[v]==0:
                    q.put(v)
                    component[v]=i
        i = i+1

    #recover original adjacency relation
    for i in range(len(p)-1):
        u,v = p[i],p[i+1]
        adj_list[v].pop()      #remove reversed edges
        adj_list[u].append(v)  #reinsert removed edges

    return bridges


def safe_sequences(G : graph.st_DAG) -> list :
    
    sequences = []

    for (u,v) in G.edge_list:
        left_extension  = find_all_bridges(G.get_adj_list_R(), u, G.source)
        right_extension = find_all_bridges(G.get_adj_list()  , v, G.sink  )

        for i in range(len(left_extension)): #reverse edges of left extension, recall G^R
            x,y = left_extension[i]
            left_extension[i] = (y,x)

        seq = left_extension[::-1] + [ (u,v,) ] + right_extension

        sequences.append(seq)

    return sequences


def safe_paths(G : graph.st_DAG) -> list :
    
    paths = []
    for e in G.edge_list:
        path = []
        u,v = e
        
        while G.unique_in_neighbor(u):
            x = G.in_neighbors(u)[0]
            path.append( (x,u) )
            u = x

        path = path[::-1]
        path.append(e)
        
        while G.unique_out_neighbor(v):
            x = G.out_neighbors(v)[0]
            path.append( (v,x) )
            v = x

        paths.append(path)

    return paths
