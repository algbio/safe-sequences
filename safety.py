import graph
import logging
import dominators
from queue import Queue

logger    = logging.getLogger(__name__)

def find_idom(adj_list, s, t) -> list:

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
    component    = [0] * n
    q            = Queue(maxsize = n+1)
    component[s] = 1
    first_node   = 0
    first_bridge = None
    q.put(s)

    while component[t]==0: #do while :(

        if i!=1:
            #find first node u of P with component[u]=0. all in all we pay |P| time for this
            while component[p[first_node]] != 0:
                first_node += 1

            first_bridge = ( p[first_node-1] ,p[first_node] )
            break

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

    return first_bridge


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


def all_diff(X : list) -> bool:
    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            if X[i] == X[j]:
                return False
    return True


def find_unitig_of_arc(G : graph.st_DAG, e : tuple):
    u,v = e
    #assert(G.is_edge(e))
    unitig = [(u,v)]
    while G.has_unique_out_neighbor(v) and G.has_unique_in_neighbor(v):
        x = G.out_neighbors(v)[0]
        unitig.append((v,x))
        v = x
    while G.has_unique_in_neighbor(u) and G.has_unique_out_neighbor(u):
        x = G.in_neighbors(u)[0]
        unitig = [(x,u)] + unitig #zz...
        u = x
    return u,v,unitig


def is_core(G : graph.st_DAG, u: int, v: int) -> bool:
    return (G.out_degree(v) < 1 or G.in_degree(v) != 1) and (G.in_degree(u) < 1 or G.out_degree(u) != 1)


def maximal_safe_sequences(G : graph.st_DAG, X = []) -> list :

    sequences = []
    processed_arcs = set()

    for e in X:

        if e in processed_arcs:
            continue

        L,R,unitig = find_unitig_of_arc(G,e) #every arc-unitig has an identifying pair of leftmost and rightmost nonunivocal vertices (or at least one of them is the source or the sink)
        
        for arc in unitig:
            processed_arcs.add(arc)

        if is_core(G,L,R): #this unitig (or arc in the compressed graph) is a core
            left_extension  = find_all_bridges(G.get_adj_list_R(), L, G.source)
            right_extension = find_all_bridges(G.get_adj_list()  , R, G.sink  )

            for i in range(len(left_extension)): #reverse edges of left extension (recall the definition of G^R)
                x,y = left_extension[i]
                left_extension[i] = (y,x)

            seq = left_extension[::-1] + unitig + right_extension

            sequences.append(seq)

    return sequences



def maximal_safe_sequences_via_dominators(G : graph.st_DAG, X = set()) -> list :

    s_idoms = dict()
    t_idoms = dict()
    for (u,v) in G.edge_list:
        s_idom = find_idom(G.get_adj_list_R(), u, G.source)
        t_idom = find_idom(G.get_adj_list()  , v,   G.sink)
        s_idoms[(u,v)] = tuple(reversed(s_idom)) if s_idom != None else G.source
        t_idoms[(u,v)] = t_idom                  if t_idom != None else G.sink

    T_s = dominators.Arc_Dominator_Tree(G.n, G.source, s_idoms, G.edge_list, X, G.id+str("_s-domtree"))
    T_t = dominators.Arc_Dominator_Tree(G.n, G.sink  , t_idoms, G.edge_list, X, G.id+str("_t-domtree"))

    leaves_s_X = [ node for node,children in T_s.children_X.items() if len(children)==0 ] #those nodes in X that do not s-dominate other nodes with respect to X

    cores = []
    for leaf in leaves_s_X: # O(m): the paths are unitary, so no node can belong to two distinct s- (or t-) unitary paths
        s_unitary_path = T_s.find_unitary_path_X(leaf,   "up")
        t_unitary_path = T_t.find_unitary_path_X(leaf, "down")
        
        if len(t_unitary_path) > len(s_unitary_path):
            continue

        i=0
        good_sequence = True
        while good_sequence and i<len(t_unitary_path):
            if s_unitary_path[i] != t_unitary_path[i]:
                good_sequence = False
            else:
                i=i+1
                
        if i == len(t_unitary_path) and T_t.is_leaf_X(t_unitary_path[len(t_unitary_path)-1]):
            assert(good_sequence==True)
            cores.append(leaf)

    maximal_safe_sequences = []
    for core in cores: # O(length of all maximal safe sequences), with no duplicates
        s_doms = T_s.get_dominators(core)
        t_doms = T_t.get_dominators(core)
        maximal_safe_sequences.append( s_doms[::-1] + t_doms[1:] )
    
    return maximal_safe_sequences
