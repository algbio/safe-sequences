from itertools import count
from graphviz  import Digraph
import networkx as nx
import graph
import random
import logging
import numpy as np

logger = logging.getLogger(__name__)
inf    = 1 << 32

def read_graph(graph_raw):
    #Input format is: ['#Graph id\n', 'n\n', 'u_1 v_1 w_1\n', ..., 'u_k v_k w_k\n']
    id = graph_raw[0][7:]
    id = id[:len(id)-1]
    n  = int(graph_raw[1])
    G  = graph.st_DAG(n+2, 0, n+1, id) #+2 because of Source and Sink

    if n == 0:
        logging.warning("Graph %s has 0 vertices.", G.id)
        return G

    for edge in graph_raw[2:]:
        u,v,w = list(map(lambda x : int(x), edge.split(" ")))
        G.add_edge(u+1,v+1,w)

    #add edges (S,s) and (t,T) for all sources s and sinks t in the original DAG
    sources = G.get_original_sources()
    sinks   = G.get_original_sinks()
    for s in sources:
        G.add_edge(G.source,      s, G.outflow(s))
    for t in sinks:
        G.add_edge(       t, G.sink, G.inflow(t))

    G.w = max_edge_antichain(G)

    return G


def read_graphs(filename):
    f      = open(filename, "r")
    lines  = f.readlines()
    f.close()
    graphs = []

    #Assume: every file contains at least one graph
    i,j = 0,1 
    while True:
        if lines[j].startswith("#"):
            graphs.append(read_graph(lines[i:j]))
            i = j
        j += 1
        if (j==len(lines)):
            graphs.append(read_graph(lines[i:j]))
            break

    return graphs
    

def ER_st_DAG(n:int, p:float) -> graph.st_DAG :
    G = graph.st_DAG(n+2, 0, n+1, "ER_"+str(p))
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            if random.random() <= p:
                G.add_edge(i,j,1)
    sources = G.get_original_sources()
    sinks   = G.get_original_sinks()
    for s in sources:
        G.add_edge(0,     s, G.outflow(s))
    for t in sinks:
        G.add_edge(t, G.n-1, G.inflow(t))
    G.w = max_edge_antichain(G)
    return G
    

def is_0_flow_everywhere(G : graph.st_DAG) -> bool:
    is_0_everywhere = True
    for edge in G.flow:
        if G.flow[edge]!=0:
            is_0_everywhere=False
            break
    return is_0_everywhere


def visualize(G : graph.st_DAG, weighted_paths=[], safe_sequences=[], tag = ''):
    dot = Digraph(format='pdf')
    dot.graph_attr['rankdir'] = 'LR'        # Display the graph in landscape mode
    dot.node_attr['shape']    = 'rectangle' # Rectangle nodes

    E = G.edge_list
    colors = ['red','blue','green','purple','brown','cyan','yellow','pink','grey']

    for (u,v) in E:
        dot.edge(str(u),str(v),label=str(G.flow[(u,v)]))
    
    i=0
    l = []
    for path,weight,slack in weighted_paths:
        path = list(zip(path, path[1:]))
        pathColor = colors[i % len(colors)]
        for (u,v) in path:
            dot.edge(str(u), str(v), fontcolor=pathColor, color=pathColor, penwidth='2.0', label=str(weight))
        i=i+1
        l.append((pathColor,weight,slack))

    with dot.subgraph(name='legend') as legend:
        legend.attr(label='Legend', style='dashed')

        # Invisible nodes with color
        for (c,w,s) in l:
            legend.node(c, str(w) + " " + str(s), style='filled', fillcolor='lightblue', shape='box')

    for sequence in safe_sequences:
        pathColor = colors[i % len(colors)]
        for (u,v) in sequence:
            dot.edge(str(u), str(v), style='dashed', fontcolor=pathColor, color=pathColor, penwidth='1.0')
        i=i+1

    dot.render(filename=G.id+tag,directory='.', view=True)


def min_cost_flow(G, s, t):
    
    flowNetwork = nx.DiGraph()
    
    flowNetwork.add_node(s, demand = -inf)
    flowNetwork.add_node(t, demand = inf)
            
    for v in G.nodes():
        if v != s and v != t:
            flowNetwork.add_node(v, demand = 0)
    
    flowNetwork.add_edge(s, t, weight = 0)

    counter = count(1) # Start an iterator given increasing integers starting from 1
    edgeMap = dict()
    
    for (x,y) in G.edges():
        z1 = str(next(counter))
        z2 = str(next(counter))
        edgeMap[(x,y)] = z1
        l = G[x][y]['l']
        u = G[x][y]['u']
        c = G[x][y]['c']
        flowNetwork.add_node(z1, demand = l)
        flowNetwork.add_node(z2, demand = -l)
        flowNetwork.add_edge(x, z1, weight = c, capacity = u)
        flowNetwork.add_edge(z1, z2, weight = 0, capacity = u)
        flowNetwork.add_edge(z2, y, weight = 0, capacity = u)

    flowCost, flowDictNet = nx.network_simplex(flowNetwork)
    
    flowDict = dict()
    for x in G.nodes():
        flowDict[x] = dict()

    for (x,y) in G.edges():
        flowDict[x][y] = flowDictNet[x][edgeMap[(x,y)]]

    return flowCost, flowDict


def max_edge_antichain(G_original : graph.st_DAG, get_antichain = False, weight_function = {}) -> list :

    G_nx       = nx.DiGraph()
    new_source = 0
    new_sink   = G_original.n+1
    G          = graph.st_DAG(G_original.n+2, new_source, new_sink, G_original.id+str("_tmp"))
    demand     = dict()

    G_nx.add_node(new_source)
    G_nx.add_node(new_sink)

    for (u,v) in G_original.edge_list:
        G.add_edge(u+1,v+1,1)
        G_nx.add_node(u+1)
        G_nx.add_node(v+1)
        demand[(u+1,v+1)] = weight_function[(u,v)] if weight_function else 1
        G_nx.add_edge(u+1, v+1, l = demand[(u+1,v+1)], u=inf, c=0)

    for v in G.get_nodes_but_st():
        G_nx.add_edge(G.source,        v, l=0, u=inf, c=1)
        G_nx.add_edge(       v,   G.sink, l=0, u=inf, c=0)
        G.add_edge(G.source,v,1)
        G.add_edge(v,G.sink,1)
        demand[(G.source,v)] = 0
        demand[(v,G.sink)]   = 0

    flowCost, flow = min_cost_flow(G_nx, G.source, G.sink)

    def DFS_find_reachable_from_source(u,visited):
        if visited[u]!=0:
            return
        assert(u!=G.sink)
        visited[u] = 1
        for v in G.out_neighbors(u):
            if flow[u][v] > demand[(u,v)]:
                DFS_find_reachable_from_source(v, visited)
        for v in G.in_neighbors(u):
            DFS_find_reachable_from_source(v,visited)

    def DFS_find_saturating(u,visited):
        if visited[u] != 1:
            return
        visited[u] = 2
        for v in G.out_neighbors(u):
            if flow[u][v] > demand[(u,v)]:
                DFS_find_saturating(v, visited)
            elif flow[u][v] == demand[(u,v)] and demand[(u,v)]>=1 and visited[v]==0:
                antichain.append((u-1,v-1))
        for v in G.in_neighbors(u):
            DFS_find_saturating(v,visited)

    if get_antichain:
        antichain = []
        visited   = [0] * G.n
        DFS_find_reachable_from_source(G.source, visited)
        DFS_find_saturating(G.source, visited)
        if weight_function:
            assert(flowCost == sum(map(lambda edge : weight_function[edge], antichain)))
        else:
            assert(flowCost == len(antichain))
        return flowCost,antichain
    
    return flowCost


class GRB_TimeOut(Exception):
    def __init__(self, message:str):
        super(GRB_TimeOut, self).__init__('TimeOut: ' + message)


class GRB_Infeasible(Exception):
    def __init__(self, message:str):
        super(GRB_Infeasible, self).__init__('Infeasibility: ' + message)



def metrics(G : graph.st_DAG):
    flow_values = list(G.flow.values())

    average = np.mean(flow_values)
    std_dev = np.std(flow_values)
    minimum = np.min(flow_values)
    maximum = np.max(flow_values)
    sum_values = np.sum(flow_values)
    median = np.median(flow_values)
    range_values = maximum - minimum
    variance = np.var(flow_values)
    percentile_25 = np.percentile(flow_values, 25)
    percentile_75 = np.percentile(flow_values, 75)

    # Print results
    print(G.id)
    print(f"Average: {average}")
    print(f"Standard Deviation: {std_dev}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Sum: {sum_values}")
    print(f"Median: {median}")
    print(f"Range (Max - Min): {range_values}")
    print(f"Variance: {variance}")
    print(f"25th Percentile: {percentile_25}")
    print(f"75th Percentile: {percentile_75}")
    print()
    return