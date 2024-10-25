class st_DAG:

    def __init__(self, n:int, source:int, target:int, id:str):
        self.id        = id
        self.n         = n
        self.m         = 0
        self.w         = 0
        self.graph     = [[] for _ in range(self.n)]
        self.graph_R   = [[] for _ in range(self.n)]
        self.edge_list = []
        self.flow      = dict()

        self.source = source
        self.sink   = target

    def get_adj_list(self):
        return self.graph

    def get_adj_list_R(self):
        return self.graph_R

    #maybe make weight argument optional
    def add_edge(self,u,v,w):
        self.graph[u].append(v)
        self.graph_R[v].append(u)
        self.edge_list.append((u,v))
        self.flow[(u,v)] = w
        self.m += 1
    
    def is_source(self,u):
        return len(self.graph_R[u]) == 0 and u!=self.source and u!=self.sink

    def is_sink(self,u):
        return len(self.graph[u]) == 0   and u!=self.source and u!=self.sink

    def get_sources(self):
        return list( filter ( self.is_source, [i for i in range(self.n)] ) )

    def get_sinks(self):
        return list( filter ( self.is_sink,   [i for i in range(self.n)] ) )

    def out_neighbors(self,u):
        return self.graph[u]

    def in_neighbors(self,u):
        return self.graph_R[u]
    
    def out_degree(self,u) -> int:
        return len(self.out_neighbors(u))

    def in_degree(self,u) -> int:
        return len(self.in_neighbors(u))

    def unique_out_neighbor(self,u) -> bool:
        return self.out_degree(u) == 1

    def unique_in_neighbor(self,u) -> bool:
        return self.in_degree(u) == 1
    
    def outflow(self,u) -> int:
        out_neighbors = self.out_neighbors(u)
        return sum(map(lambda v : self.flow[(u,v)], out_neighbors))
    
    def inflow(self,v) -> int:
        in_neighbors = self.in_neighbors(v)
        return sum(map(lambda u : self.flow[(u,v)], in_neighbors))

    def excess(self,u) -> int:
        return self.inflow(u) - self.outflow(u)

    def get_nodes(self) -> list:
        return list(range(self.n))

    def get_nodes_but_st(self) -> list:
        return list(range(1,self.n-1))

    def print(self):
        print(">>>Graph {} n={} m={}".format(self.id, self.n, self.m))
        for u in range(self.n):
            l = "{} -> ".format(u)
            for v in self.out_neighbors(u):
                l += "({},{}); ".format(v,self.flow[(u,v)])
            print(l)

    def __str__(self):
        G = ">>>Graph {} n={} m={}\n".format(self.id, self.n, self.m)
        for u in range(self.n):
            l = "{} -> ".format(u)
            for v in self.out_neighbors(u):
                l += "({},{}); ".format(v,self.flow[(u,v)])
            G += l + "\n"
        return G

class EdgeList:
    def __init__(self, n:int):
        self.n = n
        self.m = 0
        self.edge_list = []

    def add_edge(self,u,v):
        self.edge_list.append((u,v))
        self.m += 1    

    def append_edges(self, edges):
        self.edge_list = edges + self.edge_list

    def prepend_edges(self, edges):
        self.edge_list = self.edge_list + edges
