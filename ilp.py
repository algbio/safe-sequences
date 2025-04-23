from gurobipy import GRB
import gurobipy as gp
import logging
import graph
import utils

logger    = logging.getLogger(__name__)
TOLERANCE = 0.1  #tolerance allowed for Gurobi numerical values

def tail(grb_edge): return int(grb_edge.split("[")[1].split(",")[0])
def head(grb_edge): return int(grb_edge.split("[")[1].split(",")[1])


class Encode_LeastSquares:

    def __init__(self,n,E,source,sink,F,edge_width,R,P2F,epsilon,timeout,threads):
        self.n          = n
        self.m          = len(E)
        self.source     = source
        self.target     = sink
        self.E          = E
        self.F          = F
        self.R          = R
        self.vars2fix   = P2F
        self.epsilon    = epsilon
        self.k          = edge_width #the largest edge antichain is a lower bound for the size of any flow decomposition
        self.w_max      = max(map(lambda edge : F[edge], self.E))
        self.edge_vars  = {}
        self.pi_vars    = {}
        self.weights    = {}
        self.spc_vars   = {}

        self.timeout    = timeout
        self.threads    = threads

        self.model      = self.create_solver()

        self.final_k    = None

    def create_solver(self):
        env = gp.Env(empty=True)
        env.setParam('OutputFlag'   ,              0)
        env.setParam('LogToConsole' ,              0)
        env.setParam('TimeLimit'    ,   self.timeout)
        env.setParam('Threads'      ,   self.threads)
        env.start()
        model = gp.Model("MFD_LeastSquares", env=env)
        if not model:
            logger.error("FATAL, could not create Gurobi model")
            exit(0)
        return model

    def clear(self):
        self.model       = self.create_solver()
        self.edge_vars   = {}
        self.pi_vars     = {}
        self.weights     = {}
        self.spc_vars    = {}

    def solve(self):
        self.model.optimize()

    def encode(self):

        # Create variables
        edge_indexes    = [ (u,v,i) for i in range(self.k) for (u, v) in self.E        ]
        path_indexes    = [ (    i) for i in range(self.k)                             ]
        subpath_indexes = [ (i,j  ) for i in range(self.k) for j in range(len(self.R)) ]

        self.edge_vars = self.model.addVars(   edge_indexes, vtype=GRB.BINARY ,  name='e'                     )
        self.pi_vars   = self.model.addVars(   edge_indexes, vtype=GRB.INTEGER,  name='p', lb=0, ub=self.w_max)
        self.weights   = self.model.addVars(   path_indexes, vtype=GRB.INTEGER,  name='w', lb=1, ub=self.w_max)
        self.spc_vars  = self.model.addVars(subpath_indexes, vtype=GRB.BINARY ,  name='r'                     )

        #The identifiers of the constraints come from https://arxiv.org/pdf/2201.10923 page 14-15

        for i in range(self.k):
            self.model.addConstr( self.edge_vars.sum(self.source,'*',i) == 1, "10a_i={}".format(i) )
            self.model.addConstr( self.edge_vars.sum('*',self.target,i) == 1, "10b_i={}".format(i) )

        for i in range(self.k):
            for v in range(1,self.n-1): #find all wedges u->v->w for v in V\{s,t}
                self.model.addConstr( self.edge_vars.sum('*',v,i) - self.edge_vars.sum(v,'*',i) == 0, "10c_v={}_i={}".format(v,i) )

        for (u,v) in self.E:
            for i in range(self.k):
                self.model.addConstr( self.pi_vars[u,v,i] <= self.edge_vars[u,v,i] * self.w_max                         , "10e_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.pi_vars[u,v,i] <= self.weights[i]                                            , "10f_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.pi_vars[u,v,i] >= self.weights[i] - (1 - self.edge_vars[u,v,i]) * self.w_max , "10g_u={}_v={}_i={}".format(u,v,i) )
                
        #Example of a subpath constraint: R=[ [(1,3),(3,5)], [(0,1)] ], means that we have 2 paths to cover, the first one is 1-3-5. the second path is just a single edge 0-1
        def EncodeSubpathConstraints():
            for i in range(self.k):
                for j in range(len(self.R)):
                    edgevars_on_subpath = list(map(lambda e: self.edge_vars[e[0],e[1],i], self.R[j]))
                    self.model.addConstr( sum(edgevars_on_subpath) >= len(self.R[j]) * self.spc_vars[i,j] )
            for j in range(len(self.R)):
                self.model.addConstr( self.spc_vars.sum('*',j) >= 1 )

        def Fix_Variables():
            for i in range(len(self.vars2fix)):
                for (u,v) in self.vars2fix[i]:
                    self.model.addConstr( self.edge_vars[u,v,i] == 1 )
        
        if self.R!=[]:
            EncodeSubpathConstraints()

        if self.vars2fix!=[]:
            Fix_Variables()

        self.model.setObjective( sum( (self.F[(u,v)] - self.pi_vars.sum(u,v,'*') )**2 for (u,v) in self.E), GRB.MINIMIZE )

    def print_solution(self,solution):
        opt, dif, paths = solution
        print("\n#####SOLUTION#####\n","> FD size   :",   opt,"\n> LeastSquares difference :", dif,"\n> Weight-Path decomposition:")
        for p in paths:
            print(*p)
    
    def build_solution(self):
        paths = []
        for i in range(self.k):
            path = []
            u    = self.source
            while u != self.target:
                edges = list(filter(lambda grb_var : grb_var.X>1-TOLERANCE , self.edge_vars.select(u,'*',i) ))
                assert(len(edges)==1)
                v = head(edges[0].VarName)
                path.append(v)
                u = v
            paths.append( (self.weights[i].X, path[:-1]) )

        return (self.k, self.model.ObjVal, paths)

    def solve_once(self):   
        logger.info(">>> Solving once")            
        self.encode()
        self.solve()
        logger.info("Gurobi solver status " + str(self.model.status))
        if self.model.status == GRB.TIME_LIMIT:
            raise utils.GRB_TimeOut("ilp.LeastSquares.solve_once while finding feasible solution")
        elif self.model.status == GRB.INFEASIBLE:
            raise utils.GRB_Infeasible("ilp.LeastSquares.solve_once while finding feasible solution")
        elif self.model.status == GRB.OPTIMAL:
            _,_,p = self.build_solution()
            weights,paths = list(map(lambda x : x[0], p)), list(map(lambda x : x[1], p))
            #return (paths,weights)
            return self.model.ObjVal
        else:
            logger.error("FATAL: solver finished with status " + str(self.model.status) + ", which is not: OPT, TIME_LIMIT, INFEASIBLE.")

    def optimize_linear(self):
        #Assumption: the initial value of self.k is sufficiently high so that the ILP solver starts with a feasible solution

        self.encode()
        self.solve()
        if self.model.status == GRB.TIME_LIMIT:
            self.final_k = self.k
            raise utils.GRB_TimeOut("ilp.LeastSquares.optimize_linear in the feasibility stage of the optimization loop")
        elif self.model.status == GRB.INFEASIBLE:
            self.final_k = self.k
            raise utils.GRB_Infeasible("ilp.LeastSquares.optimize_linear in the feasibility stage of the optimization loop")
        prev_obj = self.model.ObjVal

        if prev_obj == 0:
            logger.info(">>> No optimization step required, found perfect solution in init solving with " + str(self.k))
            self.final_k = self.k
            return self.final_k

        self.clear()
        self.k += 1

        logger.info(">>> Optimality, starting with " + str(self.k))
        while True: #optimality criteria: find the k for which the ratio of slacks between two consecutive iterations becomes sufficiently small

            self.encode()
            self.solve()
            logger.info("Gurobi solver status " + str(self.model.status))

            if self.model.status == GRB.TIME_LIMIT:
                self.final_k = self.k - 1
                raise utils.GRB_TimeOut("ilp.LeastSquares.optimize_linear while optimizing")
            elif self.model.status == GRB.INFEASIBLE:
                self.final_k = self.k - 1
                raise utils.GRB_Infeasible("ilp.LeastSquares.optimize_linear while optimizing")
            #assert(self.model.ObjVal <= prev_obj) not necessarily true, as the lower bound of the weight variables is 1

            if self.model.ObjVal >= prev_obj: #if we do not improve by allowing more paths we stop
                self.final_k = self.k-1
                break

            if self.model.ObjVal==0: #if we found a perfect solution we stop
                self.final_k = self.k
                break

            assert(prev_obj != 0)

            if 1-self.model.ObjVal/prev_obj < self.epsilon: # if the relative improvement is small we also stop
                self.final_k = self.k - 1
                break
            
            prev_obj = self.model.ObjVal
            self.clear()
            self.k += 1
        
        return self.final_k



class Encode_Robust:

    def __init__(self,n,E,source,sink,F,edge_width,R,P2F,epsilon,timeout,threads):
        self.n          = n
        self.m          = len(E)
        self.source     = source
        self.target     = sink
        self.E          = E
        self.F          = F
        self.R          = R
        self.vars2fix   = P2F
        self.epsilon    = epsilon
        self.k          = edge_width #the largest edge antichain is a lower bound for the size of any flow decomposition
        self.w_max      = max(map(lambda edge : F[edge], self.E))
        self.edge_vars  = {}
        self.phi_vars   = {}
        self.gam_vars   = {}
        self.weights    = {}
        self.slacks     = {}
        self.spc_vars   = {}

        self.timeout    = timeout
        self.threads    = threads

        self.model      = self.create_solver()

        self.final_k    = None

    def create_solver(self):
        env = gp.Env(empty=True)
        env.setParam('OutputFlag'   ,              0)
        env.setParam('LogToConsole' ,              0)
        env.setParam('TimeLimit'    ,   self.timeout)
        env.setParam('Threads'      ,   self.threads)
        env.start()
        model = gp.Model("MFD_Robust", env=env)
        if not model:
            logger.error("FATAL, could not create Gurobi model")
            exit(0)
        return model

    def clear(self):
        self.model       = self.create_solver()
        self.edge_vars   = {}
        self.phi_vars    = {}
        self.gam_vars    = {}
        self.weights     = {}
        self.slacks      = {}
        self.spc_vars    = {}

    def solve(self):
        self.model.optimize()

    def encode(self):

        # Create variables
        edge_indexes    = [ (u,v,i) for i in range(self.k) for (u, v) in self.E        ]
        path_indexes    = [ (    i) for i in range(self.k)                             ]
        subpath_indexes = [ (i,j  ) for i in range(self.k) for j in range(len(self.R)) ]

        self.edge_vars = self.model.addVars(   edge_indexes, vtype=GRB.BINARY ,  name='e'                     )
        self.spc_vars  = self.model.addVars(subpath_indexes, vtype=GRB.BINARY ,  name='r'                     )
        self.phi_vars  = self.model.addVars(   edge_indexes, vtype=GRB.INTEGER,  name='p', lb=0, ub=self.w_max)
        self.gam_vars  = self.model.addVars(   edge_indexes, vtype=GRB.INTEGER,  name='g', lb=0, ub=self.w_max)
        self.weights   = self.model.addVars(   path_indexes, vtype=GRB.INTEGER,  name='w', lb=1, ub=self.w_max)
        self.slacks    = self.model.addVars(   path_indexes, vtype=GRB.INTEGER,  name='s', lb=0, ub=self.w_max)

        #The identifiers of the constraints come from https://www.biorxiv.org/content/10.1101/2023.03.20.533019v1.full.pdf page 13

        for i in range(self.k):
            self.model.addConstr( self.edge_vars.sum(self.source,'*',i) == 1, "14a_i={}".format(i) )
            self.model.addConstr( self.edge_vars.sum('*',self.target,i) == 1, "14b_i={}".format(i) )

        for i in range(self.k):
            for v in range(1,self.n-1): #find all wedges u->v->w for v in V\{s,t}
                self.model.addConstr( self.edge_vars.sum('*',v,i) - self.edge_vars.sum(v,'*',i) == 0, "14c_v={}_i={}".format(v,i) )

        for (u,v) in self.E:
            f_uv    = self.F[(u,v)]
            phi_sum = self.phi_vars.sum(u,v,'*')
            gam_sum = self.gam_vars.sum(u,v,'*')
            self.model.addConstr( f_uv - phi_sum <=  gam_sum, "14d_u={}_v={}".format(u,v) )
            self.model.addConstr( f_uv - phi_sum >= -gam_sum, "14e_u={}_v={}".format(u,v) )
            for i in range(self.k):
                self.model.addConstr( self.phi_vars[u,v,i] <= self.w_max * self.edge_vars[u,v,i]                        , "14f_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.gam_vars[u,v,i] <= self.w_max * self.edge_vars[u,v,i]                        , "14i_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.phi_vars[u,v,i] <= self.weights[i]                                           , "14g_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.gam_vars[u,v,i] <= self.slacks [i]                                           , "14j_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.phi_vars[u,v,i] >= self.weights[i] - (1 - self.edge_vars[u,v,i]) * self.w_max, "14h_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.gam_vars[u,v,i] >= self.slacks [i] - (1 - self.edge_vars[u,v,i]) * self.w_max, "14k_u={}_v={}_i={}".format(u,v,i) )

        #Example of a subpath constraint: R=[ [(1,3),(3,5)], [(0,1)] ], means that we have 2 paths to cover, the first one is 1-3-5. the second path is just a single edge 0-1
        def EncodeSubpathConstraints():
            for i in range(self.k):
                for j in range(len(self.R)):
                    edgevars_on_subpath = list(map(lambda e: self.edge_vars[e[0],e[1],i], self.R[j]))
                    self.model.addConstr( sum(edgevars_on_subpath) >= len(self.R[j]) * self.spc_vars[i,j] )
            for j in range(len(self.R)):
                self.model.addConstr( self.spc_vars.sum('*',j) >= 1 )

        def Fix_Variables():
            for i in range(len(self.vars2fix)):
                for (u,v) in self.vars2fix[i]:
                    self.model.addConstr( self.edge_vars[u,v,i] == 1 )
        
        if self.R!=[]:
            EncodeSubpathConstraints()

        if self.vars2fix!=[]:
            Fix_Variables()

        self.model.setObjective( self.slacks.sum(), GRB.MINIMIZE )

    def print_solution(self,solution):
        opt, slack, paths = solution
        print("\n#####SOLUTION#####\n","> FD size   :",   opt,"\n> Slack sum :", slack,"\n> Weight-Slack-Path decomposition:")
        for p in paths:
            print(*p)
    
    def build_solution(self):
        paths = []
        for i in range(self.k):
            path = []
            u    = self.source
            while u != self.target:
                edges = list(filter(lambda grb_var : grb_var.X>1-TOLERANCE , self.edge_vars.select(u,'*',i) ))
                assert(len(edges)==1)
                v = head(edges[0].VarName)
                path.append(v)
                u = v
            paths.append( (self.weights[i].X, self.slacks[i].X, path[:-1]) )

        return (self.k, self.model.ObjVal, paths)

    def solve_once(self):
        logger.info(">>> Solving once")            
        self.encode()
        self.solve()
        logger.info("Gurobi solver status " + str(self.model.status))
        if self.model.status == GRB.TIME_LIMIT:
            raise utils.GRB_TimeOut("ilp.Robust.solve_once while finding feasible solution")
        elif self.model.status == GRB.INFEASIBLE:
            raise utils.GRB_Infeasible("ilp.Robust.solve_once while finding feasible solution")
        elif self.model.status == GRB.OPTIMAL:
            _,_,p = self.build_solution()
            weights,paths,slacks = list(map(lambda x : x[0], p)), list(map(lambda x : x[2], p)), list(map(lambda x : x[1], p))
            return (paths,weights,slacks,self.model.ObjVal)
            #return self.model.ObjVal
        else:
            logger.error("FATAL: solver finished with status " + str(self.model.status) + ", which is not: OPT, TIME_LIMIT, INFEASIBLE.")

    def optimize_linear(self):
        #Assumption: the initial value of self.k is sufficiently high so that the ILP solver starts with a feasible solution

        self.encode()
        self.solve()
        if self.model.status == GRB.TIME_LIMIT:
            self.final_k = self.k
            raise utils.GRB_TimeOut("ilp.Robust.optimize_linear in the feasibility stage of the optimization loop")
        elif self.model.status == GRB.INFEASIBLE:
            self.final_k = self.k
            raise utils.GRB_Infeasible("ilp.Robust.optimize_linear in th feasibility stage of the optimization loop")
        previous_slack = self.model.ObjVal

        if previous_slack == 0:
            logger.info(">>> No optimization step required, found perfect solution in init solving with " + str(self.k))
            self.final_k = self.k
            return self.final_k

        self.clear()
        self.k += 1

        logger.info(">>> Optimality, starting with " + str(self.k))
        while True: #optimality criteria: find the k for which the ratio of slacks between two consecutive iterations becomes sufficiently small

            self.encode()
            self.solve()
            logger.info("Gurobi solver status " + str(self.model.status))

            if self.model.status == GRB.TIME_LIMIT:
                self.final_k = self.k-1
                raise utils.GRB_TimeOut("ilp.Robust.optimize_linear while optimizing")
            elif self.model.status == GRB.INFEASIBLE:
                self.final_k = self.k-1
                raise utils.GRB_Infeasible("ilp.Robust.optimize_linear while optimizing")
            #assert(self.model.ObjVal <= previous_slack) not necessarily true, as the lower bound of the weight variables is 1

            if self.model.ObjVal >= previous_slack: #if we do not improve by allowing more paths we stop
                self.final_k = self.k-1
                break
            
            if self.model.ObjVal==0: #if we found a perfect solution we stop
                self.final_k = self.k
                break
            
            assert(previous_slack != 0)

            if 1-self.model.ObjVal/previous_slack < self.epsilon: # if the relative improvement is small we also stop
                self.final_k = self.k-1
                break
            
            previous_slack = self.model.ObjVal
            self.clear()
            self.k += 1
        
        return self.final_k

    
def robust(G : graph.st_DAG, epsilon=0.25, timeout=300, threads=4, path_constraints=[], vars_to_fix=[], optimize=False):

    logger.info("Robust BEGIN on graph %s", G.id)

    if len(G.edge_list)==0:
        logger.warning("Empty edge list in graph %s. Returning empty list of paths.", G.id)
        #return []
        return -1
    
    encoder = Encode_Robust(G.n, G.edge_list, G.source, G.sink, G.flow, G.w, path_constraints, vars_to_fix, epsilon, timeout, threads)

    if optimize:
        _,_,_,x = encoder.optimize_linear() #return (paths,weights,slacks,self.model.ObjVal)
        return x
    else:
        return encoder.solve_once()
        #return (paths,weights)
        #return x

def leastsquares(G : graph.st_DAG, epsilon=0.25, timeout=300, threads=4, path_constraints=[], vars_to_fix=[], optimize=False):

    logger.info("LeastSquares BEGIN on graph %s", G.id)

    if len(G.edge_list)==0:
        logger.warning("Empty edge list in graph %s. Returning empty list of paths.", G.id)
        #return []
        return -1
    
    encoder = Encode_LeastSquares(G.n, G.edge_list, G.source, G.sink, G.flow, G.w, path_constraints, vars_to_fix, epsilon, timeout, threads)

    if optimize:
        return encoder.optimize_linear()
    else:
        return encoder.solve_once()
        #return (paths,weights)