import os
import random
import logging
from datetime import datetime
import safety
import ilp
import time
import utils
import argparse
import numpy as np
from collections import defaultdict

input_file  = None
output_file = None
THREADS     = None
TIMEOUT     = None
EPSILON     = None
VERBOSE     = None
MODE        = None

random.seed(73)
current_time = datetime.now()
dt_day       = current_time.strftime("%d-%m")
dt_time      = current_time.strftime("%H-%M-%S")
log_file     = "log_{}_{}.out".format(dt_day,dt_time)
logging.basicConfig(filename=log_file,format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s', datefmt='%H-%M-%S', level=logging.DEBUG)
logger       = logging.getLogger(__name__)



def demo_LQ():
    
    graphs = utils.read_graphs(input_file)
    f      = open("LQ_"+output_file+"_final.out","w")
    f.write("{}\nThreads:{}, Timeout:{}, Mode:{}\n".format(input_file,THREADS,TIMEOUT,MODE))

    t_rb_seqs_heur    = 0
    t_safe_seqs_heur  = 0
    t_ilp_seqs_heur   = 0
    solved_seqs_heur  = 0
    fixed_vars_s      = 0

    for G in graphs:

        print("__demo_final__ Running on " + str(G.id) + "/" + str(len(graphs)) + " with n=" + str(G.n) + ", m=" + str(G.m) + " and w=" + str(G.w))

        if utils.is_0_flow_everywhere(G):
            logger.info("Found 0 flow everywhere, skipping graph " + str(G.id))
            continue

        f.write("#Graph {}\n{}, {}, {}\n".format(G.id,G.n,G.m,G.w))
        
        logger.info("   Starting demo on graph \'{}\' with id {}".format(input_file,G.id))

        obj1=None
        obj2=None

        #Vanilla
        try:
            start  = time.time()
            obj1   = ilp.leastsquares(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS)
            end    = time.time()
            t_rb_default   = end-start
            solved_default = True
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' in vanilla mode.".format(e,G.id,input_file))
            t_rb_default   = 0
            solved_default = False
        except utils.GRB_Infeasible as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' in vanilla mode.".format(e,G.id,input_file))
            t_rb_default   = 0
            solved_default = False
        logger.info("\tvanilla: {}, {}".format(solved_default,t_rb_default) )


        try:
            #X = set( filter( lambda edge : G.flow[edge] >= np.percentile(list(G.flow.values()), 25), G.edge_list) )
            X = set(G.edge_list)
            t0        = time.time()

            safe_seqs = safety.maximal_safe_sequences_via_dominators(G, X)

            longest_safe_sequence = dict()
            len_of_longest_ss     = defaultdict(lambda:0)
            for i in range(len(safe_seqs)):
                safe_seq = safe_seqs[i]
                length = len(safe_seq)
                for edge in safe_seq:
                    if edge not in longest_safe_sequence:
                        longest_safe_sequence[edge] = i
                        len_of_longest_ss[edge]     = length
                    elif len(safe_seqs[longest_safe_sequence[edge]]) < length:
                        longest_safe_sequence[edge] = i
                        len_of_longest_ss[edge]     = length

            _, edge_antichain = utils.max_edge_antichain(G, get_antichain=True, weight_function=len_of_longest_ss)
            sequences_to_fix  = list(map(lambda edge : safe_seqs[longest_safe_sequence[edge]] , edge_antichain))

            t1   = time.time()
            obj2 = ilp.leastsquares(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, vars_to_fix=sequences_to_fix)
            t2   = time.time()

            t_rb_seqs_heur   = t2-t0
            t_safe_seqs_heur = t1 - t0
            t_ilp_seqs_heur  = t2-t1
            solved_seqs_heur = True
            fixed_vars_s     = sum(map(lambda sequence : len(sequence), sequences_to_fix))
        
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subsequence constraints.".format(e,G.id,input_file))
            t_rb_seqs_heur = t_ilp_seqs_heur = 0
            solved_seqs_heur                 = False
            fixed_vars_s = t_safe_seqs_heur  = 0
        except utils.GRB_Infeasible as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subsequence constraints.".format(e,G.id,input_file))
            t_rb_seqs_heur = t_ilp_seqs_heur = 0
            solved_seqs_heur                 = False
            fixed_vars_s = t_safe_seqs_heur  = 0
        logger.info("\tfixing safe sequences: {}, {} , {}, {}, {}".format(solved_seqs_heur, t_rb_seqs_heur, t_safe_seqs_heur, t_ilp_seqs_heur, fixed_vars_s) )

        if obj1!=None:
            if obj2==None:
                print("Safety lost against vanilla LQ")
                logger.info("\t\t: LQ: Safety lost against vanilla on graph " + str(G.id))
            else:
                if (obj1!=obj2):
                    #print(obj1,obj2)
                    print("\nPROBLEM\n")

        f.write("solved default              : {}\n".format(solved_default              ))
        f.write("total time default          : {}\n".format('%.6f' % t_rb_default       ))
        f.write("solved sequences heur       : {}\n".format(solved_seqs_heur            ))
        f.write("total time sequences heur   : {}\n".format('%.6f' % t_rb_seqs_heur     ))
        f.write("preprocess sequences heur   : {}\n".format('%.6f' % t_safe_seqs_heur   ))
        f.write("ilp time seqs heur          : {}\n".format('%.6f' % t_ilp_seqs_heur    ))
        f.write("fixed vars seqs             : {}\n".format(fixed_vars_s                ))

    f.close()
    return



def demo_RB():
    
    graphs = utils.read_graphs(input_file)
    f      = open("RB_"+output_file+"_final.out","w")
    f.write("{}\nThreads:{}, Timeout:{}, Mode:{}\n".format(input_file,THREADS,TIMEOUT,MODE))

    t_rb_seqs_heur    = 0
    t_safe_seqs_heur  = 0
    t_ilp_seqs_heur   = 0
    solved_seqs_heur  = 0
    fixed_vars_s      = 0

    for G in graphs:

        print("__demo_final__ Running on " + str(G.id) + " out of " + str(len(graphs)) + " with n=" + str(G.n) + ", m=" + str(G.m) + " and w=" + str(G.w))

        if utils.is_0_flow_everywhere(G):
            logger.info("Found 0 flow everywhere, skipping graph " + str(G.id))
            continue

        f.write("#Graph {}\n{}, {}, {}\n".format(G.id,G.n,G.m,G.w))
        
        logger.info("   Starting demo on graph \'{}\' with id {}".format(input_file,G.id))

        obj1 = None
        obj2 = None

        #Vanilla
        try:
            start  = time.time()
            obj1   = ilp.robust(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS)
            end    = time.time()
            t_rb_default   = end-start
            solved_default = True
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' in vanilla mode.".format(e,G.id,input_file))
            t_rb_default   = 0
            solved_default = False
        except utils.GRB_Infeasible as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' in vanilla mode.".format(e,G.id,input_file))
            t_rb_default   = 0
            solved_default = False
        logger.info("\tvanilla: {}, {}".format(solved_default,t_rb_default) )


        try:
            #X = set( filter( lambda edge : G.flow[edge] >= np.percentile(list(G.flow.values()), 25), G.edge_list) )
            X = set(G.edge_list)
            t0        = time.time()

            safe_seqs = safety.maximal_safe_sequences_via_dominators(G, X)

            longest_safe_sequence = dict()
            len_of_longest_ss     = defaultdict(lambda:0)
            for i in range(len(safe_seqs)):
                safe_seq = safe_seqs[i]
                length = len(safe_seq)
                for edge in safe_seq:
                    if edge not in longest_safe_sequence:
                        longest_safe_sequence[edge] = i
                        len_of_longest_ss[edge]     = length
                    elif len(safe_seqs[longest_safe_sequence[edge]]) < length:
                        longest_safe_sequence[edge] = i
                        len_of_longest_ss[edge]     = length

            _, edge_antichain = utils.max_edge_antichain(G, get_antichain=True, weight_function=len_of_longest_ss)
            
            sequences_to_fix  = list(map(lambda edge : safe_seqs[longest_safe_sequence[edge]] , edge_antichain))

            t1   = time.time()
            obj2 = ilp.robust(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, vars_to_fix=sequences_to_fix)
            t2   = time.time()

            t_rb_seqs_heur   = t2-t0
            t_safe_seqs_heur = t1-t0
            t_ilp_seqs_heur  = t2-t1
            solved_seqs_heur = True
            fixed_vars_s     = sum(map(lambda sequence : len(sequence), sequences_to_fix))
        
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subsequence constraints.".format(e,G.id,input_file))
            t_rb_seqs_heur = t_ilp_seqs_heur = 0
            solved_seqs_heur                 = False
            fixed_vars_s = t_safe_seqs_heur  = 0
        except utils.GRB_Infeasible as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subsequence constraints.".format(e,G.id,input_file))
            t_rb_seqs_heur = t_ilp_seqs_heur = 0
            solved_seqs_heur                 = False
            fixed_vars_s = t_safe_seqs_heur  = 0
        logger.info("\tfixing safe sequences: {}, {} , {}, {}, {}".format(solved_seqs_heur, t_rb_seqs_heur, t_safe_seqs_heur, t_ilp_seqs_heur, fixed_vars_s) )

        if obj1!=None:
            if obj2==None:
                print("Safety lost against vanilla RB")
                logger.info("\t\t: RB: Safety lost against vanilla on graph " + str(G.id))
            else:
                if (obj1!=obj2):
                    #print(obj1,obj2)
                    print("\nPROBLEM\n")

        f.write("solved default              : {}\n".format(solved_default              ))
        f.write("total time default          : {}\n".format('%.6f' % t_rb_default       ))
        f.write("solved sequences heur       : {}\n".format(solved_seqs_heur            ))
        f.write("total time sequences heur   : {}\n".format('%.6f' % t_rb_seqs_heur     ))
        f.write("preprocess sequences heur   : {}\n".format('%.6f' % t_safe_seqs_heur   ))
        f.write("ilp time seqs heur          : {}\n".format('%.6f' % t_ilp_seqs_heur    ))
        f.write("fixed vars seqs             : {}\n".format(fixed_vars_s                ))

    f.close()
    return


def demo_optimize_RB():
    
    graphs = utils.read_graphs(input_file)
    f      = open("OPT_RB_"+output_file+"_final.out","w")
    f.write("{}\nThreads:{}, Timeout:{}, Mode:{}, Epsilon:{}\n".format(input_file,THREADS,TIMEOUT,MODE,EPSILON))

    t_rb_seqs_heur    = 0
    t_safe_seqs_heur  = 0
    t_ilp_seqs_heur   = 0
    solved_seqs_heur  = 0
    fixed_vars_s      = 0
    w_van             = 0
    w_paths           = 0
    w_seqs            = 0

    for G in graphs:

        print("__demo_final__ Running on " + str(G.id) + "/" + str(len(graphs)) + " with n=" + str(G.n) + ", m=" + str(G.m) + " and w=" + str(G.w))

        if utils.is_0_flow_everywhere(G):
            logger.info("Found 0 flow everywhere, skipping graph " + str(G.id))
            continue

        f.write("#Graph {}\n{}, {}, {}\n".format(G.id,G.n,G.m,G.w))
        
        logger.info("   Starting demo on graph \'{}\' with id {}".format(input_file,G.id))

        #Vanilla
        try:
            start  = time.time()
            w_van  = ilp.robust(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, optimize=True)
            end    = time.time()
            t_rb_default   = end-start
            solved_default = True
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' in vanilla mode.".format(e,G.id,input_file))
            t_rb_default   = 0
            solved_default = False
            w_van          = -1
        except utils.GRB_Infeasible as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' in vanilla mode.".format(e,G.id,input_file))
            t_rb_default   = 0
            solved_default = False
            w_van          = -1
        logger.info("\tvanilla: {}, {}, {}".format(solved_default,t_rb_default, w_van) )


        try:
            start        = time.time()

            safe_seqs    = safety.maximal_safe_sequences(G, G.edge_list)

            longest_safe_sequence = dict()
            len_of_longest_ss     = dict()
            for i in range(len(safe_seqs)):
                safe_seq = safe_seqs[i]
                length = len(safe_seq)
                for edge in safe_seq:
                    if edge not in longest_safe_sequence:
                        longest_safe_sequence[edge] = i
                        len_of_longest_ss[edge]     = length
                    elif len(safe_seqs[longest_safe_sequence[edge]]) < length:
                        longest_safe_sequence[edge] = i
                        len_of_longest_ss[edge]     = length

            _, edge_antichain = utils.max_edge_antichain(G, get_antichain=True, weight_function=len_of_longest_ss)
            sequences_to_fix  = list(map(lambda edge : safe_seqs[longest_safe_sequence[edge]] , edge_antichain))

            time_safety   = time.time()

            time_lp_start = time.time()
            w_seqs        = ilp.robust(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, vars_to_fix=sequences_to_fix, optimize=True)
            time_lp_end   = time.time()
            
            end           = time.time()

            t_rb_seqs_heur   = end-start
            t_safe_seqs_heur = time_safety - start
            t_ilp_seqs_heur  = time_lp_end-time_lp_start
            solved_seqs_heur = True
            fixed_vars_s     = sum(map(lambda sequence : len(sequence), sequences_to_fix))
        
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subsequence constraints.".format(e,G.id,input_file))
            t_rb_seqs_heur = t_ilp_seqs_heur = 0
            solved_seqs_heur                 = False
            w_seqs                           = -1
            fixed_vars_s = t_safe_seqs_heur  = 0
        except utils.GRB_Infeasible as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subsequence constraints.".format(e,G.id,input_file))
            t_rb_seqs_heur = t_ilp_seqs_heur = 0
            solved_seqs_heur                 = False
            w_seqs                           = -1
            fixed_vars_s = t_safe_seqs_heur  = 0
        logger.info("\tfixing safe sequences: {}, {} , {}, {}, {}, {}".format(solved_seqs_heur, t_rb_seqs_heur, t_safe_seqs_heur, t_ilp_seqs_heur, fixed_vars_s, w_seqs) )


        f.write("solved default              : {}\n".format(solved_default              ))
        f.write("total time default          : {}\n".format('%.6f' % t_rb_default       ))
        f.write("solved sequences heur       : {}\n".format(solved_seqs_heur            ))
        f.write("total time sequences heur   : {}\n".format('%.6f' % t_rb_seqs_heur     ))
        f.write("preprocess sequences heur   : {}\n".format('%.6f' % t_safe_seqs_heur   ))
        f.write("ilp time seqs heur          : {}\n".format('%.6f' % t_ilp_seqs_heur    ))
        f.write("fixed vars seqs             : {}\n".format(fixed_vars_s                ))
        f.write("final width default         : {}\n".format(w_van                       ))
        f.write("final width sequences       : {}\n".format(w_seqs                      ))

    f.close()
    return


def skeleton():
    pass

def main():
    
    global input_file
    global output_file
    global log_file
    global THREADS
    global EPSILON
    global TIMEOUT
    global VERBOSE
    global MODE

    parser = argparse.ArgumentParser(description='Process inputs.')

    parser.add_argument('-i', '--input'  , required=True             , help='Input file path'                                                               )
    parser.add_argument('-t', '--threads', type=int                  , help='Number of threads (default: 4)'                                  , default=4   )
    parser.add_argument('-g', '--timeout', type=int                  , help='Timeout in seconds (default: 300)'                               , default=300 )
    parser.add_argument('-e', '--epsilon', type=float                , help='Relative optima improvement for Gurobi (must be between 0 and 1)', default=0.25)
    parser.add_argument('-c', '--clear'  , action='store_true'       , help='Enable clear mode'                                                             )
    parser.add_argument('-v', '--verbose', action='store_true'       , help='Enable verbose mode'                                                           )
    parser.add_argument('-m', '--mode'   , choices=['0','1','2','3'] , help='Optimization mode'                                                             )

    args = parser.parse_args()

    input_file  = args.input
    output_file = '{}_{}_{}'.format(input_file.replace("/","_"),dt_day,dt_time)
    THREADS     = args.threads
    TIMEOUT     = args.timeout
    EPSILON     = args.epsilon
    VERBOSE     = args.verbose
    MODE        = args.mode

    print(f"Input file : {input_file}")
    print(f"Num threads: {THREADS}")
    print(f"GRB epsilon: {EPSILON}")
    print(f"Timeout    : {TIMEOUT} seconds")
    print(f"Verbose    : {VERBOSE}")
    print(f"Mode       : {MODE}")
    print(f"Clear      : {args.clear}")

    if MODE == '0':
        demo_RB()
    elif MODE == '1':
        demo_LQ()
    elif MODE == '2':
        demo_optimize_RB()
    elif MODE== '3':
        skeleton()
    else:
        print("ERROR: bad mode to execute - must be 0 (for Robust) and 1 (for LeastSquares).")

    #Cleaner
    if args.clear:
        if os.path.exists(log_file):
            os.remove(log_file)
            print(f"File '{log_file}' has been deleted.")
        else:
            print(f"File '{log_file}' does not exist.")

    print("__main__ completed")


if __name__ == "__main__":
    main()

