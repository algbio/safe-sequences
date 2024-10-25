import os
import random
import logging
from datetime import datetime
import safety
import ilp
import time
import utils
import argparse


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

    t_lq_paths_heur   = 0
    t_safe_paths_heur = 0
    t_ilp_paths_heur  = 0
    solved_paths_heur = 0
    fixed_vars_p      = 0
    t_lq_seqs_heur    = 0
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

        #Vanilla
        try:
            start  = time.time()
            ilp.leastsquares(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS)
            end    = time.time()
            t_rb_default   = end-start
            solved_default = True
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' in vanilla mode.".format(e,G.id,input_file))
            t_rb_default   = 0
            solved_default = False
        logger.info("\tvanilla: {}, {}".format(solved_default,t_rb_default) )


        #Fixing heuristic
        try:
            start        = time.time()

            safe_paths   = safety.safe_paths(G)

            longest_safe_path = dict()
            len_of_longest_sp = dict()
            for i in range(len(safe_paths)):
                safe_path = safe_paths[i]
                length    = len(safe_path)
                for edge in safe_path:
                    if edge not in longest_safe_path:
                        longest_safe_path[edge] = i
                        len_of_longest_sp[edge] = length
                    elif len(safe_paths[longest_safe_path[edge]]) < length:
                        longest_safe_path[edge] = i
                        len_of_longest_sp[edge] = length

            _, edge_antichain = utils.max_edge_antichain(G, get_antichain=True, weight_function=len_of_longest_sp)
            paths_to_fix      = list(map(lambda edge : safe_paths[longest_safe_path[edge]] , edge_antichain))

            time_safety   = time.time()

            time_lp_start = time.time()
            ilp.leastsquares(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, vars_to_fix=paths_to_fix)
            time_lp_end   = time.time()
            
            end           = time.time()

            t_lq_paths_heur   = end-start
            t_safe_paths_heur = time_safety - start
            t_ilp_paths_heur  = time_lp_end-time_lp_start
            solved_paths_heur = True
            fixed_vars_p      = sum(map(lambda path : len(path), paths_to_fix))
        
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subpath constraints.".format(e,G.id,input_file))
            t_lq_paths_heur = t_ilp_paths_heur = 0
            solved_paths_heur                  = False
            fixed_vars_p = t_safe_paths_heur   = 0
        logger.info("\tfixing safe paths: {}, {} , {}, {}, {}".format(solved_paths_heur, t_lq_paths_heur, t_safe_paths_heur, t_ilp_paths_heur, fixed_vars_p) )

        try:
            start        = time.time()

            safe_seqs    = safety.safe_sequences(G)

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
            ilp.leastsquares(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, vars_to_fix=sequences_to_fix)
            time_lp_end   = time.time()
            
            end           = time.time()

            t_lq_seqs_heur   = end-start
            t_safe_seqs_heur = time_safety - start
            t_ilp_seqs_heur  = time_lp_end-time_lp_start
            solved_seqs_heur = True
            fixed_vars_s     = sum(map(lambda sequence : len(sequence), sequences_to_fix))
        
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subsequence constraints.".format(e,G.id,input_file))
            t_lq_seqs_heur = t_ilp_seqs_heur = 0
            solved_seqs_heur                 = False
            fixed_vars_s = t_safe_seqs_heur  = 0
        logger.info("\tfixing safe sequences: {}, {} , {}, {}, {}".format(solved_seqs_heur, t_lq_seqs_heur, t_safe_seqs_heur, t_ilp_seqs_heur, fixed_vars_s) )


        f.write("solved default              : {}\n".format(solved_default              ))
        f.write("total time default          : {}\n".format('%.6f' % t_rb_default       ))

        f.write("solved paths heur           : {}\n".format(solved_paths_heur           ))
        f.write("total time paths heur       : {}\n".format('%.6f' % t_lq_paths_heur    ))
        f.write("solved sequences heur       : {}\n".format(solved_seqs_heur            ))
        f.write("total time sequences heur   : {}\n".format('%.6f' % t_lq_seqs_heur     ))
        f.write("preprocess paths heur       : {}\n".format('%.6f' % t_safe_paths_heur  ))
        f.write("preprocess sequences heur   : {}\n".format('%.6f' % t_safe_seqs_heur   ))
        f.write("ilp time paths heur         : {}\n".format('%.6f' % t_ilp_paths_heur   ))
        f.write("ilp time seqs heur          : {}\n".format('%.6f' % t_ilp_seqs_heur    ))
        f.write("fixed vars paths            : {}\n".format(fixed_vars_p                ))
        f.write("fixed vars seqs             : {}\n".format(fixed_vars_s                ))

    f.close()
    return



def demo_RB():
    
    graphs = utils.read_graphs(input_file)
    f      = open("RB_"+output_file+"_final.out","w")
    f.write("{}\nThreads:{}, Timeout:{}, Mode:{}\n".format(input_file,THREADS,TIMEOUT,MODE))

    t_rb_paths_heur   = 0
    t_safe_paths_heur = 0
    t_ilp_paths_heur  = 0
    solved_paths_heur = 0
    fixed_vars_p      = 0
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

        #Vanilla
        try:
            start  = time.time()
            ilp.robust(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS)
            end    = time.time()
            t_rb_default   = end-start
            solved_default = True
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' in vanilla mode.".format(e,G.id,input_file))
            t_rb_default   = 0
            solved_default = False
        logger.info("\tvanilla: {}, {}".format(solved_default,t_rb_default) )


        #Fixing heuristic
        try:
            start        = time.time()

            safe_paths   = safety.safe_paths(G)

            longest_safe_path = dict()
            len_of_longest_sp = dict()
            for i in range(len(safe_paths)):
                safe_path = safe_paths[i]
                length    = len(safe_path)
                for edge in safe_path:
                    if edge not in longest_safe_path:
                        longest_safe_path[edge] = i
                        len_of_longest_sp[edge] = length
                    elif len(safe_paths[longest_safe_path[edge]]) < length:
                        longest_safe_path[edge] = i
                        len_of_longest_sp[edge] = length

            _, edge_antichain = utils.max_edge_antichain(G, get_antichain=True, weight_function=len_of_longest_sp)
            paths_to_fix      = list(map(lambda edge : safe_paths[longest_safe_path[edge]] , edge_antichain))

            time_safety   = time.time()

            time_lp_start = time.time()
            ilp.robust(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, vars_to_fix=paths_to_fix)
            time_lp_end   = time.time()
            
            end           = time.time()

            t_rb_paths_heur   = end-start
            t_safe_paths_heur = time_safety - start
            t_ilp_paths_heur  = time_lp_end-time_lp_start
            solved_paths_heur = True
            fixed_vars_p      = sum(map(lambda path : len(path), paths_to_fix))
        
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subpath constraints.".format(e,G.id,input_file))
            t_rb_paths_heur = t_ilp_paths_heur = 0
            solved_paths_heur                  = False
            fixed_vars_p = t_safe_paths_heur  = 0
        logger.info("\tfixing safe paths: {}, {} , {}, {}, {}".format(solved_paths_heur, t_rb_paths_heur, t_safe_paths_heur, t_ilp_paths_heur, fixed_vars_p) )

        try:
            start        = time.time()

            safe_seqs    = safety.safe_sequences(G)

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
            ilp.robust(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, vars_to_fix=sequences_to_fix)
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
            fixed_vars_s = t_safe_seqs_heur  = 0
        logger.info("\tfixing safe sequences: {}, {} , {}, {}, {}".format(solved_seqs_heur, t_rb_seqs_heur, t_safe_seqs_heur, t_ilp_seqs_heur, fixed_vars_s) )


        f.write("solved default              : {}\n".format(solved_default              ))
        f.write("total time default          : {}\n".format('%.6f' % t_rb_default       ))

        f.write("solved paths heur           : {}\n".format(solved_paths_heur           ))
        f.write("total time paths heur       : {}\n".format('%.6f' % t_rb_paths_heur    ))
        f.write("solved sequences heur       : {}\n".format(solved_seqs_heur            ))
        f.write("total time sequences heur   : {}\n".format('%.6f' % t_rb_seqs_heur     ))
        f.write("preprocess paths heur       : {}\n".format('%.6f' % t_safe_paths_heur  ))
        f.write("preprocess sequences heur   : {}\n".format('%.6f' % t_safe_seqs_heur   ))
        f.write("ilp time paths heur         : {}\n".format('%.6f' % t_ilp_paths_heur   ))
        f.write("ilp time seqs heur          : {}\n".format('%.6f' % t_ilp_seqs_heur    ))
        f.write("fixed vars paths            : {}\n".format(fixed_vars_p                ))
        f.write("fixed vars seqs             : {}\n".format(fixed_vars_s                ))

    f.close()
    return


def demo_optimize_RB():
    
    graphs = utils.read_graphs(input_file)
    f      = open("OPT_RB_"+output_file+"_final.out","w")
    f.write("{}\nThreads:{}, Timeout:{}, Mode:{}, Epsilon:{}\n".format(input_file,THREADS,TIMEOUT,MODE,EPSILON))

    t_rb_paths_heur   = 0
    t_safe_paths_heur = 0
    t_ilp_paths_heur  = 0
    solved_paths_heur = 0
    fixed_vars_p      = 0
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
        logger.info("\tvanilla: {}, {}, {}".format(solved_default,t_rb_default, w_van) )


        #Fixing heuristic
        try:
            start        = time.time()

            safe_paths   = safety.safe_paths(G)

            longest_safe_path = dict()
            len_of_longest_sp = dict()
            for i in range(len(safe_paths)):
                safe_path = safe_paths[i]
                length    = len(safe_path)
                for edge in safe_path:
                    if edge not in longest_safe_path:
                        longest_safe_path[edge] = i
                        len_of_longest_sp[edge] = length
                    elif len(safe_paths[longest_safe_path[edge]]) < length:
                        longest_safe_path[edge] = i
                        len_of_longest_sp[edge] = length

            _, edge_antichain = utils.max_edge_antichain(G, get_antichain=True, weight_function=len_of_longest_sp)
            paths_to_fix      = list(map(lambda edge : safe_paths[longest_safe_path[edge]] , edge_antichain))

            time_safety   = time.time()

            time_lp_start = time.time()
            w_paths       = ilp.robust(G, epsilon=EPSILON, timeout=TIMEOUT, threads=THREADS, vars_to_fix=paths_to_fix, optimize=True)
            time_lp_end   = time.time()
            
            end           = time.time()

            t_rb_paths_heur   = end-start
            t_safe_paths_heur = time_safety - start
            t_ilp_paths_heur  = time_lp_end-time_lp_start
            solved_paths_heur = True
            fixed_vars_p      = sum(map(lambda path : len(path), paths_to_fix))
        
        except utils.GRB_TimeOut as e:
            logger.info("\t{}. Graph {} in dataset \'{}\' with fixing subpath constraints.".format(e,G.id,input_file))
            t_rb_paths_heur = t_ilp_paths_heur = 0
            solved_paths_heur                  = False
            w_paths                            = -1
            fixed_vars_p = t_safe_paths_heur   = 0
        logger.info("\tfixing safe paths: {}, {} , {}, {}, {}, {}".format(solved_paths_heur, t_rb_paths_heur, t_safe_paths_heur, t_ilp_paths_heur, fixed_vars_p, w_paths) )

        try:
            start        = time.time()

            safe_seqs    = safety.safe_sequences(G)

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
        logger.info("\tfixing safe sequences: {}, {} , {}, {}, {}, {}".format(solved_seqs_heur, t_rb_seqs_heur, t_safe_seqs_heur, t_ilp_seqs_heur, fixed_vars_s, w_seqs) )


        f.write("solved default              : {}\n".format(solved_default              ))
        f.write("total time default          : {}\n".format('%.6f' % t_rb_default       ))

        f.write("solved paths heur           : {}\n".format(solved_paths_heur           ))
        f.write("total time paths heur       : {}\n".format('%.6f' % t_rb_paths_heur    ))
        f.write("solved sequences heur       : {}\n".format(solved_seqs_heur            ))
        f.write("total time sequences heur   : {}\n".format('%.6f' % t_rb_seqs_heur     ))
        f.write("preprocess paths heur       : {}\n".format('%.6f' % t_safe_paths_heur  ))
        f.write("preprocess sequences heur   : {}\n".format('%.6f' % t_safe_seqs_heur   ))
        f.write("ilp time paths heur         : {}\n".format('%.6f' % t_ilp_paths_heur   ))
        f.write("ilp time seqs heur          : {}\n".format('%.6f' % t_ilp_seqs_heur    ))
        f.write("fixed vars paths            : {}\n".format(fixed_vars_p                ))
        f.write("fixed vars seqs             : {}\n".format(fixed_vars_s                ))

        f.write("final width default         : {}\n".format(w_van                       ))
        f.write("final width paths           : {}\n".format(w_paths                     ))
        f.write("final width sequences       : {}\n".format(w_seqs                      ))

    f.close()
    return


def compute_safety(): # A skeleton function
    graphs = utils.read_graphs(input_file)

    for G in graphs:
        _ = safety.safe_paths(G)
        _ = safety.safe_sequences(G)
        
        # Use safety information here

    return


def main():
    
    global input_file
    global output_file
    global log_file
    global THREADS
    global EPSILON
    global TIMEOUT
    global MODE

    parser = argparse.ArgumentParser(description='Process inputs.')

    parser.add_argument('-i', '--input'  , required=True             , help='Input file path'                                                                                                        )
    parser.add_argument('-t', '--threads', type=int                  , help='Number of threads for Gurobi (default: 4)'                                                                , default=4   )
    parser.add_argument('-g', '--timeout', type=int                  , help='Timeout for Gurobi in seconds (default: 300)'                                                             , default=300 )
    parser.add_argument('-e', '--epsilon', type=float                , help='Relative optima improvement for Gurobi in consecutive iterations; must be between 0 and 1 (default: 0.25)', default=0.25)
    parser.add_argument('-c', '--clear'  , action='store_true'       , help='Clears log file before exiting'                                                                                         )
    parser.add_argument('-m', '--mode'   , choices=['0','1','2','3'] , help='Mode to run. 0: demo used in the paper for MinPathError; 1: demo used in the paper for LeastSquares; 2: same as 1 but actually optimizes on the solution size and the cumulative errors; 3: skeleton function for safety (utilize as you see fit).')

    args = parser.parse_args()

    input_file  = args.input
    output_file = '{}_{}_{}'.format(input_file.replace("/","_"),dt_day,dt_time)
    THREADS     = args.threads
    TIMEOUT     = args.timeout
    EPSILON     = args.epsilon
    MODE        = args.mode

    print(f"Input file : {input_file}")
    print(f"Num threads: {THREADS}")
    print(f"GRB epsilon: {EPSILON}")
    print(f"Timeout    : {TIMEOUT} seconds")
    print(f"Mode       : {MODE}")
    print(f"Clear      : {args.clear}")

    if MODE == '0':
        demo_RB()
    elif MODE == '1':
        demo_LQ()
    elif MODE == '2':
        demo_optimize_RB()
    elif MODE == '3':
        compute_safety()
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

