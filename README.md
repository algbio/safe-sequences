# Safe-sequences

**04/11/2024**

## About safe sequences

The driving motivation for this project is to extend the well known notion of safe paths to safe sequences and provide implementations to efficiently compute safe sequences via dominators. [Safety](https://link.springer.com/chapter/10.1007/978-3-319-31957-5_11) has found applications to numerous bioinformatic problems, namely the genome-guided RNA transcript assembly problem, where this current implementation focuses on. More specifically, given a set of RNA-seq reads, a directed acyclic graph (DAG) is constructed from their alignments to a reference genome. The graph nodes correspond to e.g exons, the arcs correspond to reads overlapping two consecutive exons, and the node or arc weights corresponding their read coverage. The RNA transcripts then correspond to a set of source-to-sink weighted paths in the DAG that best explain the nodes, arcs and their weights, under various definitions of optimality. For example, under the assumption of perfect data, the weights on the arcs of the graph induce a flow, and thus the Minimum Flow Decomposition problem (MFD) becomes a natural and suitable abstraction.

Integer Linear Programming (ILP) is an ubiquous framework to tackle NP-hard problems. As such, one can rely on ILP solvers to solve MFD and many of its variants in an attempt to find optimal sets of paths that explain the graph's structure and flow. Due to the inherent complexity of MFD and its variants, even optimized solvers (e.g., Gurobi) cannot compute solutions within a practical timeframe. Therefore, incorporating safety information into the Linear Program may help accelerate the solver's execution. The current implementation is tailored to experimentally assess these speed ups.

## Requirements

To solve linear programs we use [Gurobi](https://www.gurobi.com/). We recommend checking [this](https://www.gurobi.com/academia/academic-program-and-licenses/).
To compute maximum weight edge antichains we use a classical reduction to the Minimum Flow problem and then use [NetworkX](https://networkx.org/) to find these antichains for us. All of the above can be easily installed in Python.

## Usage

### Input format

Every graph has a header line of the form "#Graph" followed by a space and an identifier of the graph, i.e., a string. Followed by the header, in a new line, the number of nodes of the graph is written. Then, every next line should be of the form "u v f(u,v)", where u,v are nodes of the graph represented by natural numbers and f(u,v) is the weight on the arc (u,v) (so there should be as many of these lines as arcs in the graph). The values f(u,v) should be integral but not necessarily need to induce a flow.

```text
#Graph id1
n_1
u v f(u,v)
v w f(v,w)
...
u w f(u,w)
#Graph id2
n_2
a b f(a,b)
b c f(b,c)
...
g l f(g,l)
...
#Graph id8
n_100
...
```

Every input file should contain at least one graph. The graph must be a DAG without parallel arcs.
We do not give any correctness guarantees when different graphs carry the same identifiers.
Check the example/test.graph file in the repository for a concrete example.

### Example
```bash
    python3 main.py -i example/test.graph -m 0 -c
```

### Options:
```text
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file path
  -t THREADS, --threads THREADS
                        Number of threads for Gurobi (default: 4)
  -g TIMEOUT, --timeout TIMEOUT
                        Timeout for Gurobi in seconds (default: 300)
  -e EPSILON, --epsilon EPSILON
                        Relative optima improvement for Gurobi in two consecutive iterations; must be between 0 and 1 (default: 0.25)
  -c, --clear           Clears log file before exiting
  -m {0,1,2,3}, --mode {0,1,2,3}
                        Mode to run. 0: demo used in the paper for MinPathError; 1: demo used in the paper for LeastSquares; 2:
                        same as 1 but actually optimizes on the solution size and the cumulative errors; 3: skeleton function for
                        safety (utilize as you see fit).
```

### Output and results analysis'
When running with modes 0, 1, and 2, some metrics of interest are written to a file (e.g., Gurobi's running time equipped with and without safety information, maximum edge antichain of the graph, and more). To produce the LateX tables use the `stats.py` module. In the stats module please note that the `width-ranges` parameter may need adjustment.

## Contact

Please contact the authors for any problem related to the code, namely errors and suggestions to improve.