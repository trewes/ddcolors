# ddcolors
My reimplementation of van Hoeve's graph coloring algorithm using decision diagrams.
This was the work and topic of my Bachelor Thesis.
Based on the algorithm described in "W.-J. van Hoeve. Graph Coloring with Decision Diagrams. Mathematical Programming".

## Dependencies

* Exactcolors from https://github.com/heldstephan/exactcolors
  * Installation path can be specified, or a version is downloaded and compiled
* CPLEX 12.*  as the LP and MIP solver 


## Build instructions

* Download the code and enter the directory.
* Create a build directory: `mkdir build`.
* Move to the build directory: `cd build`.
* Run cmake: `cmake -DCPLEX_ROOT_DIR=</path/to/cplex> -DEXACTCOLORS_ROOT_DIR=</path/to/exactcolors> ..`.
  * The path to your CPLEX installation must be such that `</path/to/cplex>/cplex/include/ilcplex/cplex.h` exists.
  * The path to exactcolors can be given, otherwise cmake looks in the same folder that ddcolors is located in.
* Run make `make` or cmake `cmake --build .`
* Run the executable: `./ddcolors -h`.

## Structure

* In Graph.h the basic data structure is defined; and several functions for coloring heuristics 
and functions to get a permuted graph or to remove vertices.
* In DecisionDiagram, the type of a decision diagram is defined as well the functions that are executed on such a decision diagram.
* In DDColors, all subroutines from DecisionDiagram are put together to use iterative refinement or exact compilation to
 compute the chromatic number of a graph.

