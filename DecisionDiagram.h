#ifndef DDCOLORS_DECISIONDIAGRAM_H
#define DDCOLORS_DECISIONDIAGRAM_H

/*
 * DecisionDiagram.h
 * Purpose: Implementation of classes and data structure to represent a decision diagram and
 * perform various function with it as needed by the exact compilation and iterative refinement framework.
 * Most notably implements the nodes as a data struct and the decision diagram as a vector of vector of nodes.
 */

#include <vector>
#include <set>
#include <iostream>
#include <tuple>
#include <algorithm>
#include <map>
#include <cmath>

#include <iomanip>
#include <cfenv>
//#pragma STDC FENV_ACCESS ON //unused/has no effect

#include "Graph.h"
#include "unordered_map"

extern "C" {
#include "lp.h"            //needed for the interface to CPLEX to build and solve linear and integer programs
#include <ilcplex/cplex.h> //needed for CPX_INFBOUND
}

//In the path decomposition algorithms, when subtracting from the residual flow value floating point errors may occur
//set a limit to how small the residual flow can get to avoid such errors
#define double_eps std::pow(10, -8)

//allows to print sets and vectors in a nice format via std::cout
template<typename T>
std::ostream &operator<<(std::ostream &s, const std::set<T> &set);
template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec);

typedef int LayerIndex;
typedef int NodeIndex;
typedef int EdgeIndex;
typedef unsigned int Vertex;
typedef std::set<Vertex> StateInfo;

/*
 * Node
 * Purpose : Store the required data of a node needed in the algorithms performed on a decision diagram
 *
 * layer : the layer which the node belongs to
 * index : the index of the node in its layer
 * state_info : the set of eligible vertices
 * one_arc : index of node in next layer that the 1-arc is connected with
 * zero_arc : index of node in next layer that the 0-arc is connected with
 * one_arc_flow : the amount of flow on the 1-arc after having computed a flow solution
 * zero_arc_flow : the amount of flow on the 0-arc after having computed a flow solution
 *
 * Node() : construct an uninitialised node
 * Node(layer, index, state_info) : constructs a node initialising layer, index and state_info
 */

struct Node{
    Node();
    Node(LayerIndex layer, NodeIndex index, StateInfo state_info = {});

    LayerIndex layer = -1;
    NodeIndex index = -1;
    StateInfo state_info;
    NodeIndex one_arc = -1;
    NodeIndex zero_arc = -1;

    double one_arc_flow = 0;
    double zero_arc_flow = 0;

    bool equivalent(const Node &other) const;
};

/*
 * DecisionDiagram
 * Purpose : represent and store the data of a decision diagram to a graph with n vertices, this is simply
 * a vector of n+1 layers, each layer being a dynamically sized vector of nodes
 *
 * num_vars : returns the number of decision variables belonging to the problem, i.e. num layers - 1
 * num_nodes : returns the number of nodes of the decision diagram
 * num_arcs : returns the number of arcs of the decision diagram
 * get_width : returns the size of the largest layer of the decision diagram
 * print_decision_diagram : prints the decision diagram layer by layer
 *
 * initial_decision_diagram : builds the initial diagram with a single node in each layer and the 1-arc and 0-arc
 *                            both going into the next node. Has n+1 layers and nodes where n = |V(G)|
 * exact_decision_diagram : builds the exact decision diagram using top-down compilation and checking for equivalence
 *                          of nodes using the state information. Terminates if more than one million nodes are required
 */

typedef std::vector<std::vector<Node> > DecisionDiagram;

unsigned int num_vars(const DecisionDiagram &dd);
unsigned int num_nodes(const DecisionDiagram &dd);
unsigned int num_arcs(const DecisionDiagram &dd);
unsigned int get_width(const DecisionDiagram &dd);

void print_decision_diagram(const DecisionDiagram &dd, bool use_tag = false);

DecisionDiagram initial_decision_diagram(const Graph &g);

DecisionDiagram exact_decision_diagram(const Graph &g);
DecisionDiagram exact_decision_diagram(const Graph &g, const NeighborList &neighbors);

/*
 * PathLabelConflict
 * Purpose : represent a single conflict along a node path and arc label specified path
 *
 * path : the node path given by the node indices in each layer
 * label : the arc labels along the path
 * conflict : the two vertices on the path with label 1 that are adjacent and thus form a conflict with that path
 *
 * PathLabelConflict(path, label, conflict) : creates a PathLabelConflict object with given data fields
 */

typedef std::vector<NodeIndex> Path;
enum ArcType{ zeroArc, oneArc };
typedef std::vector<ArcType> Label;
typedef std::tuple<int, int> Conflict;

struct PathLabelConflict{
    Path path;
    Label label;
    Conflict conflict;

    PathLabelConflict(Path path, Label label, Conflict conflict);
};


/*
 * Iterative refinement functions
 * Purpose : functions implementing the main subroutines as described by van Hoeve in his paper
 *
 * RedirectArcs : option to use no redirection of arcs when separating a conflict or redirects to node with
 *                most similar state info, i.e. smallest intersection of state information
 *
 * intersection_size : helper function for redirecting arcs, returns the intersection size of two sets
 * separate_edge_conflict : given a decision diagram dd and a conflict with node path and arc labels (plc)
 *                          separate that conflict in dd. Redirect arcs if option is set
 *
 *
 * Model : option whether flow on edge was computed via integer or linear program
 * ConflictResolution : option which conflicts to look for, returning either a single, multiple or the conflicts on
 *                      paths with a certain amount of flow only
 * PathDecomposition : option to change the way the paths in the decomposition are chosen
 *
 * detect_edge_conflict : with computed flow solution given as flow on the edges of the decision diagram,
 *                        decompose the flow into paths and look for a conflict on such a path.
 *                        returns the empty conflict if none is found
 *
 *
 * Formulation : option to slightly change the IP/LP formulation, like making certain variables integral,
 *               adding extra constraints, imposing bounds on the variables or a bound on the objective value
 *
 * compute_flow_solution : builds up the linear or integer program from the given decision diagram and
 *                         computes an optimal solution. the flow values of the solution are stored
 *                         on the edges of the diagram, i.e. the u.one_arc_flow of a node u
 *
 *
 * find_longest_path : uses the fact that the graph of a decision diagram is acyclic. This allows one to compute
 *                     a longest path in polynomial time using a topological ordering of the graph. additionally,
 *                     the ordering in the layer already yields a topological ordering. then applies Dijkstra
 * conflict_on_longest_path : given a longest or any path, returns a minimal conflict on that path.
 *                            if no conflict exists on that path, returns the empty conflict
 */
enum RedirectArcs {OriginalArcs, MostSimilarNode};

int intersection_size(const StateInfo &a, const StateInfo &b);
void separate_edge_conflict(DecisionDiagram &dd, const NeighborList &neighbors, const PathLabelConflict &plc,
                            RedirectArcs redirect_arcs = OriginalArcs);

enum Model {IP, LP};
enum ConflictResolution {SingleConflict, MultipleConflicts, LargestFlowConflict};
enum PathDecomposition {PreferOneArcs, AvoidConflicts, PreferZeroArcs};

std::vector<PathLabelConflict>
detect_edge_conflict(DecisionDiagram dd, const NeighborList &neighbors, double flow_val, Model model = IP,
                     ConflictResolution find_conflicts = SingleConflict,
                     PathDecomposition path_decomposition = PreferOneArcs);

enum Formulation {Normal, BoundConstraints, ExtraConstraints, VarColorBound, ObjectiveValueBound, OneArcsContinuous};

double compute_flow_solution(DecisionDiagram &dd, Model model = IP, int coloring_upper_bound = -1,
                      Formulation formulation = Normal);

void find_longest_path(const DecisionDiagram &dd, Path &path, Label &label);

PathLabelConflict conflict_on_longest_path(const DecisionDiagram &dd, const NeighborList &neighbors);



//TODO remove??
//attempt at implementing the min state compilation ordering,
// always chooses the vertex appearing in the minimal number of parent states
// results were not beneficial
DecisionDiagram exact_decision_diagram_test(const Graph &g, const NeighborList &neighbors);

//attempt at skipping the building step of the dd but not beneficial either, about the same result
double exact_dd_compute_flow_solution(const NeighborList &neighbors, Model model = IP, int coloring_upper_bound = -1,
                      Formulation formulation = Normal);


#endif //DDCOLORS_DECISIONDIAGRAM_H
