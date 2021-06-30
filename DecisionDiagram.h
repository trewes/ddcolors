#ifndef DDCOLORS_DECISIONDIAGRAM_H
#define DDCOLORS_DECISIONDIAGRAM_H

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

extern "C" {
#include "lp.h"
}

//TODO decide on correct epsilon or better think it why it is necessary (-> simply floating point errors...)
#define double_eps std::pow(10, -7)

template<typename T>
std::ostream &operator<<(std::ostream &s, const std::set<T> &set);

template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec);

typedef int LayerIndex;
typedef int NodeIndex;
typedef int EdgeIndex;
typedef std::set<unsigned int> StateInfo;

struct Node{
    Node();

    Node(LayerIndex layer, NodeIndex index, StateInfo state_info = {});

    LayerIndex layer = -1;
    NodeIndex index = -1;
    StateInfo state_info;
    NodeIndex one_arc = -1;  //index of node in next layer
    NodeIndex zero_arc = -1;

    double one_arc_flow = 0;
    double zero_arc_flow = 0;

    bool equivalent(const Node &other) const;
};


typedef std::vector<std::vector<Node> > DecisionDiagram;

unsigned int num_vars(const DecisionDiagram &dd);

unsigned int num_nodes(const DecisionDiagram &dd);

unsigned int num_arcs(const DecisionDiagram &dd);

unsigned int get_width(const DecisionDiagram &dd);

void print_decision_diagram(const DecisionDiagram &dd, bool use_tag = false);

typedef std::vector<NodeIndex> Path;
enum ArcType{
    zeroArc, oneArc
};
typedef std::vector<ArcType> Label;
typedef std::tuple<int, int> Conflict;

struct PathLabelConflict{
    Path path;
    Label label;
    Conflict conflict;

    PathLabelConflict(Path path, Label label, Conflict conflict);

    bool operator==(const PathLabelConflict &other) const {
        return (path == other.path and label == other.label and conflict == other.conflict);
    }
};


DecisionDiagram exact_decision_diagram(const Graph &g);

DecisionDiagram exact_decision_diagram(const Graph &g, const NeighborList &neighbors);

DecisionDiagram initial_decision_diagram(const Graph &g);

enum RedirectArcs {OriginalArcs, MostSimilarNode};

int intersection_size(const StateInfo &a, const StateInfo &b); //helper function for redirecting arcs
void separate_edge_conflict(DecisionDiagram &dd, const NeighborList &neighbors, const PathLabelConflict &plc,
                            RedirectArcs redirect_arcs = OriginalArcs);

enum Model {IP, LP};
enum ConflictResolution {SingleConflict, MultipleConflicts, LargestFlowConflict};
enum PathDecomposition {PreferOneArcs, AvoidConflicts, PreferZeroArcs};
enum Relaxation {IP_Only, LP_First, Switch_LP_IP};

std::vector<PathLabelConflict>
detect_edge_conflict(DecisionDiagram dd, const NeighborList &neighbors, double flow_val, Model model = IP,
                     ConflictResolution find_conflicts = SingleConflict,
                     PathDecomposition path_decomposition = PreferOneArcs);

enum Formulation {Normal, BoundConstraints, ExtraConstraints, AllConstraints};

double compute_flow_solution(DecisionDiagram &dd, Model model = IP, int coloring_upper_bound = -1,
                      Formulation formulation = Normal);

void find_longest_path(const DecisionDiagram &dd, Path &path, Label &label);

PathLabelConflict conflict_on_longest_path(const DecisionDiagram &dd, const NeighborList &neighbors);

std::vector<PathLabelConflict>
experimental_conflict_on_longest_path(const DecisionDiagram &dd, const NeighborList &neighbors);


DecisionDiagram exact_decision_diagram_test(const Graph &g, const NeighborList &neighbors);


//TODO for benchmarks: consider unlimited/more than 100 longest path iterations, seems quite beneficial (sometimes!)

//TODO Idea: make some variables in LP integral, think about a heuristic which ones.
//TODO One idea is: outgoing edges of root, or, all 1-arcs made integral (+ bound of 1) because obj value only depends on those


//TODO and maybe use flow solution on previous dd to warmstart next LP/IP solving
// for that, keep flow on edges that existed before (those may point to different nodes though) and on new edges that come out of new nodes set flow to 0/1/don't define anything
//this might be what's done already since flow is saved in dd too (although decreased through edge detection) and new nodes are initialised with edge with zero flow
//TODO or, very advanced, keep the same LP/IP instance alive for the whole algorithm but change variables and constraints and solve again

//TODO what a heuristic for the best conflict to choose? or, the best flow to look for (giving the conflicts)!

//TODO check and maybe throw in all cases where we do unsigned int -1. do that quite often

//TODO turn decision diagram into a class, will mae it easier to use options object

//TODO spped up preprocessing / peeling and removing vertices routine

//TODO every 100th or so iteration compute the integer flow for hopefully and probably better primal heuristic results

#endif //DDCOLORS_DECISIONDIAGRAM_H
