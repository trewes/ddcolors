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
#include <ilcplex/cplex.h> //needed for CPX_INFBOUND
}

//catch resp. avoid floating point errors occurring during the path decomposition
#define double_eps std::pow(10, -8)

template<typename T>
std::ostream &operator<<(std::ostream &s, const std::set<T> &set);

template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec);

typedef int LayerIndex;
typedef int NodeIndex;
typedef int EdgeIndex;
typedef unsigned int Vertex;
typedef std::set<Vertex> StateInfo;

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

std::vector<PathLabelConflict>
detect_edge_conflict(DecisionDiagram dd, const NeighborList &neighbors, double flow_val, Model model = IP,
                     ConflictResolution find_conflicts = SingleConflict,
                     PathDecomposition path_decomposition = PreferOneArcs);

enum Formulation {Normal, BoundConstraints, ExtraConstraints, VarColorBound, ObjectiveValueBound, OneArcsContinuous};

double compute_flow_solution(DecisionDiagram &dd, Model model = IP, int coloring_upper_bound = -1,
                      Formulation formulation = Normal);

void find_longest_path(const DecisionDiagram &dd, Path &path, Label &label);

PathLabelConflict conflict_on_longest_path(const DecisionDiagram &dd, const NeighborList &neighbors);

std::vector<PathLabelConflict>
experimental_conflict_on_longest_path(const DecisionDiagram &dd, const NeighborList &neighbors);


//attempt at implementing the min state compilation ordering,
// always chooses the vertex appearing in the minimal number of parent states
// resulte were not beneficial
DecisionDiagram exact_decision_diagram_test(const Graph &g, const NeighborList &neighbors);

//attempt at skipping the building step of the dd but not beneficial either, about the same result
double exact_dd_compute_flow_solution(const NeighborList &neighbors, Model model = IP, int coloring_upper_bound = -1,
                      Formulation formulation = Normal);


#endif //DDCOLORS_DECISIONDIAGRAM_H
