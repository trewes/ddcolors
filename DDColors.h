#ifndef DDCOLORS_DDCOLORS_H
#define DDCOLORS_DDCOLORS_H

#include <vector>
#include <set>
#include <unordered_map>
#include <chrono>
#include <stdexcept>
#include <numeric>

#include "Graph.h"
#include "DecisionDiagram.h"
#include "scip_interface.h"

extern "C" {
#include "color.h"
}

typedef unsigned int Color;
typedef std::set<Color> ColorClass;
typedef std::vector<ColorClass> Coloring;

std::ostream &operator<<(std::ostream &s, std::chrono::duration<double> duration);

struct Statistics{
    size_t decision_diagram_size = 0;
    size_t decision_diagram_arcs = 0;
    size_t decision_diagram_width = 0;
    unsigned int num_longest_path_conflicts = 0;
    unsigned int num_conflicts = 0;
    unsigned int num_lp_solved = 0;
    unsigned int num_ip_solved = 0;
    std::chrono::steady_clock::time_point start_time;
    //time is excluding initialisation of graph but including finding vertex ordering and dsatur bound
    std::chrono::duration<double> execution_time{};
    std::chrono::duration<double> tt_lower_bound{};
    std::chrono::duration<double> tt_upper_bound{};

    void print() const;

    static void pretty_time(std::chrono::duration<double> duration);
//    std::chrono::duration<double> lp_time;
//    std::chrono::duration<double> ip_time;
//    std::chrono::duration<double> primal_heuristic;
//    std::chrono::duration<double> conflict_detection;
//    std::chrono::duration<double> conflict_separation;
};


struct Options{
    enum Algorithm {BasicRefinement, HeuristicRefinement, ExactCompilation, HeuristicOnly, ExactFractionalNumber};
    Algorithm algorithm = HeuristicRefinement;
    RedirectArcs redirect_arcs = OriginalArcs;
    ConflictResolution find_conflicts = MultipleConflicts;
    int largest_conflicts_limit = 5;
    Graph::OrderType vertex_ordering = Graph::Max_Connected_degree;
    enum Relaxation {IP_Only, LP_First, Mixed};
    Relaxation relaxation = LP_First;
    int ip_frequency = 10;
    int num_longest_path_iterations = 100;
    bool print_stats = true;
    bool print_time = true;
    bool preprocess_graph = false;
    int preprocessing_hint = 0;
    PathDecomposition path_decomposition = PreferOneArcs;//AvoidConflicts;
    Formulation formulation = Normal;
    bool ordering_random_tiebreaks = false;
    bool use_clique_in_ordering = false;
    int clique_num_branches = -1;
};

class DDColors{
public:
    DDColors(const char *filename, Options options = Options());

    DDColors(Graph in_graph, Options options = Options());

    ~DDColors();

    int run();

    int basic_iterative_refinement();

    //pass dd here as a copy so we keep flow on original dd
    Coloring primal_heuristic(DecisionDiagram dd);

    static Coloring primal_heuristic(DecisionDiagram dd, const NeighborList &neighbors);

    int heuristic_iterative_refinement();


    //Extra stuff
    //preferably used in LP case as otherwise either found conflicts or primal heuristic yield different results
    //currently the same conflicts are found but the coloring size of the primal heuristic might differ
    std::pair<std::vector<PathLabelConflict>, int>
    find_conflict_and_primal_heuristic(DecisionDiagram &dd, double flow_val, Model model = IP);

    void preprocessing_graph(int nbranches = -1);


private:
    //helper function to call when initilaising a new DDColors object
    void initialise();

    //handle updating of information when a new better bound is found
    void update_lower_bound(int flow_bound);

    void update_upper_bound(int coloring_size);

    void print_bounds(const DecisionDiagram& dd) const;

    void longest_path_refinement(DecisionDiagram &dd);

    std::set<Vertex> clique = {};


    Graph graph;
    NeighborList neighbors;
    const Options opt;
    Statistics stats;

    int lower_bound = 0;
    int heuristic_bound = std::numeric_limits<int>::max();
    int coloring_bound = std::numeric_limits<int>::max();
public:
    int upper_bound = std::numeric_limits<int>::max();
    double fractional_chromatic = 0;

};

//helper function to try and swap colors during last phase of the primal heuristic
bool heuristic_try_color_swap(Vertex vertex, const NeighborList &neighbors, Coloring &coloring);


#endif //DDCOLORS_DDCOLORS_H
