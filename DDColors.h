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

typedef std::set<unsigned int> ColorClass;
typedef std::vector<ColorClass> Coloring;

std::ostream& operator<<(std::ostream& s, std::chrono::duration<double> duration);

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
    std::chrono::duration<double> execution_time;
    std::chrono::duration<double> tt_lower_bound;
    std::chrono::duration<double> tt_upper_bound;
    void print() const;
    static void pretty_time(std::chrono::duration<double> duration);
    //TODO maybe add stats for how much time was spent where
//    std::chrono::duration<double> lp_time;
//    std::chrono::duration<double> ip_time;
//    std::chrono::duration<double> primal_heuristic;
//    std::chrono::duration<double> conflict_detection;
//    std::chrono::duration<double> conflict_separation;
};



struct Options{
    enum Algorithm {BasicRefinement, HeuristicRefinement, ExactCompilation, HeuristicOnly};
    Algorithm algorithm = HeuristicRefinement;
    RedirectArcs redirect_arcs = OriginalArcs;
    ConflictResolution find_conflicts = MultipleConflicts;
    //vertex ordering also specifies which heuristic is used for the upper bound
    Graph::OrderType vertex_ordering = Graph::Dsatur;
    Relaxation relaxation = LP_First;
    int num_longest_path_iterations = 100;
    bool print_stats = true;
    bool print_time = true;
    bool preprocess_graph = false;
    int preprocessing_hint = 0;
    //TODO mention this in BA
    PathDecomposition path_decomposition = PreferOneArcs;//AvoidConflicts;
    bool use_upperbound_in_IP = false;
    bool safe_LP_bounds  = false;
    Formulation formulation = Normal;
    bool ordering_random_tiebreaks = false;
    bool use_clique_in_ordering = false;//TODO add option of this to program call
    int multiple_dsatur = 1;//TODO use this or not
    //TODO option for safe lp bounds and changing the LP
};

//TODO maybe start timing when calling run() and not in initialisation
class DDColors{
public:
    DDColors(const char *filename, Options options = Options());
    DDColors(Graph in_graph, Options options = Options());
    ~DDColors();

    int run();

    int basic_iterative_refinement();

    //pass dd here as a copy so we keep flow on original dd
    Coloring primal_heuristic(DecisionDiagram dd);
    static Coloring primal_heuristic(DecisionDiagram dd, const NeighborList& neighbors);

    int heuristic_iterative_refinement();


    int heuristic_iterative_refinement_NUMERICALLY_SAFE();


    //Extra stuff
    //preferably used in LP case as otherwise either found conflicts or primal heuristic yield different results
    //currently the same conflicts are found but the coloring size of the primal heuristic might differ
    std::pair<std::vector<PathLabelConflict>, int>
    find_conflict_and_primal_heuristic(DecisionDiagram &dd, double flow_val, Model model = IP);

    void preprocessing_graph();

private:
    //helper function to call when initilaising a new DDColors object
    void initialise();

    //handle updating of information when a new better bound is found
    void update_lower_bound(int flow_bound);
    void update_upper_bound(int coloring_size);
    void print_bounds() const;


    Graph graph;
    NeighborList neighbors;
    const Options opt;
    Statistics stats;

    int lower_bound = 0;
    int upper_bound = std::numeric_limits<int>::max();
    int heuristic_bound = std::numeric_limits<int>::max();
    int coloring_bound = std::numeric_limits<int>::max();

    void longest_path_refinement(DecisionDiagram &dd);

    std::set<Vertex> clique = {};//TODO
};

//helper function to try and swap colors during last phase of the primal heuristic
bool heuristic_try_color_swap(Vertex vertex, const NeighborList &neighbors, Coloring &coloring);


#endif //DDCOLORS_DDCOLORS_H
