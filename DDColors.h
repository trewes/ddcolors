#ifndef DDCOLORS_DDCOLORS_H
#define DDCOLORS_DDCOLORS_H

#include <vector>
#include <set>
#include <chrono>
#include <stdexcept>
#include <numeric>

#include "DecisionDiagram.h"
#include "scip_interface.h"

typedef std::set<unsigned int> ColorClass;
typedef std::vector<ColorClass> Coloring;


struct Statistics{
    size_t decision_diagram_size = 0;
    size_t decision_diagram_arcs = 0;
    size_t decision_diagram_width = 0;
    unsigned int num_longest_path_conflicts = 0;
    unsigned int num_conflicts = 0;
    unsigned int num_lp_solved = 0;
    unsigned int num_ip_solved = 0;
    std::chrono::steady_clock::time_point start_time;
    std::chrono::duration<double> execution_time;
    //time is excluding initialisation of graph but including finding vertex ordering and dsatur bound
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
    enum Algorithm {BasicRefinement, HeuristicRefinement, ExactCompilation};
    Algorithm algorithm = HeuristicRefinement;
    RedirectArcs redirect_arcs = OriginalArcs;
    ConflictResolution find_conflicts = MultipleConflicts;
    Graph::OrderType vertex_ordering = Graph::Dsatur;
    Relaxation relaxation = LP_First;
    int num_longest_path_iterations = 100;
    bool print_stats = true;
    bool print_time = true;
    bool dsatur_only = false;
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

    int heuristic_iterative_refinement();


    int heuristic_iterative_refinement_NUMERICALLY_SAFE();


    //Extra stuff
    //preferably used in LP case as otherwise either found conflicts or primal heuristic yield different results
    //currently the same conflicts are found but the coloring size of the primal heuristic might differ
    std::pair<std::vector<PathLabelConflict>, int>
    find_conflict_and_primal_heuristic(DecisionDiagram &dd, double flow_val, Model model = IP);

private:
    Graph graph;
    NeighborList neighbors;
    const Options opt;
    Statistics stats;
    int dsatur_bound = std::numeric_limits<int>::max();

    void longest_path_refinement(DecisionDiagram &dd);

    void preprocessing_graph();
};



#endif //DDCOLORS_DDCOLORS_H
