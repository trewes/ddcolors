#ifndef DDCOLORS_DDCOLORS_H
#define DDCOLORS_DDCOLORS_H

/*
 * DDColors.h
 * Purpose : implementation of a general interface allowing to specify how exactly the chromatic number is going to be
 * calculated and handling the initialisation as well as applying the settings to the algorithm.
 * Auxiliary structs include Statistics for general runtime information of the algorithm and Options to
 * set certain parameters of the algorithm.
 */

#include <vector>
#include <set>
#include <unordered_map>
#include <chrono>
#include <stdexcept>
#include <numeric>

#include "Graph.h"
#include "DecisionDiagram.h"

extern "C" {
#include "color.h"
}

typedef unsigned int Color;
typedef std::set<Color> ColorClass;
typedef std::vector<ColorClass> Coloring;

std::ostream &operator<<(std::ostream &s, std::chrono::duration<double> duration);

/*
 * Statistics
 * Purpose: Used as a field in DDColors to store various statistics accumulated during the execution of the algorithm
 *
 * the member variables and their purpose/represented data are obvious from their name
 * times are always computed excluding the reading in of the graph but including any possible preprocessing
 *
 * print: outputs and describes most fields of the struct in a simple way
 * pretty_time: prints the given time in a human readable format like Xh Ym Zs ...
 *
 */

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
};

/*
 * Options
 * Purpose: Used as a field in DDColors and employed to set specifics of how the algorithm is executed
 *
 * algorithm : switches the method to compute the chromatic number (or other things)
 *             HeuristicRefinement default, iterative refinement with relaxed decision diagrams
 *             ExactCompilation builds exact decision diagram and solves IP to get chromatic number
 *             HeuristicOnly only performs the coloring heuristics to return an upper bound
 *             ExactFractionalNumber builds exact decision diagram and solves LP to get the fractional! chromatic number
 * redirect_arcs : option to use no redirection of arcs when separating a conflict or redirects to node with
 *                 most similar state info, i.e. smallest intersection of state information
 * find_conflicts : option which conflicts to look for, returning either a single, multiple or the conflicts on
 *                  paths with a certain amount of flow only
 * largest_conflicts_limit : specifies how large the flow must be to consider the conflict on that path,
 *                           number defines 1e^-X and flow has to be larger than that
 * vertex_ordering : specifies which vertex ordering obtained by what coloring heuristic/ordering is used
 * relaxation : allows to change that during the iterative refinement, the IP is always solved and not only after
 *              LP finds no conflicts. if set to mixing, every kth LP is solved as IP instead
 * ip_frequency : sets the frequency of how often the IP is solved instead of the LP in previous setting
 * num_longest_path_iterations : sets the number of iterations we use the longest path computation to look fo a conflict
 * preprocess_graph : enables or disables preprocessing the graph via peeling and removing dominated vertices
 * preprocessing_hint : can set a hint lower bound to use as the peeling degree, can be lager than a simple clique
 * path_decomposition : specifies how the path decomposition should operate, it can either always prefer to chose the
 *                      1- or 0-arc if it has flow or chooses 1-arc unless the 0-arc can be chosen to avoid a conflict
 * formulation : option to slightly change the IP/LP formulation, like making certain variables integral,
 *               adding extra constraints, imposing bounds on the variables or a bound on the objective value
 * ordering_random_tiebreaks : enables or disables the use of randomness in choices during the coloring heuristics
 * use_clique_in_ordering : a found clique can be passed to the orderings to fix and color those vertices first,
 *                          often this can lead to better upper bounds, might transfer to a better variable ordering too
 * clique_num_branches : specifies the number of branches to explore before returning some clique
 *                       (if number is large enough, the algorithm will find a largest clique)
 */

struct Options{
    enum Algorithm {BasicRefinement, HeuristicRefinement, ExactCompilation, HeuristicOnly, ExactFractionalNumber};
    Algorithm algorithm = HeuristicRefinement;
    RedirectArcs redirect_arcs = OriginalArcs;
    ConflictResolution find_conflicts = MultipleConflicts;
    int largest_conflicts_limit = 5;
    Graph::OrderType vertex_ordering = Graph::MaxConnectedDegree;
    enum Relaxation {IP_Only, LP_First, Mixed};
    Relaxation relaxation = LP_First;
    int ip_frequency = 10;
    int num_longest_path_iterations = 100;
    bool preprocess_graph = false;
    int preprocessing_hint = 0;
    PathDecomposition path_decomposition = PreferOneArcs;//AvoidConflicts;
    Formulation formulation = Normal;
    bool ordering_random_tiebreaks = false;
    bool use_clique_in_ordering = false;
    int clique_num_branches = -1;
    int verbosity_frequency = 0;
    int size_limit = 1;//in million nodes
#ifdef EXTENDED_EXACTCOLORS
    int MIP_emphasis = 0;
    int num_cores = 1;
#endif
};

/*
 * DDColors
 * Purpose: class that manages the high-level execution of the algorithm to calculate the chromatic number
 *
 */

class DDColors{
public:

    /*
     * Member variables
     * graph : the underlying graph which to find the chromatic number for
     * neighbors : the adjacency data of the vertices of the graph in a more useful format than just the edge list
     * opt : an options object that gets passed during the initialisation and sets how and what algorithm is executed
     * stats : a stats object in which all the execution information of the algorithm is stored
     * lower_bound : stores the lower bound if one is obtained during any algorithm
     * upper_bound : stores the lower bound if one is obtained during any algorithm
     * heuristic_bound : stores the upper bound obtained by the heuristic coloring algorithms
     * coloring_bound : stores the upper bound obtained by applying the primal heuristic to a flow solution on a diagram
     * fractional_chromatic : stores the fractional chromatic number if options are set to compute it
     * clique : stores a clique if option is set to use one for heuristic coloring and variable ordering
     *
     *
     * DDColors : builds a DDColors object with the underlying graph either read in from a given file
     *            or    passed as an argument
     * ~DDColors : defulat destructor, only takes care of releasing the environment used by the MIP solver
     * run : begins execution of the algorithm given by the settings object
     * get_upper_bound : returns the member field upper_bound
     * get_fractional_chromatic_number : returns the member field fractional_chromatic
     */

    DDColors(const char *filename, Options options = Options());
    DDColors(Graph in_graph, Options options = Options());
    ~DDColors();

    int run();

    int get_upper_bound() const {return upper_bound;}

    double get_fractional_chromatic_number() const {return fractional_chromatic;}

private:
    Graph graph;
    NeighborList neighbors;
    const Options opt;
    Statistics stats;
    int lower_bound = 0;
    int upper_bound = std::numeric_limits<int>::max();
    int heuristic_bound = std::numeric_limits<int>::max();
    int coloring_bound = std::numeric_limits<int>::max();
    double fractional_chromatic = 0;
    std::vector<Vertex> clique = {};

    /*
     * Several methods implementing the various algorithms and subroutines thereof to compute the chromatic number
     * either through the iterative framework on relaxed decision diagrams or on exact decision diagrams
     *
     * initialise : after graph is read in, does the first preprocessing steps like calling the preprocessing on the
     *              graph and computing the heuristic bound as well as the permuted graph to the set variable ordering
     * longest_path_refinement : calls the longest path routine before beginning iterative refinement, for as many times
     *                           as given in settings, looks for and separates an edge conflict on a longest path
     * basic_iterative_refinement : implements the iterative refinement algorithm without the primal heuristic
     * primal_heuristic : the primal heuristic using the computed flow on the edges of the decision diagram
     * find_conflict_and_primal_heuristic : combines finding an edge conflict and the primal heuristic into
     *                                      a single function, useful since both use a similar path decomposition and
     *                                      both invalidate the flow on edges, this way we don't need to work on a copy
     * heuristic_try_color_swap : implements the recolor technique for the primal heuristics. if a new color class
     *                            would be created, check if we can swap the colors of two vertices to avoid
     *                            having to use a new color class. This can improve the heuristic
     * heuristic_iterative_refinement : implements the iterative refinement algorithm with the primal heuristic
     *                                  to allow to terminate when lower_bound = upper_bound
     * preprocessing_graph : preprocesses the underlying graph via peeling and removing dominated vertices,
     *                       uses COLORclique_ostergard to get a clique size as lower bound and use as peeling degree
     * update_lower_bound, update_upper_bound, print_bounds : self explanatory by their name
     */

    void initialise();

    void longest_path_refinement(DecisionDiagram &dd);

    int basic_iterative_refinement();

    //pass dd here as a copy so we keep flow on original dd
    Coloring primal_heuristic(DecisionDiagram dd);

    static bool heuristic_try_color_swap(Vertex vertex, const NeighborList &neighbors, Coloring &coloring);

    int heuristic_iterative_refinement();

    //preferably used in LP case as otherwise either found conflicts or primal heuristic yield different results
    //currently the same conflicts are found but the coloring size of the primal heuristic might differ
    std::pair<std::vector<PathLabelConflict>, int>
    find_conflict_and_primal_heuristic(DecisionDiagram &dd, double flow_val, Model model = IP);

    void preprocessing_graph(int nbranches = -1);

    void update_lower_bound(int flow_bound);

    void update_upper_bound(int coloring_size);

    void print_bounds(const DecisionDiagram &dd, double obj_value) const;

};

#endif //DDCOLORS_DDCOLORS_H
