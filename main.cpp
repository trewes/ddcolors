#include <iostream>
#include <getopt.h>

#include "DecisionDiagram.h"
#include "DDColors.h"

#include <unistd.h>

//#include <cstdio>
//#include <gurobi_c.h>

#include "scip_interface.h"



void print_help(){
    std::cout<<"Usage: ddcolors.exe [options] graph"<<std::endl;
    std::cout<<"Options:"<<std::endl;
    std::cout<<"-h|--help               :Prints this help message."<<std::endl;
    std::cout<<"-i|--info               :Enables output of information gathered during the algorithm."<<std::endl;
    std::cout<<"-t|--time               :Enables output of execution time."<<std::endl;
    std::cout<<"-m|--method arg         :Which algorithm is to be run to get the chromatic number: i for iterative refinement (default),"<<std::endl;//b for basic as hidden option
    std::cout<<"                         e for exact compilation and h for heuristic only (note that this gives only an upper bound)."<<std::endl;
    std::cout<<"-o|--ordering arg       :Option of which vertex ordering is to be used, standard is Dsatur."<<std::endl;
    std::cout<<"                         l for lexicographic, m for max connected degree, d for dsatur ordering (default)."<<std::endl;
    std::cout<<"-r|--relaxation arg     :IP solver relaxation, standard is LP first."<<std::endl;
    std::cout<<"                         i for IP only, l for LP first and then IP (default), s for switching between IP and LP."<<std::endl;
    std::cout<<"-l|--longest_paths arg  :The number of iterations the longest path preprocessing should be applied, standard is 100."<<std::endl;
    std::cout<<"-p|--preprocessing      :Enables the use of peeling and removing dominated vertices of the graph in a preprocessing step."<<std::endl;
    std::cout<<"-a|--arc_redirection    :Redirects arcs to the most similar node in the next level."<<std::endl;
    std::cout<<"-c|--conflicts arg      :Detect and resolve conflicts in a flow solution after chosen method."<<std::endl;
    std::cout<<"                         s for single conflict, m for multiple conflicts (default), f for conflicts with largest flow."<<std::endl;
    std::cout<<"-d|--decomposition arg  :Specifies the heuristic used in the path decomposition of the flow."<<std::endl;
    std::cout<<"                         o is for PreferOneArcs (default), a for paths avoiding conflicts and zero for PreferZeroArcs."<<std::endl;
    std::cout<<"-u|--use_bound          :Makes the solver use the upper bound as bound on the variables in the IP/LP to solve."<<std::endl;
    std::cout<<"-s|--safe_bound         :Uses technique described be Neumaier to get safe lower bounds from the LP solutions."<<std::endl;
    std::cout<<"-f|--formulation arg    :Allows slight changes to the way the IP/LP is formulated. n id for normal (default)"<<std::endl;
    std::cout<<"                         e is for constraint of setting first 1-arc to 1, b for add the variable bounds as constraints, -a for both."<<std::endl;
}

#define MAX_PNAME_LEN 128
static const char* color_verion_string = "some version string";
int program_header(int ac, char **av) {
    int   rval     = 0;
    time_t starttime;
    char my_hostname[MAX_PNAME_LEN];
    pid_t my_pid   = getpid();
    int   i;

    rval = gethostname (my_hostname, MAX_PNAME_LEN - 1);
    if(rval){
        throw std::runtime_error("gethostname failed");
    }

    printf("##############################################################\n");

    printf("Running  :");
    for (i = 0; i < ac; ++i) {printf(" %s", av[i]);}
    printf ("\n");

    printf("Machine  : %s, pid: %lld\n",
           my_hostname, (long long) my_pid);

    (void) time(&starttime);

    printf("Date     : %s", ctime(&starttime));

    printf("Changeset: %s\n", color_verion_string);
    fflush (stdout);
    printf("##############################################################\n");
}

//std::map<char , Options> init_settings(){
//    std::map<char , Options> settings;
//
//    settings['A'].vertex_ordering = Graph::Lexicographic;
//    settings['B'];//already the standard setting
//    settings['C'].vertex_ordering = Graph::Max_Connected_degree;
//    settings['D'].vertex_ordering = Graph::Max_Connected_degree;
//    settings['D'].redirect_arcs = MostSimilarNode;
//    settings['E'].vertex_ordering = Graph::Lexicographic;
//    settings['E'].algorithm = Options::ExactCompilation;
//    settings['F'].algorithm = Options::ExactCompilation;
//    settings['G'].vertex_ordering = Graph::Max_Connected_degree;
//    settings['G'].algorithm = Options::ExactCompilation;
//    settings['H'].dsatur_only = true;//ignores all other options
//    return settings;
//}


int main(int argc, char *argv[]) {

    if(0){
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        NeighborList nlist = g.get_neighbor_list();
        for(const auto &node : nlist){
            std::cout << node << std::endl;
        }

        DecisionDiagram dd = exact_decision_diagram(g, nlist);
        print_decision_diagram(dd, false);

        dd.clear();
        dd = initial_decision_diagram(g);
        std::cout << "initial diagram:" << std::endl;
        print_decision_diagram(dd, false);
        std::cout << "after separating conflict (1,3):" << std::endl;
        separate_edge_conflict(dd, nlist, PathLabelConflict({0, 0, 0}, {oneArc, oneArc, oneArc}, std::make_tuple(1, 3)),
                               OriginalArcs);
        print_decision_diagram(dd, false);
        int flow_val = compute_flow_solution(dd, IP);
        for(auto &plc : detect_edge_conflict(dd, nlist, flow_val, IP, MultipleConflicts, PreferOneArcs)){
            int j, k;
            std::tie(j, k) = plc.conflict;
            std::cout << "found conflict " << j << ", " << k << std::endl;
            for(int u : plc.path){
                std::cout << u << " ";
            }
            std::cout << std::endl;
            for(bool u : plc.label){
                std::cout << u << " ";
            }
            std::cout << std::endl;
            separate_edge_conflict(dd, nlist, plc, OriginalArcs);
        }
//        std::cout << "after separating a second time with conflict (3,4):" << std::endl;
//        separate_edge_conflict(dd, nlist, , OriginalArcs);
//        print_decision_diagram(dd);
    }

    if(0){
        std::cout << "Hello" << std::endl;
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        NeighborList nlist = g.get_neighbor_list();
        DecisionDiagram dd = initial_decision_diagram(g);
        separate_edge_conflict(dd, nlist, PathLabelConflict({0, 0, 0}, {oneArc, oneArc, oneArc}, std::make_tuple(1, 3)),
                               OriginalArcs);
        print_decision_diagram(dd, false);
        int flow_val = compute_flow_solution(dd, IP);

        std::vector<PathLabelConflict> conflict_info = detect_edge_conflict(dd, nlist, flow_val, IP, MultipleConflicts,
                                                                            PreferOneArcs);
        if(conflict_info.empty()){
            std::cout << "Failed to detect any edge conflict with given decision diagram and flow input." << std::endl;
        } else{
            for(const PathLabelConflict &plc : conflict_info){
                int j, k;
                std::tie(j, k) = plc.conflict;
                std::cout << "found conflict " << j << ", " << k << std::endl;
                for(int u : plc.path){
                    std::cout << u << " ";
                }
                std::cout << std::endl;
                for(bool u : plc.label){
                    std::cout << u << " ";
                }
                std::cout << std::endl;
            }
        }
        std::cout << "Done with computing the flow" << std::endl;
    }

    if(0){
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        NeighborList nlist = g.get_neighbor_list();

        for(const auto &node : nlist){
            for(int neighbor : node){
                std::cout << neighbor << " ";
            }
            std::cout << std::endl;
        }

        DecisionDiagram dd = initial_decision_diagram(g);

        print_decision_diagram(dd, false);

        std::cout << "Try to detect edge conflict." << std::endl;

        for(auto &layer : dd){
            layer[0].one_arc_flow = 1;
        }

        double flow_value = 1.0;
        std::vector<PathLabelConflict> conflict_info = detect_edge_conflict(dd, nlist, flow_value, IP,
                                                                            MultipleConflicts, PreferOneArcs);
        if(conflict_info.empty()){
            std::cout << "Failed to detect any edge conflict with given decision diagram and flow input." << std::endl;
        } else{
            for(const PathLabelConflict &plc : conflict_info){
                int j, k;
                std::tie(j, k) = plc.conflict;
                std::cout << "found conflict " << j << ", " << k << std::endl;
                for(int u : plc.path){
                    std::cout << u << " ";
                }
                std::cout << std::endl;
                for(bool u : plc.label){
                    std::cout << u << " ";
                }
                std::cout << std::endl;

                separate_edge_conflict(dd, nlist, plc, OriginalArcs);

            }

        }

    }

    if(0){
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        DDColors ddcol(g, {});
        int lower_bound = ddcol.basic_iterative_refinement();

        std::cout << "lower bound is " << lower_bound << std::endl;
    }

    if(0){
        std::string filename = "../Graphs/mulsol.i.1.col";
        DDColors ddcol(filename.c_str(), {});
        //int lower_bound = ddcol.basic_iterative_refinement();
        int lower_bound = ddcol.heuristic_iterative_refinement();
        std::cout << "lower bound is " << lower_bound << std::endl;
    }

    if(0){
        Graph g("../Graphs/myciel4.col");
        NeighborList nlist = g.get_neighbor_list();
        DecisionDiagram dd = exact_decision_diagram(g, nlist);
        print_decision_diagram(dd, false);
        int lower_bound = compute_flow_solution(dd);

        std::cout << "lower bound is " << lower_bound << std::endl;
    }

    if(0){
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        Graph perm_graph = g.perm_graph({2, 1, 3, 4});
        for(auto n : perm_graph.get_neighbor_list()){
            std::cout << n << std::endl;
        }
    }

    if(0){
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        DDColors ddcol(g, {});
        int lower_bound = ddcol.run();

        std::cout << "lower bound is " << lower_bound << std::endl;
    }

    if(0){
        Options opt{};
        opt.relaxation = Switch_LP_IP; //opt.algorithm = Options::HeuristicRefinement;
        //opt.redirect_arcs = RedirectArcs::MostSimilarNode;
        //opt.vertex_ordering = Graph::OrderType::Lexicographic;
        opt.num_longest_path_iterations = 1;
        std::string filename = "../Graphs/myciel3.col";//"../Graphs/DSJC125.1.col";
        DDColors ddcol(filename.c_str(), opt);
        int lower_bound = ddcol.run();

        std::cout << "lower bound is " << lower_bound << std::endl;
    }

    if(0){
        Graph g("../Graphs/mulsol.i.1.col");
        Permutation ordering;
        Coloring c = g.dsatur(ordering);
        std::cout << "DSatur upper bound is " << c.size() << std::endl;
    }

    if(0){
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
//        g = Graph("../Graphs/mulsol.i.1.col");
        DecisionDiagram dd = initial_decision_diagram(g);
        separate_edge_conflict(dd, g.get_neighbor_list(),
                               PathLabelConflict({0, 0, 0}, {oneArc, oneArc, oneArc}, std::make_tuple(1, 3)));
        PathLabelConflict plc = conflict_on_longest_path(dd, g.get_neighbor_list());
        separate_edge_conflict(dd, g.get_neighbor_list(), plc);
        plc = conflict_on_longest_path(dd, g.get_neighbor_list());
        separate_edge_conflict(dd, g.get_neighbor_list(), plc);
        plc = conflict_on_longest_path(dd, g.get_neighbor_list());
        separate_edge_conflict(dd, g.get_neighbor_list(), plc);
        plc = conflict_on_longest_path(dd, g.get_neighbor_list());

    }

    if(0){
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        Options opt;
        opt.num_longest_path_iterations = 0;
        opt.vertex_ordering = Graph::Lexicographic;
        opt.find_conflicts = SingleConflict;
        opt.algorithm = Options::BasicRefinement;
        opt.relaxation = LP_First;
//        opt.algorithm = Options::ExactCompilation;
        DDColors ddcolors("../Graphs/myciel4.col", opt);
//        DDColors ddcolors(g, opt);
        std::cout << "chromatic number is: " << ddcolors.run() << std::endl;
    }

    if(0){
        std::vector<int> edge_list;
        int cycle_length = 7;
        for(int i = 1; i < cycle_length; i++){
            edge_list.push_back(i);
            edge_list.push_back(i + 1);
        }
        edge_list.push_back(cycle_length);
        edge_list.push_back(1);
//        Graph g(7, 7, {1, 3, 3, 5, 5, 7, 7, 2, 2, 4, 4, 6, 6, 1});
//        Graph g(cycle_length, cycle_length, edge_list);
        Graph g("../Graphs/myciel4.col");
        Options opt;
        opt.algorithm = Options::ExactCompilation;
//        opt.vertex_ordering = Graph::Max_Conneceted_degree;
        opt.num_longest_path_iterations = 0;
//        opt.relaxation = IP_Only;
//        opt.find_conflicts = LargestFlowConflict;

        DDColors(g, opt).run();
    }

    if(0){ //testing integrality gap for all small graphs
        Options opt;
        opt.algorithm = Options::ExactCompilation;
        opt.vertex_ordering = Graph::Lexicographic;
        opt.print_stats = false;
        opt.print_time = false;

        std::ifstream file("../Graphs/graph6/graph9c.g6");
        if(not file){
            throw std::runtime_error("Cannot open file.");
        }

        std::string line;
        std::getline(file, line);//read in the file line by line
        int graph_count = 0;
        do{
            graph_count++;
            std::cout << "graph " << graph_count << std::endl;
            Graph g(line);
            DDColors(g, opt).run();
        } while(std::getline(file, line));
    }

    if(0){
//        Graph g(4, 3, {1, 3, 3, 4, 4, 2});

        Graph g("../Graphs/myciel4.col");
        Options opt;
//        opt.algorithm = Options::ExactCompilation;
        opt.vertex_ordering = Graph::Lexicographic;
        opt.num_longest_path_iterations = 0;

        DDColors(g, opt).run();
    }

    if(0){
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
//        Graph g("../Graphs/inithx.i.3.col");
        Options opt;
        opt.vertex_ordering = Graph::Lexicographic;
        opt.num_longest_path_iterations = 0;
//        opt.algorithm = Options::ExactCompilation;
//        opt.relaxation = IP_Only;
//        opt.find_conflicts = LargestFlowConflict;
//        opt.dsatur_only = true;

        std::cout << DDColors(g, opt).run() << std::endl;
    }


//    if(1){
//        struct COLORlp {
//            GRBmodel *model;
//            double dbl_cutoff;
//        };
//
//        static     GRBenv *grb_env = NULL;
//        //init env
//        int rval = 0;
//        if (!grb_env) {
//            rval = GRBloadenv (&grb_env, NULL);
//            std::cout << "code: " << rval << std::endl;
//        }
//        COLORlp* p;
//    }


    if(0){
        SCIP_RETCODE retcode;

        retcode = runCircle();

        /* evaluate return code of the SCIP process */
        if(retcode != SCIP_OKAY){
            /* write error back trace */
            SCIPprintError(retcode);
            return -1;
        }
    }


    if(0){
        test();
    }


    if(0){
//        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        Graph g("../Graphs/mulsol.i.1.col");
        Permutation perm;
        int bound = g.max_connected_degree_coloring(perm).size();
        g = g.perm_graph(perm);
        DecisionDiagram dd = exact_decision_diagram(g);
//        compute_flow_solution(dd, IP, -1);
        compute_flow_solution_SCIP(dd, IP, -1);
        validate_flow_solution_SCIPexact(dd, 0, IP, -1);
    }


    if(0){
        std::vector<Vertex> edge_list;
        int length = 7;
        for(int i = 1; i <= 3; i++){
            for(int j = i + 1; j <= length; j++){
                edge_list.push_back(i);
                edge_list.push_back(j);
            }
        }

        for(int i = 5; i <= length; i++){
            edge_list.push_back(i);
            edge_list.push_back(i + 1);
        }

        Graph g(length + 1, edge_list.size() / 2, edge_list);
        g.print();
        g.remove_dominated_vertices();
        g.print();
        g.remove_dominated_vertices();
        g.print();



//        h = g.peel_graph_ordered(4);
//        h.print();
//        h = h.peel_graph_ordered(5);
//        h.print();

    }

    if(0){
        Graph g("../Graphs/inithx.i.3.col");
        std::cout << "original: " << g.ncount() << ", " << g.ecount() << std::endl;
        int bound = 30;
        Options opt;
//        opt.algorithm = Options::ExactCompilation;
//        opt.vertex_ordering = Graph::Lexicographic;
        opt.num_longest_path_iterations = 0;

        DDColors(g, opt).preprocessing_graph();

//        DDColors(g, opt).run();
        g.peel_graph_ordered(bound);
        g.peel_graph_ordered(bound);
        g.peel_graph_ordered(bound);
//        DDColors(g, opt).run();
        g.peel_graph_ordered(bound);
//        DDColors(g, opt).run();
        g.peel_graph_ordered(bound);
        DDColors(g, opt).run();

    }

    if(0){
        Graph g("../Graphs/3-Insertions_5.col");

        std::cout << "original: " << g.ncount() << ", " << g.ecount() << std::endl;
        int bound = 4;

        unsigned int num_nodes_previous_iteration = 0;
        while(true){
            num_nodes_previous_iteration = g.ncount();
            g.remove_dominated_vertices();
            g.peel_graph_ordered(bound);
            if(num_nodes_previous_iteration == g.ncount())
                break;
        }

        Options opt;
//        opt.algorithm = Options::ExactCompilation;
//        DDColors(g, opt).run();
    }

    if(0){
        char const *filename = (0) ? "../Graphs/queen9_9.col" : argv[1];
//        Graph g("../Graphs/zeroin.i.3.col");
        Graph g(filename);

        Options opt;
        opt.multiple_dsatur = 1;
        opt.algorithm = Options::HeuristicOnly;

        std::cout << "Heuristic found coloring of size " << DDColors(g, opt).run() << std::endl;
    }

    if(0){
        char const *filename = (1) ? "../Graphs/queen9_9.col" : argv[1];
//        Graph g("../Graphs/zeroin.i.3.col");
        Graph g(filename);

        Permutation perm;
        Coloring c;

        c = g.dsatur(perm);
        g = g.perm_graph(perm_inverse(perm));
        DecisionDiagram dd = exact_decision_diagram(g);
        std::cout << num_nodes(dd) << " nodes, " << num_arcs(dd) << " arcs, " << get_width(dd) << " width. dsatur "
                  << c.size() << std::endl;

        c = g.dsatur_original(perm);
        g = g.perm_graph(perm_inverse(perm));
        dd = exact_decision_diagram(g);
        std::cout << num_nodes(dd) << " nodes, " << num_arcs(dd) << " arcs, " << get_width(dd)
                  << " width. dsatur_original " << c.size() << std::endl;

        c = g.max_connected_degree_coloring(perm);
        g = g.perm_graph(perm_inverse(perm));
        dd = exact_decision_diagram(g);
        std::cout << num_nodes(dd) << " nodes, " << num_arcs(dd) << " arcs, " << get_width(dd)
                  << " width. max_connected_degree_coloring " << c.size() << std::endl;

        g = Graph(filename);
        dd = exact_decision_diagram(g);
        std::cout << num_nodes(dd) << " nodes, " << num_arcs(dd) << " arcs, " << get_width(dd)
                  << " width. lexicographic" << std::endl;
    }

    if(0){
        char const *filename = (0) ? "../Graphs/queen9_9.col" : argv[2];
        Graph g(filename);
        int hint = std::atoi(argv[1]);
        std::cout << hint << std::endl;

        DDColors ddc(g);
        ddc.preprocessing_graph();
    }

    if(0){
        char const *filename = (0) ? "../Graphs/queen9_9.col" : argv[1];
        Graph g(filename);
        NeighborList neighbors = g.get_neighbor_list();
        DecisionDiagram dd = exact_decision_diagram(g);
        double ip_flow = compute_flow_solution(dd, IP, Normal);
        int colorsize = int(DDColors::primal_heuristic(dd, neighbors).size());

        std::cout << "flow bound " << ip_flow << " vs color bound " << colorsize << std::endl;


    }

    if(0){
        char const *filename = (0) ? "../Graphs/r250.1c.col" : argv[1];
        Graph g(filename);
        std::cout << "width is " << g.constraint_graph_width() << std::endl;
    }

//return 0;
    Options dd_settings{};

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
            {"help",            no_argument,       nullptr, 'h'},
            {"info",            no_argument,       nullptr, 'i'},
            {"time",            no_argument,       nullptr, 't'},
            {"arc_redirection", no_argument,       nullptr, 'a'},
            {"preprocessing",   optional_argument, nullptr, 'p'},
            {"method",          required_argument, nullptr, 'm'},
            {"conflicts",       required_argument, nullptr, 'c'},
            {"ordering",        required_argument, nullptr, 'o'},
            {"relaxation",      required_argument, nullptr, 'r'},
            {"longest_paths",   required_argument, nullptr, 'l'},
            {"decomposition",   required_argument, nullptr, 'd'},
            {"use_bound",       no_argument,       nullptr, 'u'},
            {"safe_bound",      no_argument,       nullptr, 's'},
            {"formulation",     no_argument,       nullptr, 'f'},
    };

    while((opt = getopt_long(argc, argv, "hitap::uskc:o:r:l:m:d:f:", long_options, &option_index)) != -1){
        switch(opt){

            default:
            case '?':                                                       //getopt itself will return an error message
                return -1;
            case 'h':
                print_help();
                return -1;

            case 'i':
                dd_settings.print_stats = true;
                break;
            case 't':
                dd_settings.print_time = true;
                break;
            case 'm':
                if(*optarg == 'b'){
                    dd_settings.algorithm = Options::BasicRefinement;
                } else if(*optarg == 'i'){
                    dd_settings.algorithm = Options::HeuristicRefinement;
                } else if(*optarg == 'e'){
                    dd_settings.algorithm = Options::ExactCompilation;
                } else if(*optarg == 'h'){
                    dd_settings.algorithm = Options::HeuristicOnly;
                } else{
                    std::cout << "The method to compute the chromatic number was not correctly specified." << std::endl;
                    std::cout << "Program failed." << std::endl;
                    return -1;
                }
                break;
            case 'a':
                dd_settings.redirect_arcs = RedirectArcs::MostSimilarNode;
                break;
            case 'c':
                if(*optarg == 's'){
                    dd_settings.find_conflicts = SingleConflict;
                } else if(*optarg == 'm'){
                    dd_settings.find_conflicts = MultipleConflicts;
                } else if(*optarg == 'f'){
                    dd_settings.find_conflicts = LargestFlowConflict;
                } else{
                    std::cout << "The conflict detection method was not correctly specified." << std::endl;
                    std::cout << "Program failed." << std::endl;
                    return -1;
                }
                break;
            case 'o':
                if(*optarg == 'l'){
                    dd_settings.vertex_ordering = Graph::OrderType::Lexicographic;
                } else if(*optarg == 'm'){
                    dd_settings.vertex_ordering = Graph::OrderType::Max_Connected_degree;
                } else if(*optarg == 'd'){
                    dd_settings.vertex_ordering = Graph::OrderType::Dsatur;
                } else if(*optarg == 'o'){
                    dd_settings.vertex_ordering = Graph::OrderType::Dsatur_original;
                } else{
                    std::cout << "The ordering type was not correctly specified." << std::endl;
                    std::cout << "Program failed." << std::endl;
                    return -1;
                }
                break;

            case 'r':
                if(*optarg == 'i'){
                    dd_settings.relaxation = IP_Only;
                } else if(*optarg == 'l'){
                    dd_settings.relaxation = LP_First;
                } else if(*optarg == 's'){
                    dd_settings.relaxation = Switch_LP_IP;
                } else{
                    std::cout << "The IP relaxation method was not correctly specified." << std::endl;
                    std::cout << "Program failed." << std::endl;
                    return -1;
                }
                break;
            case 'l':
                dd_settings.num_longest_path_iterations = std::stoi(optarg, nullptr, 10);
                break;
            case 'p':
                dd_settings.preprocess_graph = true;
                if(optarg != nullptr and *optarg != '\0'){
                    dd_settings.preprocessing_hint = std::stoi(optarg, nullptr, 10);
                }
                break;
            case 'd':
                if(*optarg == 'o'){
                    dd_settings.path_decomposition = PreferOneArcs;
                } else if(*optarg == 'a'){
                    dd_settings.path_decomposition = AvoidConflicts;
                } else if(*optarg == 'z'){
                    dd_settings.path_decomposition = PreferZeroArcs;
                } else{
                    std::cout << "The decomposition type was not correctly specified." << std::endl;
                    std::cout << "Program failed." << std::endl;
                    return -1;
                }
                break;
            case 'u':
                dd_settings.use_upperbound_in_IP = true;
                break;
            case 's':
                dd_settings.safe_LP_bounds = true;
                break;
            case 'f':
                if(*optarg == 'n'){
                    dd_settings.formulation = Normal;
                } else if(*optarg == 'b'){
                    dd_settings.formulation = BoundConstraints;
                } else if(*optarg == 'e'){
                    dd_settings.formulation = ExtraConstraints;
                } else if(*optarg == 'a'){
                    dd_settings.formulation = AllConstraints;
                } else{
                    std::cout << "The IP/LP formulation type was not correctly specified." << std::endl;
                    std::cout << "Program failed." << std::endl;
                    return -1;
                }
                break;
            case 'k':
                dd_settings.use_clique_in_ordering = true;//TODO add long command/help
                break;
        }
    }

    //optind is the index in argv after going through all the options, now the arguments are given
    char const *filename = (argc > 1) ? argv[optind] : "../Graphs/r250.1c.col";

    try{
        program_header(argc, argv);
        std::cout << "Begin ddcolors:" << std::endl;
        Graph g(filename);
        DDColors ddcolors(g, dd_settings);

        int chromatic_number = ddcolors.run();

        if(dd_settings.algorithm == Options::HeuristicOnly){
            std::cout << "Heuristic bound of input graph is " << chromatic_number << std::endl;
        } else
            std::cout << "Chromatic number of input graph is " << chromatic_number << std::endl;
    }
    catch(const std::runtime_error &e){
        std::cout << e.what() << std::endl;
        std::cout << "Program failed." << std::endl;
        return -1;
    }
    catch(...){
        std::cout << "Program failed for an unknown reason." << std::endl;
        return -1;
    }
    std::cout << "Finished." << std::endl;
    return 0;
}
