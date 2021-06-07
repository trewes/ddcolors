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
    std::cout<<"-s|--stats              :Enables output of statistics gathered during the algorithm."<<std::endl;
    std::cout<<"-t|--time               :Enables output of execution time."<<std::endl;
    std::cout<<"-b|--basic_refinement   :Uses refinement algorithm without primal heuristic or Dsatur bound."<<std::endl;
    std::cout<<"-e|--exact_compilation  :Uses exact method to compile the decision diagram ."<<std::endl;
    std::cout<<"-a|--arc_redirection    :Redirects arcs to the most similar node in the next level."<<std::endl;
    std::cout<<"-d|--dsatur_only        :Only runs Dsatur to get and return upper bound. Warning: algorithm does not return a lower bound."<<std::endl;
    std::cout<<"-c|--conflicts arg      :Detect and resolve conflicts in a flow solution after chosen method."<<std::endl;
    std::cout<<"                         s for single conflict, m for multiple conflicts (default), f for conflicts with largest flow."<<std::endl;
    std::cout<<"-o|--ordering arg       :Option of which vertex ordering is to be used, standard is Dsatur."<<std::endl;
    std::cout<<"                         l for lexicographic, m for max connected degree, d for dsatur ordering (default)."<<std::endl;
    std::cout<<"-r|--relaxation arg     :IP solver relaxation, standard is LP first."<<std::endl;
    std::cout<<"                         i for IP only, l for LP first and then IP (default), s for switching between IP and LP."<<std::endl;
    std::cout<<"-p|--preprocessing arg  :The number of iterations the longest path preprocessing should be applied, standard is 100."<<std::endl;
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


int main(int argc, char* argv[]) {

    if(0) {
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        NeighborList nlist = g.get_neighbor_list();
        for (const auto& node : nlist) {
            std::cout << node << std::endl;
        }

        DecisionDiagram dd = exact_decision_diagram(g, nlist);
        print_decision_diagram(dd, false);

        dd.clear();
        dd = initial_decision_diagram(g);
        std::cout << "initial diagram:" << std::endl;
        print_decision_diagram(dd, false);
        std::cout << "after separating conflict (1,3):" << std::endl;
        separate_edge_conflict(dd, nlist, PathLabelConflict({0,0,0}, {oneArc,oneArc,oneArc}, std::make_tuple(1,3)), OriginalArcs);
        print_decision_diagram(dd, false);
        int flow_val = compute_flow_solution(dd);
        for (auto &plc : detect_edge_conflict(dd, nlist, flow_val, IP, MultipleConflicts)){
            int j, k;
            std::tie(j, k) = plc.conflict;
            std::cout << "found conflict " << j << ", " << k << std::endl;
            for (int u : plc.path) {
                std::cout << u << " ";
            }
            std::cout << std::endl;
            for (bool u : plc.label) {
                std::cout << u << " ";
            }
            std::cout << std::endl;
            separate_edge_conflict(dd, nlist, plc, OriginalArcs);
        }
//        std::cout << "after separating a second time with conflict (3,4):" << std::endl;
//        separate_edge_conflict(dd, nlist, , OriginalArcs);
//        print_decision_diagram(dd);
    }

    if(0) {
        std::cout << "Hello" << std::endl;
        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        NeighborList nlist = g.get_neighbor_list();
        DecisionDiagram dd = initial_decision_diagram(g);
        separate_edge_conflict(dd, nlist, PathLabelConflict({0,0,0}, {oneArc,oneArc,oneArc}, std::make_tuple(1,3)), OriginalArcs);
        print_decision_diagram(dd, false);
        int flow_val = compute_flow_solution(dd);

        std::vector<PathLabelConflict> conflict_info = detect_edge_conflict(dd, nlist, flow_val, IP, MultipleConflicts);
        if(conflict_info.empty()) {
            std::cout << "Failed to detect any edge conflict with given decision diagram and flow input." << std::endl;
        } else {
            for (const PathLabelConflict &plc : conflict_info) {
                int j, k;
                std::tie(j, k) = plc.conflict;
                std::cout << "found conflict " << j << ", " << k << std::endl;
                for (int u : plc.path) {
                    std::cout << u << " ";
                }
                std::cout << std::endl;
                for (bool u : plc.label) {
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

        for (const auto &node : nlist) {
            for (int neighbor : node) {
                std::cout << neighbor << " ";
            }
            std::cout << std::endl;
        }

        DecisionDiagram dd = initial_decision_diagram(g);

        print_decision_diagram(dd, false);

        std::cout << "Try to detect edge conflict." << std::endl;

        for(auto& layer : dd){
            layer[0].one_arc_flow = 1;
        }

        double flow_value = 1.0;
        std::vector<PathLabelConflict> conflict_info = detect_edge_conflict(dd, nlist, flow_value, IP,
                                                                            MultipleConflicts);
        if(conflict_info.empty()) {
            std::cout << "Failed to detect any edge conflict with given decision diagram and flow input." << std::endl;
        } else {
            for(const PathLabelConflict & plc : conflict_info){
                int j,k;
                std::tie(j,k) = plc.conflict;
                std::cout << "found conflict " << j << ", " << k << std::endl;
                for(int u : plc.path){
                    std::cout << u << " ";
                }std::cout << std::endl;
                for(bool u : plc.label){
                    std::cout << u << " ";
                }std::cout << std::endl;

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
        Graph perm_graph = g.perm_graph({2,1,3,4});
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
        Options opt{}; opt.relaxation = Switch_LP_IP; //opt.algorithm = Options::HeuristicRefinement;
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
        separate_edge_conflict(dd, g.get_neighbor_list(), PathLabelConflict({0,0,0}, {oneArc,oneArc,oneArc}, std::make_tuple(1,3)));
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
        Options opt; opt.num_longest_path_iterations = 0; opt.vertex_ordering = Graph::Lexicographic;
        opt.find_conflicts = SingleConflict; opt.algorithm = Options::BasicRefinement; opt.relaxation = LP_First;
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
            edge_list.push_back(i+1);
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
        opt.print_stats = false; opt.print_time = false;

        std::ifstream file("../Graphs/graph6/graph9c.g6");
        if (not file) {
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
        if( retcode != SCIP_OKAY )
        {
            /* write error back trace */
            SCIPprintError(retcode);
            return -1;
        }
    }

    if(0){
//        Graph g(4, 3, {1, 3, 3, 4, 4, 2});
        Graph h("../Graphs/miles1000.col");
        Permutation perm;
        int bound = h.dsatur(perm).size();
//        Graph g = h.perm_graph(perm);
        DecisionDiagram dd = exact_decision_diagram(h);
        compute_flow_solution(dd, IP, -1);
        validate_flow_solution_SCIP(dd, 0, IP, bound);
    }


    if(0){
        std::vector<Vertex> edge_list;
        int length = 7;
        for(int i = 1; i <= 3; i++){
            for(int j = i+1; j <= length; j++){
                edge_list.push_back(i);
                edge_list.push_back(j);
            }
        }

        for(int i = 5; i <= length; i++){
                edge_list.push_back(i);
                edge_list.push_back(i+1);
        }

        Graph g(length+1, edge_list.size()/2, edge_list);
        g.print();
        Graph h = g.peel_graph(5);
        h.print();
        h = g.peel_graph(5);
        h.print();
        h = h.peel_graph(5);
        h.print();

    }

    if(1){
        Graph g("../Graphs/inithx.i.3.col");
        int bound = 30;
        Options opt;
//        opt.algorithm = Options::ExactCompilation;
//        opt.vertex_ordering = Graph::Lexicographic;
        opt.num_longest_path_iterations = 0;

        DDColors(g, opt).run();
        std::cout << g.ncount() << ", " << g.ecount() << std::endl;
        g = g.peel_graph(bound);
        std::cout << g.ncount() << ", " << g.ecount() << std::endl;
        g = g.peel_graph(bound);
        std::cout << g.ncount() << ", " << g.ecount() << std::endl;
        g = g.peel_graph(bound);
        std::cout << g.ncount() << ", " << g.ecount() << std::endl;
        DDColors(g, opt).run();

    }

return 0;
    Options dd_settings{};

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
            {"help", no_argument, nullptr, 'h'},
            {"stats", no_argument, nullptr, 's'},
            {"time", no_argument, nullptr, 't'},
            {"basic_refinement", no_argument, nullptr, 'b'},
            {"exact_compilation", no_argument, nullptr, 'e'},
            {"arc_redirection", no_argument, nullptr, 'a'},
            {"dsatur_only", no_argument, nullptr, 'd'},
            {"conflicts", required_argument, nullptr, 'c'},
            {"ordering", required_argument, nullptr, 'o'},
            {"relaxation", required_argument, nullptr, 'r'},
            {"preprocessing", required_argument, nullptr, 'p'},
    };

    while ((opt = getopt_long(argc, argv, "hstbeadc:o:r:p:", long_options, &option_index)) != -1){
        switch (opt) {

            default:
            case '?':                                                       //getopt itself will return an error message
                return -1;
            case 'h':
                print_help();
                return -1;

            case 's':
                dd_settings.print_stats = true;
                break;
            case 't':
                dd_settings.print_time = true;
                break;
            case 'b':
                if(dd_settings.algorithm == Options::Algorithm::ExactCompilation) break; //don't override exact compilation
                dd_settings.algorithm = Options::Algorithm::BasicRefinement;
                break;
            case 'e':
                dd_settings.algorithm = Options::Algorithm::ExactCompilation;
                break;
            case 'a':
                dd_settings.redirect_arcs = RedirectArcs::MostSimilarNode;
                break;
            case 'd':
                dd_settings.dsatur_only = true;
                break;
            case 'c':
                if(*optarg=='s'){
                    dd_settings.find_conflicts = SingleConflict;
                }
                else if(*optarg=='m'){
                    dd_settings.find_conflicts = MultipleConflicts;
                }
                else if(*optarg=='f'){
                    dd_settings.find_conflicts = LargestFlowConflict;
                }
                else{
                    std::cout<<"The conflict detection method was not correctly specified."<<std::endl;
                    std::cout<<"Program failed."<<std::endl;
                    return -1;
                }
                break;
            case 'o':
                if(*optarg=='l'){
                    dd_settings.vertex_ordering = Graph::OrderType::Lexicographic;
                }
                else if(*optarg=='m'){
                    dd_settings.vertex_ordering = Graph::OrderType::Max_Connected_degree;
                }
                else if(*optarg=='d'){
                    dd_settings.vertex_ordering = Graph::OrderType::Dsatur;
                }
                else{
                    std::cout<<"The ordering type was not correctly specified."<<std::endl;
                    std::cout<<"Program failed."<<std::endl;
                    return -1;
                }
                break;

            case 'r':
                if(*optarg=='i'){
                    dd_settings.relaxation = IP_Only;
                }
                else if(*optarg=='l'){
                    dd_settings.relaxation = LP_First;
                }
                else if(*optarg=='s'){
                    dd_settings.relaxation = Switch_LP_IP;
                }
                else{
                    std::cout<<"The IP relaxation method was not correctly specified."<<std::endl;
                    std::cout<<"Program failed."<<std::endl;
                    return -1;
                }
                break;
            case 'p':
                dd_settings.num_longest_path_iterations = std::stoi(optarg, nullptr, 10);
                break;
        }
    }

    //optind is the index in argv after going through all the options, now the arguments are given
    char const* filename = (argc>1) ? argv[optind] : "../Graphs/mulsol.i.5.col";

    try{
        program_header(argc, argv);
        std::cout<<"Begin ddcolors:"<<std::endl;
        Graph g(filename);
        DDColors ddcolors(g, dd_settings);

        int chromatic_number = ddcolors.run();

        if(dd_settings.dsatur_only){
            std::cout<<"Dsatur bound of input graph is "<<chromatic_number <<std::endl;
        } else
        std::cout<<"Chromatic number of input graph is "<<chromatic_number <<std::endl;
    }
    catch(const std::runtime_error& e){
        std::cout<<e.what()<<std::endl;
        std::cout<<"Program failed."<<std::endl;
        return -1;
    }
    catch(...){
        std::cout<<"Program failed for an unknown reason."<<std::endl;
        return -1;
    }
    std::cout<<"Finished."<<std::endl;
    return 0;
}
