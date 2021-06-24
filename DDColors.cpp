#include "DDColors.h"


void Statistics::print() const{
    std::cout << "Size of decision diagram at the end: " << decision_diagram_size << " nodes, " << decision_diagram_arcs
    << " arcs and width " << decision_diagram_width << "." << std::endl;
    if(num_conflicts) {//if we never found any conflicts that means we used the exact decision diagram
        std::cout << "Conflicts found: " << num_conflicts;
        if (num_longest_path_conflicts) {
            std::cout << " with and additional " << num_longest_path_conflicts << " conflicts found in preprocessing.";
        }
        std::cout << std::endl;
    }
    std::cout << "Solved " << num_ip_solved << " IP's and " << num_lp_solved << " LP's" << std::endl;
}


void Statistics::pretty_time(std::chrono::duration<double> duration){

    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    if(hours.count()){
        std::cout << hours.count() << "h ";
    }

    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    if(minutes.count()){
        std::cout << minutes.count() << "m ";
    }

    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    duration -= seconds;
    std::cout << seconds.count() << "s ";

    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
    duration -= milliseconds;
    std::cout << milliseconds.count() << "ms";

    if(not milliseconds.count()) {
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        duration -= microseconds;
        std::cout <<" "<< microseconds.count() << "us";
    }

    std::cout<<"."<<std::endl;
}


std::ostream& operator<<(std::ostream& s, std::chrono::duration<double> duration){
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    if(hours.count()){
        s << hours.count() << "h ";
    }

    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    if(minutes.count()){
        s << minutes.count() << "m ";
    }

    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    duration -= seconds;
    s << seconds.count() << "s ";

    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
    duration -= milliseconds;
    s << milliseconds.count() << "ms";

    if(not milliseconds.count()) {
        auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        duration -= microseconds;
        s <<" "<< microseconds.count() << "us";
    }
}

DDColors::DDColors(const char *filename, Options options) : graph(filename), opt(options), stats(Statistics()){
    initialise();
}


DDColors::DDColors(Graph in_graph, Options options) : graph(std::move(in_graph)), opt(options), stats(Statistics()){
    initialise();
}


void DDColors::initialise(){
    stats.start_time = std::chrono::steady_clock::now();

    if(opt.preprocess_graph){
        preprocessing_graph();
    }

    Permutation graph_ordering;
    switch (opt.vertex_ordering) {
        case Graph::Lexicographic:
            //still use some heuristic for upper bound, which one?
            heuristic_bound = int(graph.dsatur(graph_ordering).size());
            graph_ordering = identity(graph.ncount()); //reset ordering
            break;
        case Graph::Dsatur:
            heuristic_bound = int(graph.dsatur(graph_ordering).size());
            break;
        case Graph::Dsatur_original:
            heuristic_bound = int(graph.dsatur_original(graph_ordering).size());
            break;
        case Graph::Max_Connected_degree:
            heuristic_bound = int(graph.max_connected_degree_coloring(graph_ordering).size());
            break;
    }

    if(opt.algorithm == Options::HeuristicOnly){
        return;//only need heuristic bound, that is already computed
    }

    graph = graph.perm_graph(perm_inverse(graph_ordering));
    neighbors = graph.get_neighbor_list();
    COLORlp_init_env();
}


DDColors::~DDColors() {
    //anything else?
    COLORlp_free_env();
}


int DDColors::run() {
    if(opt.algorithm == Options::HeuristicOnly){
        std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
        stats.execution_time = (end_time-stats.start_time);
        if(opt.print_time){
            std::cout << "Execution took: ";
            Statistics::pretty_time(stats.execution_time);
        }
        //careful, heuristic bound is computed according to vertex ordering (if not lexicographic)
        return heuristic_bound; //and prints no stats
    }

    lower_bound = 0;
    if(opt.algorithm == Options::BasicRefinement){
        lower_bound = basic_iterative_refinement();
    } else if(opt.algorithm == Options::HeuristicRefinement){
        lower_bound = heuristic_iterative_refinement();
    } else if(opt.algorithm == Options::ExactCompilation){
        DecisionDiagram dd = exact_decision_diagram(graph, neighbors);
        stats.decision_diagram_size = num_nodes(dd);
        stats.decision_diagram_arcs = num_arcs(dd);
        stats.decision_diagram_width = get_width(dd);
        double ip_flow = compute_flow_solution(dd, IP, heuristic_bound);//TODO does upper bound make it better?
        lower_bound = int(ip_flow);
        stats.num_ip_solved = 1;
    }
    std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
    stats.execution_time = (end_time-stats.start_time);

    if(opt.print_stats){
        stats.print();                                  //print optional statistics about the execution of the algorithm
    }
    if(opt.print_time){
        std::cout << "Execution took: ";
        Statistics::pretty_time(stats.execution_time);
    }
    return lower_bound;
}



int DDColors::basic_iterative_refinement() {

    Model model = (opt.relaxation == IP_Only) ? IP : LP;

    bool found_solution = false;
    DecisionDiagram dd = initial_decision_diagram(graph);
    lower_bound = 0;

    longest_path_refinement(dd);

    while (not found_solution){
        double flow_value = compute_flow_solution(dd, model, -1); // = obj(F)
        int flow_bound = (model == IP) ? int(std::round(flow_value)) : (int)std::ceil(flow_value);
        if(model == IP) stats.num_ip_solved++; else stats.num_lp_solved++;
        if(flow_bound > lower_bound){
            update_lower_bound(flow_bound);
        }

        std::vector<PathLabelConflict> conflict_info =
                detect_edge_conflict(dd, neighbors, flow_value, model, opt.find_conflicts, opt.path_decomposition);
        if(conflict_info.empty()) {
            if(model == LP){
                std::cout << "Switching from LP to IP since no conflict was found" << std::endl;
                model = IP;
                continue;
            } else break;
        }

        stats.num_conflicts += conflict_info.size();
        for(const PathLabelConflict & plc : conflict_info){
            separate_edge_conflict(dd, neighbors, plc, opt.redirect_arcs);
        }
        if(opt.relaxation == Switch_LP_IP and model == IP){
            model = LP;
            std::cout << "Switching back from IP to LP" << std::endl;
        }
    }

    stats.decision_diagram_size = num_nodes(dd);
    stats.decision_diagram_arcs = num_arcs(dd);
    stats.decision_diagram_width = get_width(dd);
    return lower_bound;
}


Coloring DDColors::primal_heuristic(DecisionDiagram dd) {
        int n = int(dd.size() -1);
        Coloring coloring;
        bool color_used = true;
        std::set<unsigned int> coloring_selected_vertices;

        std::set<unsigned int> all_vertices;
        for (int i = 1; i <= n; ++i){
            all_vertices.insert(all_vertices.end(), i);
        }

        while(int(coloring_selected_vertices.size()) < n and color_used) {
            ColorClass cclass;
            color_used = false;
            Path path;
            Label label;
            std::set<int> neighboring_selected_vertices;
            double path_min_flow = n;
            path.push_back(0); //index of root node

            for (int i = 0; i < n; i++) {
                int u = path.back();
                if ((dd[i][u].one_arc != -1) and dd[i][u].one_arc_flow >= dd[i][u].zero_arc_flow) {
                    path_min_flow = std::min(path_min_flow, dd[i][u].one_arc_flow);
                    path.push_back(dd[i][u].one_arc);
                    label.push_back(oneArc);
                    if ((neighboring_selected_vertices.count(i + 1) == 0) and (coloring_selected_vertices.count(i + 1) == 0)) {
                        cclass.insert(cclass.end(), i + 1);
                        coloring_selected_vertices.insert(coloring_selected_vertices.end(), i + 1);
                        neighboring_selected_vertices.insert(neighboring_selected_vertices.end(), i + 1);
                        neighboring_selected_vertices.insert(neighbors[i].begin(), neighbors[i].end());
                        color_used = true;
                    }
                } else {
                    path_min_flow = std::min(path_min_flow, dd[i][u].zero_arc_flow);
                    path.push_back(dd[i][u].zero_arc);
                    label.push_back(zeroArc);
                }
            }
            if (color_used) {
                coloring.push_back(cclass);
                if (path_min_flow > 0) {
                    for (int i = 0; i < n; i++) {
                        if (label[i]) {
                            dd[i][path[i]].one_arc_flow -= path_min_flow;
                        } else {
                            dd[i][path[i]].zero_arc_flow -= path_min_flow;
                        }
                    }
                }
            }
        }
        if(int(coloring_selected_vertices.size()) < n){
            std::set<int> not_selected_vertices;
            std::set_difference(all_vertices.begin(), all_vertices.end(),
                                coloring_selected_vertices.begin(), coloring_selected_vertices.end(),
                                std::inserter(not_selected_vertices, not_selected_vertices.begin()));
            for(auto it = not_selected_vertices.begin(); it != not_selected_vertices.end();
                                        it = not_selected_vertices.erase(it)){ //return iterator following iterator that was deleted, acts as ++
                unsigned int i = *it;
                bool found_non_conflicting_cc;
                for(ColorClass& c : coloring){
                    found_non_conflicting_cc = true;
                    for(unsigned int j : c){
                        if(neighbors[i - 1].count(j)){
                            found_non_conflicting_cc = false;
                            break;
                        }
                    }
                    if(found_non_conflicting_cc){
                        c.insert(c.end(), i);
                        coloring_selected_vertices.insert(coloring_selected_vertices.end(), i);
                        break;
                    }
                }
                if(not coloring_selected_vertices.count(i)) {
                    //could not assign i to any cc, swap colors around or create new class to then assign it a colors
                    bool recolored = heuristic_try_color_swap(i, neighbors, coloring);
                    if(not recolored){
                        coloring.push_back({i});
                    } else {
//                        std::cout << "Swapped color in heuristic" << std::endl;
                    }
                    coloring_selected_vertices.insert(coloring_selected_vertices.end(), i);
                }

            }
        }
    return coloring;
}

Coloring DDColors::primal_heuristic(DecisionDiagram dd, const NeighborList& neighbors) {
    int n = int(dd.size() -1);
    Coloring coloring;
    bool color_used = true;
    std::set<unsigned int> coloring_selected_vertices;

    std::set<unsigned int> all_vertices;
    for (int i = 1; i <= n; ++i){
        all_vertices.insert(all_vertices.end(), i);
    }

    while(int(coloring_selected_vertices.size()) < n and color_used) {
        ColorClass cclass;
        color_used = false;
        Path path;
        Label label;
        std::set<int> neighboring_selected_vertices;
        double path_min_flow = n;
        path.push_back(0); //index of root node

        for (int i = 0; i < n; i++) {
            int u = path.back();
            if ((dd[i][u].one_arc != -1) and dd[i][u].one_arc_flow >= dd[i][u].zero_arc_flow) {
                path_min_flow = std::min(path_min_flow, dd[i][u].one_arc_flow);
                path.push_back(dd[i][u].one_arc);
                label.push_back(oneArc);
                if ((neighboring_selected_vertices.count(i + 1) == 0) and (coloring_selected_vertices.count(i + 1) == 0)) {
                    cclass.insert(cclass.end(), i + 1);
                    coloring_selected_vertices.insert(coloring_selected_vertices.end(), i + 1);
                    neighboring_selected_vertices.insert(neighboring_selected_vertices.end(), i + 1);
                    neighboring_selected_vertices.insert(neighbors[i].begin(), neighbors[i].end());
                    color_used = true;
                }
            } else {
                path_min_flow = std::min(path_min_flow, dd[i][u].zero_arc_flow);
                path.push_back(dd[i][u].zero_arc);
                label.push_back(zeroArc);
            }
        }
        if (color_used) {
            coloring.push_back(cclass);
            if (path_min_flow > 0) {
                for (int i = 0; i < n; i++) {
                    if (label[i]) {
                        dd[i][path[i]].one_arc_flow -= path_min_flow;
                    } else {
                        dd[i][path[i]].zero_arc_flow -= path_min_flow;
                    }
                }
            }
        }
    }
    if(int(coloring_selected_vertices.size()) < n){
        std::set<int> not_selected_vertices;
        std::set_difference(all_vertices.begin(), all_vertices.end(),
                            coloring_selected_vertices.begin(), coloring_selected_vertices.end(),
                            std::inserter(not_selected_vertices, not_selected_vertices.begin()));
        for(auto it = not_selected_vertices.begin(); it != not_selected_vertices.end();
            it = not_selected_vertices.erase(it)){ //return iterator following iterator that was deleted, acts as ++
            unsigned int i = *it;
            bool found_non_conflicting_cc;
            for(ColorClass& c : coloring){
                found_non_conflicting_cc = true;
                for(unsigned int j : c){
                    if(neighbors[i - 1].count(j)){
                        found_non_conflicting_cc = false;
                        break;
                    }
                }
                if(found_non_conflicting_cc){
                    c.insert(c.end(), i);
                    coloring_selected_vertices.insert(coloring_selected_vertices.end(), i);
                    break;
                }
            }
            if(not coloring_selected_vertices.count(i)) {
                //could not assign i to any cc, swap colors around or create new class to then assign it a colors
                bool recolored = false;//heuristic_try_color_swap(i, neighbors, coloring);
                if(not recolored){//TODO turn on again
                    coloring.push_back({i});
                } else {
//                        std::cout << "Swapped color in heuristic" << std::endl;
                }
                coloring_selected_vertices.insert(coloring_selected_vertices.end(), i);
            }

        }
    }
    return coloring;
}

int DDColors::heuristic_iterative_refinement() {

    Model model = (opt.relaxation == IP_Only) ? IP : LP;

    DecisionDiagram dd = initial_decision_diagram(graph);
    lower_bound = 0;
    upper_bound = heuristic_bound;

    longest_path_refinement(dd);


    int loops = 0;
    bool testing = false;
    while (lower_bound < upper_bound){
        if(testing and loops == 0) {
            std::cout << "How many loops?" << std::endl;
            std::cin >> loops;
            if(loops == 0) model = IP;
        }loops--;


        double flow_value = compute_flow_solution(dd, model, upper_bound); // = obj(F)
        int flow_bound = (model == IP) ? int(std::round(flow_value)) : (int)std::ceil(flow_value);
        if(model == IP) stats.num_ip_solved++; else stats.num_lp_solved++;
        if(flow_bound > lower_bound){
            update_lower_bound(flow_bound);
            if(lower_bound == upper_bound) break;
        }

//        double lp_flow = compute_flow_solution(dd, LP, upper_bound);
//        std::cout << "lp: " << lp_flow << " vs ip: " << flow_value << " gap of " << double(flow_value)/double(lp_flow) << std::endl;

        int coloring_size;
        std::vector<PathLabelConflict> conflict_info;

        if(model == IP) {
            coloring_size = int(primal_heuristic(dd).size());
            if(coloring_size < upper_bound){
                update_upper_bound(coloring_size);
                if(lower_bound == upper_bound) break;
            }
//            std::cout << "bounds: lower " << lower_bound << " upper " << upper_bound << " coloring " << coloring_size << std::endl;
            conflict_info = detect_edge_conflict(dd, neighbors, flow_value, model, opt.find_conflicts, opt.path_decomposition);
        }
        else if(model == LP){//distinction is necessary because combined function behaves differently in IP case
            std::tie(conflict_info, coloring_size) = find_conflict_and_primal_heuristic(dd, flow_value, model);
            if(coloring_size < upper_bound){
                update_upper_bound(coloring_size);
                if(lower_bound == upper_bound) break;
            }
//            std::cout << "bounds: lower " << lower_bound << " upper " << upper_bound << std::endl;
        }

        if(conflict_info.empty()) {
            if(model == LP){
                std::cout << "Switching from LP to IP since no conflict was found" << std::endl;
                model = IP;
                continue;
            } else break;
        }

        stats.num_conflicts+=conflict_info.size();
        for(const PathLabelConflict & plc : conflict_info){
            separate_edge_conflict(dd, neighbors, plc, opt.redirect_arcs);
        }

        if(opt.relaxation == Switch_LP_IP and model == IP){
            model = LP;
            std::cout << "Switching back from IP to LP" << std::endl;
        }
    }

    stats.decision_diagram_size = num_nodes(dd);
    stats.decision_diagram_arcs = num_arcs(dd);
    stats.decision_diagram_width = get_width(dd);
    return lower_bound;
}

//TODO update
int DDColors::heuristic_iterative_refinement_NUMERICALLY_SAFE() {

    Model model = (opt.relaxation == IP_Only) ? IP : LP;

    DecisionDiagram dd = initial_decision_diagram(graph);
    lower_bound = 0;
    upper_bound = heuristic_bound;

    longest_path_refinement(dd);


    while (lower_bound < upper_bound){

        double flow_value = compute_flow_solution(dd, model, upper_bound); // = obj(F)
        if(model == IP) stats.num_ip_solved++; else stats.num_lp_solved++;

        if((model == IP) and (int(flow_value) > lower_bound)) {
            //validate computed mip solution with SCIP
            double flow_value_scip = validate_flow_solution_SCIP(dd, flow_value, model, upper_bound);
            if(flow_value != flow_value_scip){
                throw std::runtime_error("Flow values do not match");
            }
        }

        lower_bound = std::max(lower_bound, (model == IP) ? int(flow_value) : int(std::ceil(flow_value - std::pow(10, -5))));

        int coloring_bound;
        std::vector<PathLabelConflict> conflict_info;

        if(model == IP) {
            coloring_bound = int(primal_heuristic(dd).size());
            upper_bound = std::min(upper_bound, coloring_bound);

            std::cout << "bounds: lower " << lower_bound << " upper " << upper_bound << std::endl;
            if(lower_bound == upper_bound) break;

            conflict_info = detect_edge_conflict(dd, neighbors, flow_value, model, opt.find_conflicts, PreferOneArcs);
        } else if(model == LP){//distinction is necessary because combined function behaves differently in IP case
            std::tie(conflict_info, coloring_bound) = find_conflict_and_primal_heuristic(dd, flow_value, model);
            upper_bound = std::min(upper_bound, coloring_bound);

            std::cout << "bounds: lower " << lower_bound << " upper " << upper_bound << std::endl;
            if(lower_bound == upper_bound) break;
        }

        if(conflict_info.empty()) {
            if(model == LP){
                std::cout << "Switching from LP to IP since no conflict was found" << std::endl;
                model = IP;
                continue;
            } else break;
        }

        stats.num_conflicts += conflict_info.size();
        for(const PathLabelConflict & plc : conflict_info){
            separate_edge_conflict(dd, neighbors, plc, opt.redirect_arcs);
        }

        if(opt.relaxation == Switch_LP_IP and model == IP){
            model = LP;
            std::cout << "Switching back from IP to LP" << std::endl;
        }
    }

    stats.decision_diagram_size = num_nodes(dd);
    stats.decision_diagram_arcs = num_arcs(dd);
    stats.decision_diagram_width = get_width(dd);
    return lower_bound;
}


void DDColors::update_lower_bound(int flow_bound) {
    lower_bound = flow_bound;

    std::chrono::steady_clock::time_point time_to_bound = std::chrono::steady_clock::now();
    stats.tt_lower_bound = (time_to_bound-stats.start_time);
    print_bounds();
}

void DDColors::update_upper_bound(int coloring_size){
    upper_bound = coloring_size;
    coloring_bound = coloring_size;

    std::chrono::steady_clock::time_point time_to_bound = std::chrono::steady_clock::now();
    stats.tt_upper_bound = (time_to_bound-stats.start_time);
    print_bounds();
}

void DDColors::print_bounds() const{
    if(true){ //TODO some condition on when to print this
        if(opt.algorithm == Options::BasicRefinement){
            std::cout << "bounds: lower " << lower_bound << " time: tt_lb " << stats.tt_lower_bound << std::endl;
        }else if(opt.algorithm == Options::HeuristicRefinement){
            std::cout << "bounds: lower " << lower_bound << (coloring_bound < heuristic_bound ? " upper  " : " upper* ")
            << upper_bound << " time: tt_lb " << stats.tt_lower_bound << " tt_ub " << stats.tt_upper_bound << std::endl;
        }
    }
}


void DDColors::longest_path_refinement(DecisionDiagram &dd) {
    stats.num_longest_path_conflicts = opt.num_longest_path_iterations;
    enum Mode {normal, experimental};
    Mode mode = normal;
    for(int longest_path_iteration = 0; longest_path_iteration < opt.num_longest_path_iterations; longest_path_iteration++){
        if(mode == normal){
            PathLabelConflict plc = conflict_on_longest_path(dd, neighbors);
            if(plc.path.empty()){
                std::cout << "Longest path procedure did not find a conflict" << std::endl;
                stats.num_longest_path_conflicts = longest_path_iteration;
                break;
            }
            separate_edge_conflict(dd, neighbors, plc, opt.redirect_arcs);
        }
        else if(mode == experimental){
            std::vector<PathLabelConflict> conflict_info = experimental_conflict_on_longest_path(dd, neighbors);
            if(conflict_info.empty()){
                std::cout << "Longest path procedure did not find a conflict" << std::endl;
                stats.num_longest_path_conflicts = longest_path_iteration;
                break;
            }
            for(PathLabelConflict & plc : conflict_info){
                separate_edge_conflict(dd, neighbors, plc, opt.redirect_arcs);
            }
        }
    }
    if(opt.num_longest_path_iterations){
        std::cout << "Finished longest path procedure, begin LP/IP refinement" << std::endl;
    }
}


//Extra stuff
std::pair<std::vector<PathLabelConflict>, int>
DDColors::find_conflict_and_primal_heuristic(DecisionDiagram &dd, double flow_val, Model model) {
    //  conflict detection
    //********************
    //  primal heuristic
    //(except for redefinitions)

    unsigned int n = num_vars(dd);
    std::vector<PathLabelConflict> conflict_info;
    std::vector<double> conflict_flow;
    bool continue_conflict_detection = true;
    //********************************************
    bool continue_primal_heuristic = true;
    Coloring coloring;
    bool color_used = true;
    std::set<unsigned int> coloring_selected_vertices;

    std::set<unsigned int> all_vertices;
    for (unsigned int i = 1; i <= n; ++i){
        all_vertices.insert(all_vertices.end(), i);
    }


    while((continue_conflict_detection) or (continue_primal_heuristic)) {
        Path path;
        Label label;
        std::vector<unsigned int> selected_vertices;
        bool path_found_conflict = false;
        path.push_back(0); //index of root node
        double path_min_flow = n;
        //**************************************
        ColorClass cclass;
        color_used = false;
        std::set<unsigned int> neighboring_selected_vertices;


        for(unsigned int i = 0; i < n; i++) {
            int u = path.back();
            //slight change in choosing arc here: in IP case only care about having greater flow than on 0-arc
            // but don't care whether 1-arc flow is at least 1
            if ((dd[i][u].one_arc != -1) and (not (opt.path_decomposition == PreferZeroArcs and (dd[i][u].zero_arc_flow >= dd[i][u].one_arc_flow)) ) and
                    ((model == IP and int(std::round(dd[i][u].one_arc_flow) >= 1)) or
                (model == LP and (dd[i][u].one_arc_flow >= dd[i][u].zero_arc_flow) )) )
            {
                path_min_flow = std::min(path_min_flow, dd[i][u].one_arc_flow);

                bool use_zero_arc_instead = false;
                if ((not path_found_conflict) and continue_conflict_detection) {
                    for (auto j_it = selected_vertices.rbegin(); j_it != selected_vertices.rend(); j_it++) {
                        unsigned int j = j_it.operator*();
                        if (neighbors[i].count(j)) {
                            conflict_info.emplace_back(
                                    Path(path.begin() + j - 1, path.end()),
                                    Label(label.begin() + j - 1, label.end()),
                                    std::make_tuple(j, i + 1)
                            );
                            if (opt.find_conflicts == SingleConflict) {
                                continue_conflict_detection = false;
                            }
                            //check if it's possible to take the 0-arc instead to avoid the conflict
                            if(dd[i][u].zero_arc_flow >= double_eps and opt.path_decomposition == AvoidConflicts){
                                path_min_flow = std::min(path_min_flow, dd[i][u].zero_arc_flow);
                                path.push_back(dd[i][u].zero_arc);
                                label.push_back(zeroArc);
                                use_zero_arc_instead = true;
                                break;
                            }
                            //otherwise continue but remember we found a conflict on this path already
                            path_found_conflict = true;
                            break;//don't return here but continue to search for conflicts
                        }
                    }
                }
                if (not use_zero_arc_instead) {
                    path.push_back(dd[i][u].one_arc);
                    label.push_back(oneArc);
                    selected_vertices.push_back(i + 1);
                }
                //********************************************
                if(continue_primal_heuristic){// and (not use_zero_arc_instead)
                    if ((neighboring_selected_vertices.count(i + 1) == 0) and (coloring_selected_vertices.count(i + 1) == 0)) {
                        cclass.insert(cclass.end(), i + 1);
                        coloring_selected_vertices.insert(coloring_selected_vertices.end(), i + 1);
                        neighboring_selected_vertices.insert(neighboring_selected_vertices.end(), i + 1);
                        neighboring_selected_vertices.insert(neighbors[i].begin(), neighbors[i].end());
                        color_used = true;
                    }
                }

            } else {
                path_min_flow = std::min(path_min_flow, dd[i][u].zero_arc_flow);
                path.push_back(dd[i][u].zero_arc);
                label.push_back(zeroArc);
            }
        }


        if (path_found_conflict and (opt.find_conflicts == LargestFlowConflict)) {
            conflict_flow.push_back(path_min_flow);
        }

        if (model == IP) {
            if (path_min_flow > 1 + double_eps) {
                throw std::runtime_error(
                        "Happened that path min flow is larger than 1 in IP case " + std::to_string(path_min_flow));
            }
        }

        if (path_min_flow > 0) {
            for (unsigned int i = 0; i < n; i++) {
                if (label[i]) {
                    //path_min_flow is just 1 if IP is used
                    dd[i][path[i]].one_arc_flow -= path_min_flow;
                } else {
                    dd[i][path[i]].zero_arc_flow -= path_min_flow;
                }
            }
            flow_val -= path_min_flow;
        } else if(continue_conflict_detection){
            throw std::runtime_error("Min flow on path was zero despite flow not being zero!");
        }

        continue_conflict_detection = continue_conflict_detection and
                                        ((model == IP and flow_val >= (1-double_eps)) or (model == LP and flow_val > double_eps));
        //*************************************************
        if (color_used) {
            coloring.push_back(cclass);
        }
        continue_primal_heuristic = (int(coloring_selected_vertices.size()) < n and color_used);

        if(not (continue_conflict_detection or continue_primal_heuristic)) {
            break;
        }
    }

//    std::cout << "found " << conflict_info.size() << " conflicts" << std::endl;


    //exception here: first do the primal heuristic stuff and the the code for LargestFlowConflict
    if(int(coloring_selected_vertices.size()) < n){
        std::set<int> not_selected_vertices;
        std::set_difference(all_vertices.begin(), all_vertices.end(),
                            coloring_selected_vertices.begin(), coloring_selected_vertices.end(),
                            std::inserter(not_selected_vertices, not_selected_vertices.begin()));
//            std::cout << "S: " << selected_vertices << " and not S: " << not_selected_vertices << std::endl;
        for(auto it = not_selected_vertices.begin(); it != not_selected_vertices.end();
            it = not_selected_vertices.erase(it)){ //return iterator following iterator that was deleted, acts as ++
            unsigned int i = *it;
            bool found_non_conflicting_cc;
            for(ColorClass& c : coloring){
                found_non_conflicting_cc = true;
                for(unsigned int j : c){
                    if(neighbors[i - 1].count(j)){
                        found_non_conflicting_cc = false;
                        break;
                    }
                }
                if(found_non_conflicting_cc){
                    c.insert(c.end(), i);
                    coloring_selected_vertices.insert(coloring_selected_vertices.end(), i);
                    break;
                }
            }
            if(not coloring_selected_vertices.count(i)) {
                //could not assign i to any cc, swap colors around or create new class to then assign it a colors
                bool recolored = heuristic_try_color_swap(i, neighbors, coloring);
                if(not recolored){
                    coloring.push_back({i});
                } else {
//                    std::cout << "Swapped color in heuristic" << std::endl;
                }
                coloring_selected_vertices.insert(coloring_selected_vertices.end(), i);
            }

        }
    }

    //*********************************************
    if(opt.find_conflicts == LargestFlowConflict){
        if(conflict_info.empty()){
            //no conflict was found
            return std::make_pair(std::vector<PathLabelConflict>{}, coloring.size());
        }
        std::vector<int> large_flow_indices;
        double largest_conflict_flow = *std::max_element(conflict_flow.begin(), conflict_flow.end());
        auto max_pos = std::find_if(conflict_flow.begin(), conflict_flow.end(),
                                    [&largest_conflict_flow](double flow){return flow >= (largest_conflict_flow * 0.9);});

        while(max_pos != conflict_flow.end()){
            large_flow_indices.push_back(int(max_pos - conflict_flow.begin()));
            max_pos = std::find_if(std::next(max_pos), conflict_flow.end(),
                                   [&largest_conflict_flow](double flow){return flow >= (largest_conflict_flow * 0.9);});
        }
//        std::cout << "Out of " << conflict_info.size() << " selected " << large_flow_indices.size() << " paths with flow " << largest_conflict_flow << " : " << conflict_flow << " indices " << large_flow_indices << std::endl;
        std::vector<PathLabelConflict> result;
        std::transform(large_flow_indices.begin(), large_flow_indices.end(), std::back_inserter(result),
                       [&conflict_info](int conflict_index){return conflict_info[conflict_index];});
        return std::make_pair(result, coloring.size());
    }

    return std::make_pair(conflict_info, coloring.size());
}


void DDColors::preprocessing_graph(int lower_bound) {
    int n = int(graph.ncount());
    std::cout << "original graph: " << graph.ncount() << ", " << graph.ecount() << std::endl;

    int clique_size = lower_bound;
   //preprocessing loop until graph doesn't change any further
    unsigned int num_nodes_previous_iteration = graph._ncount;
    while(true){
        num_nodes_previous_iteration = graph._ncount;

        //dominated vertices are not part of a clique, can remove this before searching for one
        graph.remove_dominated_vertices();
        //look for clique to get lower bound on chromatic number
        std::set<Vertex> clique = graph.find_clique(50); //TODO how many branches
        clique_size = std::max(clique_size, int(clique.size()));

        graph.peel_graph_ordered(clique_size);

        if(graph.ncount() == 0 or clique_size == heuristic_bound){
            //can return with clique size/lower bound
            //this already solves the graph
            break;
        }

        if(num_nodes_previous_iteration == graph._ncount)
            break;
    }
    if(graph.ncount() == n) throw std::runtime_error("no change");
    stats.start_time = std::chrono::steady_clock::now();//??
}


bool heuristic_try_color_swap(Vertex vertex, const NeighborList &neighbors, Coloring &coloring) {
    for(unsigned int j = 0; j < coloring.size(); j++){
        for(unsigned int k = j+1; k < coloring.size(); k++){
            std::set<Vertex> neighbors_v_intersect_color_j;
            std::set_intersection(neighbors[vertex - 1].begin(), neighbors[vertex - 1].end(),
                                  coloring[j].begin(), coloring[j].end(),
                                  std::inserter(neighbors_v_intersect_color_j, neighbors_v_intersect_color_j.end()));

            if(neighbors_v_intersect_color_j.size() == 1){
                Vertex u = *neighbors_v_intersect_color_j.begin();
                std::set<Vertex> neighbors_u_intersect_color_k;
                std::set_intersection(neighbors[u-1].begin(), neighbors[u-1].end(),
                                      coloring[k].begin(), coloring[k].end(),
                                      std::inserter(neighbors_u_intersect_color_k, neighbors_u_intersect_color_k.end()));

                if(neighbors_u_intersect_color_k.empty()){
                    //we can swap colors around, i.e. color max_saturated vertex and u with j and k respectively
                    coloring[j].erase(u);
                    coloring[j].insert(vertex);
                    coloring[k].insert(u);
                    k = coloring.size();
                    j = coloring.size();
                    return true;//was able to swap colors
                }
            }
        }
    }
    return false; //unable to swap two colors
}
