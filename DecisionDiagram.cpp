#include "DecisionDiagram.h"

template<typename T>
std::ostream &operator<<(std::ostream &s, const std::set<T> &set) {
    s << "{";
    std::string sep;
    for(T el : set){
        s << sep << el;
        sep = ", ";
    }
    return s << "}";
}

template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec) {
    s << "[";
    std::string sep;
    for(T el : vec){
        s << sep << el;
        sep = ", ";
    }
    return s << "]";
}


Node::Node() = default;

Node::Node(int layer, NodeIndex index, std::set<unsigned int> state_info)
        : layer(layer), index(index), state_info(std::move(state_info)) {}

bool Node::equivalent(const Node &other) const {
    return this->state_info == other.state_info;
}

unsigned int num_vars(const DecisionDiagram &dd) {
    return dd.size() - 1;
}

unsigned int num_nodes(const DecisionDiagram &dd) {
    return std::accumulate(dd.begin(), dd.end(), 0,
                           [](size_t sum, const std::vector<Node> &level) { return sum + level.size(); });
}

unsigned int num_arcs(const DecisionDiagram &dd) {
    unsigned int num_arcs = 0;
    for(const auto &layer : dd){
        for(const auto &node : layer){
            if(node.one_arc != -1){
                num_arcs++;
            }
            if(node.zero_arc != -1){
                num_arcs++;
            }
        }
    }
    return num_arcs;
}

unsigned int get_width(const DecisionDiagram &dd) {
    return std::max_element(dd.begin(), dd.end(),
                            [](const std::vector<Node> &layer1, const std::vector<Node> &layer2) {
                                return layer1.size() < layer2.size();
                            })->size();
}

void print_decision_diagram(const DecisionDiagram &dd, bool use_tag) {
    for(const auto &layer : dd){
        std::cout << (use_tag ? "[DD] " : "") << "Layer: " << layer.front().layer << " of size: " << layer.size()
                  << std::endl;
        for(const auto &node : layer){
            std::cout << (use_tag ? "[DD] " : "") << "Node " << node.index << " state " << node.state_info
                      << " with 1-arc " << node.one_arc << " and 0-arc " << node.zero_arc << std::endl;
        }
    }
}

PathLabelConflict::PathLabelConflict(Path path, Label label, Conflict conflict)
        : path(std::move(path)), label(std::move(label)), conflict(std::move(conflict)) {}


DecisionDiagram exact_decision_diagram(const Graph &g) {
    NeighborList neighbors = g.get_neighbor_list();
    return exact_decision_diagram(g, neighbors);
}

DecisionDiagram exact_decision_diagram(const Graph &g, const NeighborList &neighbors) {
    DecisionDiagram dd;
    dd.resize(g.ncount() + 1);
    int num_nodes = 0;
    int num_arcs = 0;

    StateInfo initial_state;
    for(unsigned int i = 1; i <= g.ncount(); ++i){
        initial_state.insert(initial_state.end(), i);
    }

    dd[0].emplace_back(1, 0, initial_state);
    num_nodes++;
    for(unsigned int layer = 0; layer < g.ncount(); ++layer){
        for(Node &u : dd[layer]){
            if(num_nodes > std::pow(10, 6)){
                throw std::runtime_error("The exact decision diagram contains more than one million nodes.");
            }

            if(u.state_info.count(layer + 1)){
                auto new_state = u.state_info;
                new_state.erase(layer + 1);
                for(unsigned int neighbor : neighbors[layer]){
                    new_state.erase(neighbor);
                }

                //look if there exists an equivalent node
                auto one_arc_it = std::find_if(dd[layer + 1].begin(), dd[layer + 1].end(),
                                               [&new_state](const Node &v) { return v.state_info == new_state; });
                if(one_arc_it == dd[layer + 1].end()){
                    dd[layer + 1].emplace_back(layer + 2, int(dd[layer + 1].size()), new_state);
                    num_nodes++;
                    u.one_arc = int(dd[layer + 1].size() - 1);
                } else{
                    u.one_arc = int(one_arc_it - dd[layer + 1].begin());
                }
                num_arcs++;

                new_state = u.state_info;
                new_state.erase(layer + 1);
                auto zero_arc_it = std::find_if(dd[layer + 1].begin(), dd[layer + 1].end(),
                                                [&new_state](const Node &v) { return v.state_info == new_state; });
                if(zero_arc_it == dd[layer + 1].end()){
                    dd[layer + 1].emplace_back(layer + 2, int(dd[layer + 1].size()), new_state);
                    num_nodes++;
                    u.zero_arc = int(dd[layer + 1].size() - 1);
                } else{
                    u.zero_arc = int(zero_arc_it - dd[layer + 1].begin());
                }
                num_arcs++;
            } else{
                auto zero_arc_it = std::find_if(dd[layer + 1].begin(), dd[layer + 1].end(),
                                                [&u](const Node &v) { return v.state_info == u.state_info; });
                if(zero_arc_it == dd[layer + 1].end()){
                    dd[layer + 1].emplace_back(layer + 2, int(dd[layer + 1].size()), u.state_info);
                    num_nodes++;
                    u.zero_arc = int(dd[layer + 1].size() - 1);
                } else{
                    u.zero_arc = int(zero_arc_it - dd[layer + 1].begin());
                }
                num_arcs++;
            }
        }
    }
    std::cout << "Done building the exact decision diagram containing " << num_nodes << " nodes, " << num_arcs
              << " arcs and width " << get_width(dd) << "." << std::endl;
    return dd;
}

DecisionDiagram initial_decision_diagram(const Graph &g) {
    DecisionDiagram dd;
    dd.resize(g.ncount() + 1);

    StateInfo new_state;
    for(unsigned int i = 1; i <= g.ncount(); ++i){
        new_state.insert(new_state.end(), i);
    }

    dd[0].emplace_back(1, 0, new_state);
    for(unsigned int j = 0; j < g.ncount(); ++j){
        new_state.erase(j + 1);
        dd[j + 1].emplace_back(j + 2, 0, new_state);
        dd[j][0].zero_arc = 0;
        dd[j][0].one_arc = 0;
    }
    return dd;
}


int intersection_size(const StateInfo &a, const StateInfo &b) {
    if(a.empty() or b.empty()) return 0;
    return int(std::count_if(a.begin(), a.end(), [&b](const unsigned int &k) { return b.count(k); }));
}

void separate_edge_conflict(DecisionDiagram &dd, const NeighborList &neighbors, const PathLabelConflict &plc,
                            RedirectArcs redirect_arcs) {
    //check assumptions? of conflict being minimal

    Path path = plc.path;
    const Label &label = plc.label;
    const Conflict &conflict = plc.conflict;

    int j, k;
    std::tie(j, k) = conflict;
    for(int i = 0; i <= (k - 1) - j; ++i){
        int layer = j + i - 1;
        StateInfo new_state = dd[layer][path[i]].state_info;
        new_state.erase(layer + 1);
        if(label[i]){
            for(unsigned int n : neighbors[layer]){
                new_state.erase(n);
            }
        }

        int t = -1;
        auto search_equivalent_node = std::find_if(dd[layer + 1].begin(), dd[layer + 1].end(),
                                                   [&new_state](const Node &v) { return v.state_info == new_state; });
        if(search_equivalent_node == dd[layer + 1].end()){
            Node w(layer + 2, int(dd[layer + 1].size()), new_state);
            if(redirect_arcs == OriginalArcs){
                if(w.state_info.count(layer + 1 + 1)){
                    w.one_arc = dd[layer + 1][path[i + 1]].one_arc;
                }
                w.zero_arc = dd[layer + 1][path[i + 1]].zero_arc;
            } else if(redirect_arcs == MostSimilarNode){
                StateInfo next_state = w.state_info;
                next_state.erase(layer + 1 + 1);
                auto most_similar_node = std::max_element(dd[layer + 2].begin(), dd[layer + 2].end(),
                                                          [&next_state](const Node &u, const Node &v) {
                                                              return intersection_size(next_state, u.state_info) <
                                                                     intersection_size(next_state, v.state_info);
                                                          });
                w.zero_arc = int(most_similar_node - dd[layer + 2].begin());

                if(w.state_info.count(layer + 1 + 1)){
                    for(unsigned int n : neighbors[layer + 1]){
                        next_state.erase(n);
                    }
                    most_similar_node = std::max_element(dd[layer + 2].begin(), dd[layer + 2].end(),
                                                         [&next_state](const Node &u, const Node &v) {
                                                             return intersection_size(next_state, u.state_info) <
                                                                    intersection_size(next_state, v.state_info);
                                                         });
                    w.one_arc = int(most_similar_node - dd[layer + 2].begin());
                }

            }
            dd[layer + 1].push_back(w);
            t = int(dd[layer + 1].size() - 1);

        } else{
            t = int(search_equivalent_node - dd[layer + 1].begin());
        }

        if(label[i]){
            dd[layer][path[i]].one_arc = t;
        } else{
            dd[layer][path[i]].zero_arc = t;
        }
        path[i + 1] = t;
    }
}

std::vector<PathLabelConflict>
detect_edge_conflict(DecisionDiagram dd, const NeighborList &neighbors, double flow_val, Model model,
                     ConflictResolution find_conflicts, PathDecomposition path_decomposition) {
    unsigned int n = num_vars(dd);
    std::vector<PathLabelConflict> conflict_info;
    std::vector<double> conflict_flow;

    while((model == IP and flow_val >= (1 - double_eps)) or (model == LP and flow_val > double_eps)){
        Path path;
        Label label;
        std::vector<unsigned int> selected_vertices;
        bool path_found_conflict = false;
        path.push_back(0); //index of root node
        double path_min_flow = (double) n;

        for(unsigned int i = 0; i < n; i++){
            int u = path.back();
            if((dd[i][u].one_arc != -1) and
               (not(path_decomposition == PreferZeroArcs and (dd[i][u].zero_arc_flow >= dd[i][u].one_arc_flow))) and
               ((model == IP and (int(std::round(dd[i][u].one_arc_flow) >= 1))) or
                (model == LP and dd[i][u].one_arc_flow >= dd[i][u].zero_arc_flow))){
                path_min_flow = std::min(path_min_flow, dd[i][u].one_arc_flow);

                bool use_zero_arc_instead = false;
                if(not path_found_conflict){
                    for(auto j_it = selected_vertices.rbegin();
                        j_it != selected_vertices.rend(); j_it++){//TODO what about selected and neighbors?
                        unsigned int j = j_it.operator*();
                        if(neighbors[i].count(j)){
                            conflict_info.emplace_back(
                                    Path(path.begin() + j - 1, path.end()),
                                    Label(label.begin() + j - 1, label.end()),
                                    std::make_tuple(j, i + 1)
                            );
                            if(find_conflicts == SingleConflict){
                                return conflict_info;
                            }
                            //check if it's possible to take the 0-arc instead to avoid the conflict
                            if(dd[i][u].zero_arc_flow >= double_eps and path_decomposition == AvoidConflicts){
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
                if(not use_zero_arc_instead){
                    path.push_back(dd[i][u].one_arc);
                    label.push_back(oneArc);
                    selected_vertices.push_back(i + 1);
                }
            } else{
                path_min_flow = std::min(path_min_flow, dd[i][u].zero_arc_flow);
                path.push_back(dd[i][u].zero_arc);
                label.push_back(zeroArc);
            }
        }

        if(path_found_conflict and (find_conflicts == LargestFlowConflict)){
            conflict_flow.push_back(path_min_flow);
        }

        if(model == IP){
//            std::cout << "IP and min flow on path was " << path_min_flow << std::endl;
            if(path_min_flow > 1 + double_eps){
                throw std::runtime_error(
                        "Happened that path min flow is larger than 1 in IP case " + std::to_string(path_min_flow));
            }
        } else if(path_min_flow > 1 + double_eps){
            //TODO remove throwing here
            throw std::runtime_error(
                    "Happened that path min flow is larger than 1 in IP case " + std::to_string(path_min_flow));
        }

        if(path_min_flow > 0){
            for(int i = 0; i < n; i++){
                if(label[i]){
                    dd[i][path[i]].one_arc_flow -= path_min_flow;
                } else{
                    dd[i][path[i]].zero_arc_flow -= path_min_flow;
                }
            }
            flow_val -= path_min_flow;
        } else{
            throw std::runtime_error("Min flow on path was zero despite flow not being zero!");
        }
    }
//    std::cout << "found " << conflict_info.size() << " conflicts" << std::endl;
    if(conflict_info.empty()){
        return conflict_info;
    }

    if(find_conflicts == LargestFlowConflict){
        double largest_up_to = std::pow(10, -3);
        std::vector<int> large_flow_indices;
        double largest_conflict_flow = *std::max_element(conflict_flow.begin(), conflict_flow.end());

        auto max_pos = std::find_if(conflict_flow.begin(), conflict_flow.end(),
                                    [&largest_conflict_flow, &largest_up_to](double flow) {
                                        return flow >= (largest_conflict_flow - largest_up_to);
                                    });
        while(max_pos != conflict_flow.end()){
            large_flow_indices.push_back(int(max_pos - conflict_flow.begin()));
            max_pos = std::find_if(std::next(max_pos), conflict_flow.end(),
                                   [&largest_conflict_flow, &largest_up_to](double flow) {
                                       return flow >= (largest_conflict_flow - largest_up_to);
                                   });
        }
//        std::cout << "Out of " << conflict_info.size() << " selected " << large_flow_indices.size() << " paths with flow " << largest_conflict_flow << " : " << conflict_flow << " indices " << large_flow_indices << std::endl;
        std::vector<PathLabelConflict> result;
        std::transform(large_flow_indices.begin(), large_flow_indices.end(), std::back_inserter(result),
                       [&conflict_info](int conflict_index) { return conflict_info[conflict_index]; });
        return result;
    }
    return conflict_info;
}


double compute_flow_solution(DecisionDiagram &dd, Model model, int coloring_upper_bound, Formulation formulation) {

    COLORlp *flow_lp;//environment is initialised in DDColors
    COLORlp_init(&flow_lp, "flow_lp");
    COLORlp_objective_sense(flow_lp, COLORlp_MIN); //minimizing should be the default setting anyway

    unsigned int n = num_vars(dd);
    int ub = (coloring_upper_bound == -1) ? int(n) : coloring_upper_bound - 1; //have flow at most chi(G)-1 on edges


    //add variables for all edges in decision diagram but only set obj to 1 for outgoing edges of root
    EdgeIndex edge_index = 0;
    auto var_type = (model == IP) ? COLORlp_INTEGER : COLORlp_CONTINUOUS;
    for(int layer = 0; layer < n; layer++){
        for(Node &u : dd[layer]){
            int obj = layer ? 0 : 1;
//            var_type = (layer) ? COLORlp_CONTINUOUS: COLORlp_INTEGER;//TODO test
            //First add one_arc variables and then zero_arc variables. keep this ordering in mind!
            if(u.one_arc != -1){ //bound 1-arcs by 1
                COLORlp_addcol(flow_lp, 0, (int *) nullptr, (double *) nullptr, obj,
                               0.0, double(ub), var_type, nullptr);
                edge_index++;
            }
            //another variable for zero arc
            COLORlp_addcol(flow_lp, 0, (int *) nullptr, (double *) nullptr, obj,
                           0.0, double(ub), var_type, nullptr);
            edge_index++;
        }
    }
    int num_columns = edge_index;
    // introduced all variables and set objective function but have no rows/constraints

    //add rows in two phases: 1. in each level there is exactly one 1-arc with flow 1
    // 2. flow conservation and //3. that variables are integer in range is set during adding those variables/columns

    //1.
    edge_index = 0; //count/track edges while iterating over them
    for(unsigned int layer = 0; layer < n; layer++){
        std::vector<NodeIndex> vec_ind;
        for(Node &u : dd[layer]){
            //get index of the 1-arc of that node if it exists and set that cval to 1
            if(u.one_arc != -1){
                vec_ind.push_back(edge_index);
                edge_index++;
            }
            if(u.zero_arc != -1){
                edge_index++;
            } else{ //this should never happen since we don't consider the terminal node
                std::cout << "Error: there is no zero-arc!";
            }
        }
        std::vector<double> vec_val(vec_ind.size(), 1.0); //this is just a vector full of 1.0 the size of vec_ind
        int num_one_arcs = int(vec_ind.size());

        //TODO also consider this in the formulation of the dual and thus for the safe bounds. That just means that the variables z_j are >= 0 instead of unrestricted
//        COLORlp_addrow(flow_lp, num_one_arcs, vec_ind.data(), vec_val.data(), COLORlp_EQUAL, 1, nullptr);
        COLORlp_addrow(flow_lp, num_one_arcs, vec_ind.data(), vec_val.data(), COLORlp_GREATER_EQUAL, 1, nullptr);
    }

    //2.
    edge_index = 0;
    int node_index = 0;
    //map of nodes to edge index
    std::map<NodeIndex, std::vector<EdgeIndex> > old_incoming_arcs;
    Node &root = dd[0][0];
    if(root.one_arc != -1){
        old_incoming_arcs[dd[1][root.one_arc].index].push_back(edge_index);
        edge_index++;
    }
    if(root.zero_arc != -1){
        old_incoming_arcs[dd[1][root.zero_arc].index].push_back(edge_index);
        edge_index++;
    }
    std::map<NodeIndex, std::vector<EdgeIndex> > new_incoming_arcs;
    for(unsigned int layer = 1; layer < n; layer++){
        //find all incoming arcs of u and set cval to 1, for outgoing which are easy to find set cval to -1
        for(Node &u : dd[layer]){
            std::vector<NodeIndex> vec_ind;
            std::vector<double> vec_val;
            //incoming edges
            for(EdgeIndex e : old_incoming_arcs[u.index]){
                vec_ind.push_back(e);
                vec_val.push_back(1.0);
            }
            //outgoing edges
            if(u.one_arc != -1){
                new_incoming_arcs[dd[layer + 1][u.one_arc].index].push_back(edge_index);
                vec_ind.push_back(edge_index);
                vec_val.push_back(-1.0);
                edge_index++;
            }
            if(u.zero_arc != -1){
                new_incoming_arcs[dd[layer + 1][u.zero_arc].index].push_back(edge_index);
                vec_ind.push_back(edge_index);
                vec_val.push_back(-1.0);
                edge_index++;
            } else{ //this should never happen since we don't consider the terminal node
                std::cout << "Error: there is no zero-arc!" << std::endl;
            }

            int node_degree = int(vec_ind.size());
            COLORlp_addrow(flow_lp, node_degree, vec_ind.data(), vec_val.data(), COLORlp_EQUAL, 0, nullptr);
            node_index++;
        }
        old_incoming_arcs = new_incoming_arcs;
        new_incoming_arcs.clear();
    }


    if(formulation == ExtraConstraints or formulation == AllConstraints){
        std::vector<int> index = {0};
        std::vector<double> val = {1.0};
        COLORlp_addrow(flow_lp, 1, index.data(), val.data(), COLORlp_EQUAL, 1, nullptr);//or COLORlp_GREATER_EQUAL
    }
    if(formulation == BoundConstraints or formulation == AllConstraints){
        //input variable bounds as extra constraints
        std::vector<int> index = {0};
        std::vector<double> val = {1.0};
        for(int i = 0; i < num_columns; i++){
            if(i == 0 and formulation == AllConstraints) continue;
            index = {i};
            COLORlp_addrow(flow_lp, 1, index.data(), val.data(), COLORlp_LESS_EQUAL, ub, nullptr);
        }
    }


    //TODO how does this work
//    if(model == IP) COLORlp_set_cutoff(flow_lp, double(ub+1));


    COLORlp_optimize(flow_lp);
    double flow_val = 0;
    COLORlp_objval(flow_lp, &flow_val);
//todo hint at code
//    std::cout << "Done with solving the " << (model==LP ? "LP" : "IP") << ", obj: " << std::setprecision(20) << flow_val << std::endl;

    std::vector<double> solution;
    solution.reserve(num_columns);
    COLORlp_x(flow_lp, solution.data());
    solution.assign(solution.data(), solution.data() + num_columns);

    double double_safe_bound = 0.0;
    if(model == LP){

        if("Held" and false){
            std::vector<double> dual_solution;
            unsigned int dd_size = num_nodes(dd);
            dual_solution.reserve(n + dd_size - 2 + num_columns);
            COLORlp_pi(flow_lp, dual_solution.data());
            dual_solution.assign(dual_solution.data(), dual_solution.data() + n + dd_size - 2);

            double dual_obj = std::accumulate(dual_solution.begin(), dual_solution.begin() + n, double(0.0),
                                              [](double sum, double var) { return sum + var; });
            std::cout << "dual obj: " << std::setprecision(20) << dual_obj << "\n" << dual_solution << std::endl;

            const int originalRounding = fegetround();
            fesetround(FE_TOWARDZERO);

            using Ktype = long long;
            Ktype K = std::min(-(std::numeric_limits<long>::min() + 1), std::numeric_limits<long>::max()) /
                      (n * num_columns);
            Ktype K2 = std::pow(2, std::floor(std::log2(K)));
            std::cout << -(std::numeric_limits<Ktype>::min() + 1) << " and " << std::numeric_limits<Ktype>::max()
                      << " -> " << K << std::endl;
            std::vector<Ktype> int_dual_solution(dual_solution.size());
            std::transform(dual_solution.begin(), dual_solution.end(), int_dual_solution.begin(),
                           [&K](double var) { return std::floor(K * var); });
            std::cout << "max of dual solution is " << *std::max_element(dual_solution.begin(), dual_solution.end())
                      << " with sum " << std::setprecision(200) << dual_solution[0] * K << std::endl;
            std::cout << int_dual_solution << std::endl;
            Ktype int_dual_obj = std::accumulate(int_dual_solution.begin(), int_dual_solution.begin() + n, Ktype(0),
                                                 [](Ktype sum, Ktype var) { return sum + var; });
            long double scaled_int_dual_obj = int_dual_obj / (long double) K;
            std::cout << int_dual_obj << " and reduced again: " << std::setprecision(20) << scaled_int_dual_obj
                      << std::endl;
            if(int_dual_obj < K * flow_val) std::cout << "test here" << std::endl;
            else
                std::cout << "other test here" << std::endl;
            std::cout << "test alternative: " << std::setprecision(20) << int_dual_obj << " vs " << K * flow_val
                      << std::endl;
            if(flow_val < scaled_int_dual_obj)
                std::cout << std::setprecision(200) << flow_val << "!=" << scaled_int_dual_obj
                          << " Not matching obj and dual obj values: " << std::abs(flow_val - scaled_int_dual_obj)
                          << std::endl;
            if(flow_val < scaled_int_dual_obj)
                throw std::runtime_error("Not matching obj and scaled int dual obj values: " +
                                         std::to_string(std::abs(flow_val - scaled_int_dual_obj)));
            fesetround(originalRounding);
        }

        if("Neumaier" and true){

            const int originalRounding = fegetround();

            std::vector<double> dual_solution;
            unsigned int dd_size = num_nodes(dd);
            // dual solution consists of n variabes for the constraints of each layer + the variables for each node except r,t
            // + possibly dual variables for the bound constraints, if those are used
            dual_solution.reserve(n + dd_size - 2 + num_columns);
            COLORlp_pi(flow_lp, dual_solution.data());
            dual_solution.assign(dual_solution.data(),
                                 dual_solution.data() + n + dd_size - 2);//consider bound constraints! //TODO

            fesetround(FE_DOWNWARD);
            long double dual_obj = std::accumulate(dual_solution.begin(), dual_solution.begin() + n,
                                                   (long double) (0.0),
                                                   [](long double sum, double var) { return sum + var; });

            fesetround(FE_UPWARD);

            //save total index of first node on each level
            //! shifted by one since root node has no dual variable
            std::vector<int> shifted_node_indices;
            shifted_node_indices.push_back(-1); // -1 to shift the indices
            for(int i = 1; i < int(dd.size()); i++){
                int next_level_index = shifted_node_indices[i - 1] + int(dd[i - 1].size());
                shifted_node_indices.push_back(next_level_index);
            }

            std::vector<long double> residual;
            residual.reserve(num_columns);
            for(int layer = 0; layer < n; layer++){
                for(int index = 0; index < int(dd[layer].size()); index++){
                    Node &u = dd[layer][index];
                    if(u.one_arc != -1){
                        long double res_val;
                        if(layer == 0){
                            res_val = (long double) 1 * dual_solution[layer] - 1 +
                                      1 * dual_solution[n + shifted_node_indices[layer + 1] + u.one_arc];
                        } else if(layer == n - 1){
                            res_val = (long double) 1 * dual_solution[layer] -
                                      1 * dual_solution[n + shifted_node_indices[layer] + index] + 0;
                        } else{
                            res_val = (long double) 1 * dual_solution[layer] -
                                      1 * dual_solution[n + shifted_node_indices[layer] + index] +
                                      1 * dual_solution[n + shifted_node_indices[layer + 1] + u.one_arc];
                        }
                        residual.push_back(res_val);
                    }
                    if(u.zero_arc != -1){
                        long double res_val;
                        if(layer == 0){
                            res_val = (long double) 0 * dual_solution[layer] - 1 +
                                      1 * dual_solution[n + shifted_node_indices[layer + 1] + u.zero_arc];
                        } else if(layer == n - 1){
                            res_val = (long double) 0 * dual_solution[layer] -
                                      1 * dual_solution[n + shifted_node_indices[layer] + index] + 0;
                        } else{
                            res_val = (long double) 0 * dual_solution[layer] -
                                      1 * dual_solution[n + shifted_node_indices[layer] + index] +
                                      1 * dual_solution[n + shifted_node_indices[layer + 1] + u.zero_arc];
                        }
                        residual.push_back(res_val);
                    }
                }
            }

            if(num_columns != int(residual.size())){
                throw std::runtime_error("Not the correct sizes: " + std::to_string(num_columns) + " and " +
                                         std::to_string(residual.size()));
            }

            //compute delta with set upper bounds for primal variables
            long double delta = 0;
            for(long double &res_val : residual){
                delta += ub * std::max(res_val, (long double) (0.0));
            }

            fesetround(FE_DOWNWARD);

            long double safe_bound = dual_obj - delta;
            double_safe_bound = (double) safe_bound;//return this double as the safe lower bound

            fesetround(originalRounding);

//            std::cout << "'safe' lower bound is " << std::setprecision(20) << safe_bound << " with appx. dual bound " << dual_obj << " and delta " << delta  << "\nflow value was " << flow_val << std::endl;
            if(delta > std::pow(10, -6)){
                std::cout << "Warning: delta is quite big, larger than 10^-6: " << std::setprecision(25) << delta
                          << std::endl;
                throw std::runtime_error("Delta");
            }
            if(std::ceil(double_safe_bound) != std::ceil(flow_val - std::pow(10, -5))){
                std::stringstream message;
                message << std::setprecision(25) << "ceil of flow_value and double_safe_bound were different: " +
                                                    std::to_string(flow_val - std::pow(10, -5)) + " and " +
                                                    std::to_string(double_safe_bound);
                throw std::runtime_error(message.str());
            }
        }
    }

    //assign the flow to the edges of the decision diagram
    edge_index = 0;
    for(unsigned int layer = 0; layer < n; layer++){
        for(Node &u : dd[layer]){
            if(u.one_arc != -1){
                u.one_arc_flow = solution[edge_index];
                edge_index++;
            }
            if(u.zero_arc != -1){
                u.zero_arc_flow = solution[edge_index];
                edge_index++;
            }
        }
    }

    COLORlp_free(&flow_lp);

    if(model == LP){
        return double_safe_bound;
    }
    return flow_val;
}

void find_longest_path(const DecisionDiagram &dd, Path &path, Label &label) {
    int n = int(dd.size());
    //use that dd is dag and find longest path (longest in regards to 1-arcs) using topological order of the dd
    std::vector<std::vector<int> > prev(n);
    std::vector<std::vector<int> > dist(n);
    for(int layer = 0; layer < n; layer++){
        dist[layer] = std::vector<int>(dd[layer].size(), std::numeric_limits<int>::min());
        prev[layer].resize(dd[layer].size(), -1);
    }
    dist[0][0] = 0; //source vertex
    for(int layer = 0; layer < n - 1; layer++){
        for(int node = 0; node < int(dd[layer].size()); node++){
            if((dd[layer][node].one_arc != -1) and (dist[layer + 1][dd[layer][node].one_arc] < dist[layer][node] + 1)){
                dist[layer + 1][dd[layer][node].one_arc] = dist[layer][node] + 1;
                prev[layer + 1][dd[layer][node].one_arc] = node;
            }
            if(dist[layer + 1][dd[layer][node].zero_arc] < dist[layer][node]){
                dist[layer + 1][dd[layer][node].zero_arc] = dist[layer][node];
                prev[layer + 1][dd[layer][node].zero_arc] = node;
            }
        }
    }
    path.clear();
    path.reserve(n);
    label.clear();
    label.reserve(n);

    int prev_node = 0;
    path.push_back(prev_node);
    for(int layer = n - 1; layer > 0; layer--){
        prev_node = prev[layer][prev_node];
        if(dd[layer - 1][prev_node].one_arc == path.back()){
            label.push_back(oneArc);
        } else{
            label.push_back(zeroArc);
        }
        path.push_back(prev_node);
    }

    std::reverse(path.begin(), path.end());
    std::reverse(label.begin(), label.end());
}

PathLabelConflict conflict_on_longest_path(const DecisionDiagram &dd, const NeighborList &neighbors) {
    Path path;
    Label label;
    find_longest_path(dd, path, label);

    //    std::cout << path << " and " << label << std::endl;

    //look for a conflict on said longest path
    for(int j = 1; j < int(path.size()); j++){
        if(not label[j]){
            continue;
        }
        //iterating this way we find the first smallest conflict, important assumption for Algorithm 1
        for(int i = j - 1; i >= 0; i--){
            if(not label[i]){
                continue;
            }
            //both vertex i and j of original graph are chosen. check if that is a conflict or not
            if(neighbors[i].count(j + 1)){
//                std::cout << "plc: " << std::vector<int>(path.begin()+i, path.begin()+j +1) << " " << std::vector<bool>(label.begin()+i, label.begin()+j) << " " << i+1 << "," << j+1 << std::endl;
                return PathLabelConflict(Path(path.begin() + i, path.begin() + j + 1),
                                         Label(label.begin() + i, label.begin() + j), std::make_tuple(i + 1, j + 1));
            }
        }
    }
    std::cout << "no more conflicts found on longest path" << std::endl;
    return {{}, {}, std::make_tuple(-1, -1)};
}

std::vector<PathLabelConflict>
experimental_conflict_on_longest_path(const DecisionDiagram &dd, const NeighborList &neighbors) {
    int n = int(dd.size());
    //use that dd is dag and find longest path (longest in regards to 1-arcs) using topological order of the dd
    std::vector<std::vector<std::vector<int> > > prev(n);
    std::vector<std::vector<int> > dist(n);
    for(int layer = 0; layer < n; layer++){
        dist[layer] = std::vector<int>(dd[layer].size(), std::numeric_limits<int>::min());
        prev[layer].resize(dd[layer].size(), std::vector<int>{-1});
    }
    dist[0][0] = 0; //source vertex
    for(int layer = 0; layer < n - 1; layer++){
        for(int node = 0; node < int(dd[layer].size()); node++){
            if((dd[layer][node].one_arc != -1)){
                if(dist[layer + 1][dd[layer][node].one_arc] < dist[layer][node] + 1){
                    dist[layer + 1][dd[layer][node].one_arc] = dist[layer][node] + 1;
                    prev[layer + 1][dd[layer][node].one_arc].clear();
                    prev[layer + 1][dd[layer][node].one_arc].push_back(node);
                } else if(dist[layer + 1][dd[layer][node].one_arc] == dist[layer][node] + 1){
                    prev[layer + 1][dd[layer][node].one_arc].push_back(node);
                }
            }
            if(dist[layer + 1][dd[layer][node].zero_arc] < dist[layer][node]){
                dist[layer + 1][dd[layer][node].zero_arc] = dist[layer][node];
                prev[layer + 1][dd[layer][node].zero_arc].clear();
                prev[layer + 1][dd[layer][node].zero_arc].push_back(node);
            } else if(dist[layer + 1][dd[layer][node].zero_arc] == dist[layer][node]){
                prev[layer + 1][dd[layer][node].zero_arc].push_back(node);
            }
        }
    }


    std::vector<Path> paths;
//    paths.reserve(n);
    std::vector<Label> labels;
//    label.reserve(n);

    int prev_node = 0;
    paths.push_back({prev_node});
    labels.emplace_back();
    for(unsigned int layer = n - 1; layer > 0; layer--){
        for(unsigned int i = 0, num_paths = paths.size(), new_num_paths = num_paths; i < num_paths; i++){
            int last_node = paths[i].back();
            for(int more_prevs = 1;
                more_prevs < prev[layer][last_node].size(); more_prevs++){ //for each additional prev copy current path
                paths.push_back(paths[i]);
                labels.push_back(labels[i]);
            }

            prev_node = prev[layer][last_node].front();     //deal with first (if not only) prev node with the original path and thus different index
            if(dd[layer - 1][prev_node].one_arc == paths[i].back()){
                labels[i].push_back(oneArc);
            } else{
                labels[i].push_back(zeroArc);
            }
            paths[i].push_back(prev_node);

            for(int current_prev = 1; current_prev < prev[layer][last_node].size(); current_prev++){
                prev_node = prev[layer][last_node][current_prev];
                if(dd[layer - 1][prev_node].one_arc == paths[new_num_paths - 1 + current_prev].back()){
                    labels[new_num_paths - 1 + current_prev].push_back(oneArc);
                } else{
                    labels[new_num_paths - 1 + current_prev].push_back(zeroArc);
                }
                paths[new_num_paths - 1 + current_prev].push_back(prev_node);
            }
            new_num_paths += prev[layer][last_node].size() - 1;
        }
    }
    for(int i = 0; i < paths.size(); i++){
        std::reverse(paths[i].begin(), paths[i].end());
        std::reverse(labels[i].begin(), labels[i].end());
    }
//    std::cout << path << " and " << label << std::endl;

    std::cout << "found " << paths.size() << " paths" << std::endl;
    for(int path = 0; path < paths.size(); path++){
//        std::cout << "A path is: " << paths[path]  << " with length " << std::accumulate(labels[path].begin(), labels[path].end(), int(0))<< std::endl;
    }

    std::vector<PathLabelConflict> conflict_info;
    std::set<Conflict> conflicts;
    std::set<int> already_used_vertices;
    //look for a conflict on said longest paths
    for(int path = 0; path < paths.size(); path++){
        for(int j = 1; j <
                       int(paths[path].size()); j++){ //this way we find the first smallest conflict, important assumption for Algorithm 1
            if((not labels[path][j]) or already_used_vertices.count(j + 1)){
                continue;
            }
            for(int i = j - 1; i >= 0; i--){
                if((not labels[path][i]) or already_used_vertices.count(i + 1)){
                    continue;
                }
                //both vertex i and j of original graph are chosen. check if that is a conflict or not
                if(neighbors[i].count(j + 1) and (not conflicts.count(std::make_tuple(i + 1, j + 1)))){
                    already_used_vertices.insert(i + 1);
                    already_used_vertices.insert(j + 1);
//                  std::cout << "plc: " << std::vector<int>(path.begin()+i, path.begin()+j +1) << " " << std::vector<bool>(label.begin()+i, label.begin()+j) << " " << i+1 << "," << j+1 << std::endl;
                    conflict_info.emplace_back(Path(paths[path].begin() + i, paths[path].begin() + j + 1),
                                               Label(labels[path].begin() + i, labels[path].begin() + j),
                                               std::make_tuple(i + 1, j + 1)
                    );
                    conflicts.insert(std::make_tuple(i + 1, j + 1));
                    i = 0;
                    j = int(paths[path].size()); //escape both loops, continue with next path
                }
            }
        }
    }
    std::cout << conflict_info.size() << " paths had conflicts" << std::endl;
    return conflict_info;
}


DecisionDiagram exact_decision_diagram_test(const Graph &g, const NeighborList &neighbors) {
    DecisionDiagram dd;
    dd.resize(g.ncount() + 1);
    int num_nodes = 0;
    int num_arcs = 0;

    StateInfo initial_state;
    for(unsigned int i = 1; i <= g.ncount(); ++i){
        initial_state.insert(initial_state.end(), i);
    }

    dd[0].emplace_back(1, 0, initial_state);
    num_nodes++;
    for(unsigned int layer = 0; layer < g.ncount(); ++layer){
        //which variable to chosse next?
        std::map<Vertex, int> vertex_min_states;
        for(Node &u : dd[layer]){
            for(Vertex v : u.state_info){
                vertex_min_states[v]++;
            }
        }
        Vertex next_var = std::min_element(vertex_min_states.begin(), vertex_min_states.end(),
                                           [](std::pair<Vertex, int> u, std::pair<Vertex, int> w) {
                                               return u.second < w.second;
                                           })->first;
        if(next_var == 0) throw std::runtime_error("Next vertex was 0, this can not happen.");
        std::cout << "Choose vertex/var " << next_var << " appearing in " << vertex_min_states[next_var] << " states"
                  << std::endl;
        for(Node &u : dd[layer]){

            if(num_nodes > std::pow(10, 6)){
                throw std::runtime_error("The exact decision diagram contains more than one million nodes.");
            }

            if(u.state_info.count(next_var)){
                auto new_state = u.state_info;
                new_state.erase(next_var);
                for(unsigned int neighbor : neighbors[next_var - 1]){
                    new_state.erase(neighbor);
                }

                //look if there exists an equivalent node
                auto one_arc_it = std::find_if(dd[layer + 1].begin(), dd[layer + 1].end(),
                                               [&new_state](const Node &v) { return v.state_info == new_state; });
                if(one_arc_it == dd[layer + 1].end()){
                    dd[layer + 1].emplace_back(layer + 2, int(dd[layer + 1].size()), new_state);
                    num_nodes++;
                    u.one_arc = int(dd[layer + 1].size() - 1);
                } else{
                    u.one_arc = int(one_arc_it - dd[layer + 1].begin());
                }
                num_arcs++;

                new_state = u.state_info;
                new_state.erase(next_var);
                auto zero_arc_it = std::find_if(dd[layer + 1].begin(), dd[layer + 1].end(),
                                                [&new_state](const Node &v) { return v.state_info == new_state; });
                if(zero_arc_it == dd[layer + 1].end()){
                    dd[layer + 1].emplace_back(layer + 2, int(dd[layer + 1].size()), new_state);
                    num_nodes++;
                    u.zero_arc = int(dd[layer + 1].size() - 1);
                } else{
                    u.zero_arc = int(zero_arc_it - dd[layer + 1].begin());
                }
                num_arcs++;
            } else{
                auto zero_arc_it = std::find_if(dd[layer + 1].begin(), dd[layer + 1].end(),
                                                [&u](const Node &v) { return v.state_info == u.state_info; });
                if(zero_arc_it == dd[layer + 1].end()){
                    dd[layer + 1].emplace_back(layer + 2, int(dd[layer + 1].size()), u.state_info);
                    num_nodes++;
                    u.zero_arc = int(dd[layer + 1].size() - 1);
                } else{
                    u.zero_arc = int(zero_arc_it - dd[layer + 1].begin());
                }
                num_arcs++;
            }
        }
    }
    std::cout << "Done building the exact decision diagram containing " << num_nodes << " nodes, " << num_arcs
              << " arcs and width " << get_width(dd) << "." << std::endl;
    return dd;
}


