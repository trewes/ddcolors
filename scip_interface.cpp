#include "scip_interface.h"


double compute_flow_solution_SCIP(DecisionDiagram &dd, Model model, int coloring_upper_bound) {
    if(model == LP){
        std::cout << "Warning: using SCIP as an LP_solver, make sure this is what you want." << std::endl;
    }

    SCIP* scip;
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL( SCIPcreateProbBasic(scip, "scip") );
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));
    SCIP_CALL( SCIPsetLongintParam(scip, "constraints/countsols/sollimit", 1) );

    unsigned int n = num_vars(dd);
    int ub = (coloring_upper_bound == -1) ? int(n) : coloring_upper_bound - 1; //have flow at most chi(G)-1 on edges

    //add variables for all edges in decision diagram but only set obj to 1 for outgoing edges of root
    EdgeIndex edge_index = 0;
    unsigned int n_arcs = num_arcs(dd);
    std::vector<SCIP_VAR*> scip_vars(n_arcs);
    auto var_type = (model == IP) ? SCIP_VARTYPE_INTEGER: SCIP_VARTYPE_CONTINUOUS;
    for(int layer = 0; layer < n; layer++){
        for(Node &u : dd[layer]){
            int obj = layer ? 0 : 1;
            //First add one_arc variables and then zero_arc variables. keep this ordering in mind!
            if(u.one_arc != -1) { //bound 1-arcs by 1
                SCIP_CALL( SCIPcreateVarBasic(scip, &scip_vars[edge_index], ("var_"+std::to_string(edge_index)).c_str(), 0.0, double(ub), double(obj), var_type));
                SCIP_CALL( SCIPaddVar(scip, scip_vars[edge_index]) );
                edge_index++;
            }
            //another variable for zero arc
            SCIP_CALL( SCIPcreateVarBasic(scip, &scip_vars[edge_index], ("var_"+std::to_string(edge_index)).c_str(), 0.0, double(ub), double(obj), var_type));
            SCIP_CALL( SCIPaddVar(scip, scip_vars[edge_index]) );
            edge_index++;
        }
    }
    // introduced all variables and set objective function but have no rows/constraints
    //add rows in two phases: 1. in each level there is exactly one 1-arc with flow 1
    // 2. flow conservation and //3. that variables are integer in range is set during adding those variables/columns

    //1.
    edge_index = 0; //count/track edges while iterating over them
    for(unsigned int layer = 0; layer < n; layer++){
        std::vector<SCIP_VAR*> vec_ind;
        for(Node &u : dd[layer]) {
            //get index of the 1-arc of that node if it exists and set that cval to 1
            if(u.one_arc != -1){
                vec_ind.push_back(scip_vars[edge_index]);
                edge_index++;
            }
            if(u.zero_arc != -1){
                edge_index++;
            } else { //this should never happen since we don't consider the terminal node
                std::cout << "Error: there is no zero-arc!";
            }
        }
        std::vector<double> vec_val(vec_ind.size(), 1.0); //this is just a vector full of 1.0 the size of vec_ind
        int num_one_arcs = int(vec_ind.size());

        SCIP_CONS* constraint = nullptr;
        SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, ("constraint_level_"+std::to_string(layer)).c_str(), num_one_arcs, vec_ind.data(), vec_val.data(), 1, SCIPinfinity(scip)));//greater equal 1
        SCIP_CALL( SCIPaddCons(scip, constraint));
        SCIP_CALL( SCIPreleaseCons(scip, &constraint) );
    }

    //2.
    edge_index = 0;
    int node_index = 0;
    //map of nodes to edge index
    std::map<NodeIndex , std::vector<EdgeIndex> > old_incoming_arcs;
    Node & root = dd[0][0];
    if(root.one_arc != -1){
        old_incoming_arcs[dd[1][root.one_arc].index].push_back(edge_index);
        edge_index++;
    }
    if(root.zero_arc != -1){
        old_incoming_arcs[dd[1][root.zero_arc].index].push_back(edge_index);
        edge_index++;
    }
    std::map<NodeIndex , std::vector<EdgeIndex> > new_incoming_arcs;
    for(unsigned int layer = 1; layer < n; layer++){
        //find all incoming arcs of u and set cval to 1, for outgoing which are easy to find set cval to -1
        for(Node &u : dd[layer]){
            std::vector<SCIP_VAR*> vec_ind;
            std::vector<double> vec_val;
            //incoming edges
            for(EdgeIndex e : old_incoming_arcs[u.index]){
                vec_ind.push_back(scip_vars[e]);
                vec_val.push_back(1.0);
            }
            //outgoing edges
            if(u.one_arc != -1){
                new_incoming_arcs[dd[layer+1][u.one_arc].index].push_back(edge_index);
                vec_ind.push_back(scip_vars[edge_index]);
                vec_val.push_back(-1.0);
                edge_index++;
            }
            if(u.zero_arc != -1){
                new_incoming_arcs[dd[layer+1][u.zero_arc].index].push_back(edge_index);
                vec_ind.push_back(scip_vars[edge_index]);
                vec_val.push_back(-1.0);
                edge_index++;
            } else { //this should never happen since we don't consider the terminal node
                std::cout << "Error: there is no zero-arc!" << std::endl;
            }

            int node_degree = int(vec_ind.size());
            SCIP_CONS* constraint = nullptr;
            SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, ("constraint_node_"+std::to_string(node_index)).c_str(), node_degree, vec_ind.data(), vec_val.data(), 0, 0));
            SCIP_CALL( SCIPaddCons(scip, constraint));
            SCIP_CALL( SCIPreleaseCons(scip, &constraint) );
            node_index++;
        }
        old_incoming_arcs = new_incoming_arcs;
        new_incoming_arcs.clear();
    }

    //TODO how does this work
//    if(model == IP) SCIPupdateCutoffbound(scip, double(ub+1));

    SCIP_VAR* root_arcs[2] = {scip_vars[0], scip_vars[1]};
    double root_coeff[2] = {1.0,1.0};
    SCIP_CONS* constraint = nullptr;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "upper_bound_optimal_value", 2, root_arcs, root_coeff, 0, double(ub+1)));
    SCIP_CALL( SCIPaddCons(scip, constraint));
    SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

    SCIP_CALL( (SCIPwriteOrigProblem(scip, "test_scip.lp", nullptr, FALSE)));
    SCIP_CALL( SCIPsolve(scip) );

    double flow_val = 0;
    flow_val = SCIPgetPrimalbound(scip);
    SCIP_SOL* solution;
    solution = SCIPgetBestSol(scip);

    edge_index = 0;
    for(unsigned int layer = 0; layer < n; layer++) {
        for (Node &u : dd[layer]) {
            if(u.one_arc != -1){
                u.one_arc_flow = SCIPgetSolVal(scip, solution, scip_vars[edge_index]);
                SCIP_CALL( SCIPreleaseVar(scip, &scip_vars[edge_index]));
                edge_index++;
            }
            if(u.zero_arc != -1){
                u.zero_arc_flow = SCIPgetSolVal(scip, solution, scip_vars[edge_index]);
                SCIP_CALL( SCIPreleaseVar(scip, &scip_vars[edge_index]));
                edge_index++;
            }
        }
    }

    SCIP_CALL( SCIPfree(&scip) );
    return flow_val;
}


double validate_flow_solution_SCIP(DecisionDiagram &dd, double flow_value, Model model, int coloring_upper_bound) {
    if(model == LP){
        std::cout << "Warning: using SCIP as an LP_solver, make sure this is what you want." << std::endl;
    }

    SCIP* scip;
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL( SCIPcreateProbBasic(scip, "scip") );
    SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

    SCIP_CALL( SCIPsetRealParam(scip, "solvingphases/optimalvalue", flow_value) );
//    SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_FEASIBILITY, 1) );
    SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, 1) );

    unsigned int n = num_vars(dd);
    int ub = (coloring_upper_bound == -1) ? int(n) : coloring_upper_bound - 1; //have flow at most chi(G)-1 on edges

    //add variables for all edges in decision diagram but only set obj to 1 for outgoing edges of root
    EdgeIndex edge_index = 0;
    unsigned int n_arcs = num_arcs(dd);
    std::vector<SCIP_VAR*> scip_vars(n_arcs);
    std::vector<double> given_flow(n_arcs);
    auto var_type = (model == IP) ? SCIP_VARTYPE_INTEGER: SCIP_VARTYPE_CONTINUOUS;
    for(int layer = 0; layer < n; layer++){
        for(Node &u : dd[layer]){
            int obj = layer ? 0 : 1;
            //First add one_arc variables and then zero_arc variables. keep this ordering in mind!
            if(u.one_arc != -1) { //bound 1-arcs by 1
                SCIP_CALL( SCIPcreateVarBasic(scip, &scip_vars[edge_index], ("var_"+std::to_string(edge_index)).c_str(), 0.0, double(ub), double(obj), var_type));
                SCIP_CALL( SCIPaddVar(scip, scip_vars[edge_index]) );
                given_flow[edge_index] = u.one_arc_flow;
                edge_index++;
            }
            //another variable for zero arc
            SCIP_CALL( SCIPcreateVarBasic(scip, &scip_vars[edge_index], ("var_"+std::to_string(edge_index)).c_str(), 0.0, double(ub), double(obj), var_type));
            SCIP_CALL( SCIPaddVar(scip, scip_vars[edge_index]) );
            given_flow[edge_index] = u.zero_arc_flow;
            edge_index++;
        }
    }
    // introduced all variables and set objective function but have no rows/constraints
    //add rows in two phases: 1. in each level there is exactly one 1-arc with flow 1
    // 2. flow conservation and //3. that variables are integer in range is set during adding those variables/columns

    //1.
    edge_index = 0; //count/track edges while iterating over them
    for(unsigned int layer = 0; layer < n; layer++){
        std::vector<SCIP_VAR*> vec_ind;
        for(Node &u : dd[layer]) {
            //get index of the 1-arc of that node if it exists and set that cval to 1
            if(u.one_arc != -1){
                vec_ind.push_back(scip_vars[edge_index]);
                edge_index++;
            }
            if(u.zero_arc != -1){
                edge_index++;
            } else { //this should never happen since we don't consider the terminal node
                std::cout << "Error: there is no zero-arc!";
            }
        }
        std::vector<double> vec_val(vec_ind.size(), 1.0); //this is just a vector full of 1.0 the size of vec_ind
        int num_one_arcs = int(vec_ind.size());

        SCIP_CONS* constraint = nullptr;
        SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, ("constraint_level_"+std::to_string(layer)).c_str(), num_one_arcs, vec_ind.data(), vec_val.data(), 1, SCIPinfinity(scip)));//greater equal 1
        SCIP_CALL( SCIPaddCons(scip, constraint));
        SCIP_CALL( SCIPreleaseCons(scip, &constraint) );
    }

    //2.
    edge_index = 0;
    int node_index = 0;
    //map of nodes to edge index
    std::map<NodeIndex , std::vector<EdgeIndex> > old_incoming_arcs;
    Node & root = dd[0][0];
    if(root.one_arc != -1){
        old_incoming_arcs[dd[1][root.one_arc].index].push_back(edge_index);
        edge_index++;
    }
    if(root.zero_arc != -1){
        old_incoming_arcs[dd[1][root.zero_arc].index].push_back(edge_index);
        edge_index++;
    }
    std::map<NodeIndex , std::vector<EdgeIndex> > new_incoming_arcs;
    for(unsigned int layer = 1; layer < n; layer++){
        //find all incoming arcs of u and set cval to 1, for outgoing which are easy to find set cval to -1
        for(Node &u : dd[layer]){
            std::vector<SCIP_VAR*> vec_ind;
            std::vector<double> vec_val;
            //incoming edges
            for(EdgeIndex e : old_incoming_arcs[u.index]){
                vec_ind.push_back(scip_vars[e]);
                vec_val.push_back(1.0);
            }
            //outgoing edges
            if(u.one_arc != -1){
                new_incoming_arcs[dd[layer+1][u.one_arc].index].push_back(edge_index);
                vec_ind.push_back(scip_vars[edge_index]);
                vec_val.push_back(-1.0);
                edge_index++;
            }
            if(u.zero_arc != -1){
                new_incoming_arcs[dd[layer+1][u.zero_arc].index].push_back(edge_index);
                vec_ind.push_back(scip_vars[edge_index]);
                vec_val.push_back(-1.0);
                edge_index++;
            } else { //this should never happen since we don't consider the terminal node
                std::cout << "Error: there is no zero-arc!" << std::endl;
            }

            int node_degree = int(vec_ind.size());
            SCIP_CONS* constraint = nullptr;
            SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, ("constraint_node_"+std::to_string(node_index)).c_str(), node_degree, vec_ind.data(), vec_val.data(), 0, 0));
            SCIP_CALL( SCIPaddCons(scip, constraint));
            SCIP_CALL( SCIPreleaseCons(scip, &constraint) );
            node_index++;
        }
        old_incoming_arcs = new_incoming_arcs;
        new_incoming_arcs.clear();
    }


    //TODO how does this work
//    if(model == IP) SCIPupdateCutoffbound(scip, double(ub+1));


//    SCIP_VAR* root_arcs[2] = {scip_vars[0], scip_vars[1]};
//    double root_coeff[2] = {1.0,1.0};
//    SCIP_CONS* constraint = nullptr;
//    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &constraint, "upper_bound_optimal_value", 2, root_arcs, root_coeff, 0, double(ub+1)));
//    SCIP_CALL( SCIPaddCons(scip, constraint));
//    SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

//    constraint for first 1-arc being at least 1
//    COLORlp_addrow(flow_lp, 1, root_arcs, root_coeff, COLORlp_GREATER_EQUAL, double(1), nullptr);


    SCIP_Sol* given_solution;
    SCIP_CALL( SCIPcreateOrigSol(scip, &given_solution, NULL) );
    SCIP_CALL( SCIPsetSolVals(scip, given_solution, int(n_arcs), scip_vars.data(), given_flow.data()) );
    unsigned int stored;
    SCIP_CALL( SCIPaddSolFree(scip, &given_solution, &stored) );

    SCIP_CALL( (SCIPwriteOrigProblem(scip, "test_scip.lp", nullptr, FALSE)));
    SCIP_CALL( SCIPsolve(scip) );

    double flow_val = SCIPgetPrimalbound(scip);
    SCIP_SOL* solution;
    solution = SCIPgetBestSol(scip);
    std::cout << "Done with solving the " << (model==LP ? "LP" : "IP") << ", obj: " << std::setprecision(20) << flow_val << std::endl;



    edge_index = 0;
    for(unsigned int layer = 0; layer < n; layer++) {
        for (Node &u : dd[layer]) {
            if(u.one_arc != -1){
                u.one_arc_flow = SCIPgetSolVal(scip, solution, scip_vars[edge_index]);
                SCIP_CALL( SCIPreleaseVar(scip, &scip_vars[edge_index]));
                edge_index++;
            }
            if(u.zero_arc != -1){
                u.zero_arc_flow = SCIPgetSolVal(scip, solution, scip_vars[edge_index]);
                SCIP_CALL( SCIPreleaseVar(scip, &scip_vars[edge_index]));
                edge_index++;
            }
        }
    }

    SCIP_CALL( SCIPfree(&scip) );
    return flow_val;
}