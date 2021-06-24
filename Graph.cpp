#include "Graph.h"


template<typename T>
std::ostream& operator<<(std::ostream& s, const std::set<T> & set){
    s << "{";
    std::string sep;
    for(T el : set){
        s << sep << el;
        sep = ", ";
    }
    return s << "}";
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T> & vec){
    s << "[";
    std::string sep;
    for(T el : vec){
        s << sep << el;
        sep = ", ";
    }
    return s << "]";
}

template<typename T> T random_element(std::set<T> const &v){
    auto it = v.cbegin();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, v.size()-1);
    int random = distr(gen);
    std::advance(it, random);
    return *it;
}

Permutation identity(size_t size){
    Permutation identity(size);
    for(int i = 1; i <= size; i++)
        identity[i-1] = i;
    return identity;
}

Permutation perm_inverse(const Permutation& perm) {
    Permutation inverse(perm.size());
    for(int i = 1; i <= perm.size(); i++){
        inverse[perm[i-1]-1] = i;
    }
    return inverse;
}


NeighborList Graph::get_neighbor_list() const {
    NeighborList list;
    list.resize(_ncount);
    for(size_t i = 0; i < _ecount; i++){
        list[_elist[2 * i] - 1].insert(list[_elist[2 * i] - 1].end(), _elist[2 * i + 1]);
        list[_elist[2 * i + 1] - 1].insert(list[_elist[2 * i + 1] - 1].end(), _elist[2 * i]);
    }
    return list;
}

Graph::Graph(unsigned int ncount, unsigned int ecount, std::vector<Vertex> elist) :
            _ncount(ncount), _ecount(ecount), _elist(std::move(elist)){}

unsigned int Graph::ncount() const {
    return _ncount;
}

unsigned int Graph::ecount() const {
    return _ecount;
}

void Graph::print() const {
    std::cout << "Graph has " << _ncount << " vertices and " << _ecount << " edges." << std::endl;
    NeighborList neighbors = get_neighbor_list();
    for(std::set<Vertex>& adj : neighbors){
        std::cout << adj << std::endl;
    }
}

std::vector<Vertex> Graph::elist() const{
    return _elist;
}
std::vector<int> Graph::shifted_elist() const{
    std::vector<int> shifted_elist(2*_ecount);
    std::transform(_elist.begin(), _elist.end(), shifted_elist.begin(),
                   [](Vertex v){if(v == 0) throw std::runtime_error("Vertex is 0, cannot shift it."); return (int)(v-1);});
    return shifted_elist;
}



Graph::Graph(const char *filename) : _ncount(0), _ecount(0){
    std::string string_name = filename;
    if(string_name.substr(string_name.find_last_of('.')+1) == "col"){
        read_dimacs(filename);
    } else if(string_name.substr(string_name.find_last_of('.')+1) == "g6"){
        read_graph6(filename);
    } else if (string_name.substr(string_name.find_last_of('.')+1) == "b"){
        throw std::runtime_error("Binary dimacs format is currently not supported, please translate it yourself.");
    } else {
        throw std::runtime_error("File format is not accepted.");
    }
}

Graph::Graph(std::string g6_string) : _ncount(0), _ecount(0){
    read_graph6_string(std::move(g6_string));
}

Graph::~Graph(){
//    delete opt;
}

void Graph::read_dimacs(const char *filename) {
    std::ifstream file(filename);
    if (not file) {
        throw std::runtime_error("Cannot open file.");
    }

    std::string line;
    std::getline(file, line);//read in the file line by line

    while(line[0] == 'c' or line[0] == '\r' or line[0] == '\0'){          //special character to account for empty lines
        std::getline(file, line);                                                               //ignore comments
    }
    if(line[0] != 'p'){
        throw std::runtime_error("File is not in correct dimacs format.");          //a different case was expected
    }
    std::string p, edge;
    int n, e;
    std::stringstream ss(line);
    ss >> p >> edge >> n >> e;                                                       //read in number of nodes and edges

    _ncount = n;
    _ecount = e;
    _elist.reserve(2 * _ecount);

    std::getline(file, line);
    while(line[0] == 'n'){
        std::cout << "The dimacs file is specifying colors of vertices, these will be ignored for this problem." << std::endl;
        std::getline(file, line);
    }

    char first;
    int head, tail;
    do {
        while(line[0] == 'c' or line[0] == '\r' or line[0] == '\0'){          //special character to account for empty lines
            std::getline(file, line);                                                               //ignore comments
        }
        ss.clear();
        ss.str(std::string());
        ss << line;
        ss >> first >> head >> tail;
        _elist.push_back(head);
        _elist.push_back(tail);
    }
    while (std::getline(file, line));
}

void Graph::read_graph6(const char *filename) {
    std::ifstream file(filename);
    if (not file) {
        throw std::runtime_error("Cannot open file.");
    }
    std::cout << "Warning: if the graph6 file contains more than one graph, the first is read in" << std::endl;
    std::string line;
    std::getline(file, line);
    read_graph6_string(line);
}

void Graph::read_graph6_string(std::string g6_string) {// Extract vertex count from graph6
    if (g6_string.substr(0, 1) == ":") {
        g6_string = g6_string.substr(1, g6_string.size());
    }

    int vertexCount = -1;

    int r0 = g6_string[0];
    int adjIdx = 1;

    // single byte expansion (0 <= n <= 62)
    if (r0 >= 0 + 63 and r0 <= 63 + 62) {
        vertexCount = r0 - 63;
        adjIdx = 1;
    }
        // four or eight byte expansions
    else if (r0 == 126) {
        char r1 = g6_string[1];

        // eight byte expansions (258048 <= n <= 68719476735)
        if (r1 == 126) {
            adjIdx = 8;
            throw std::runtime_error("Terribly sorry, but we don\t do 8-byte expansions yet :(");
        }
            // four byte expansion (63 <= n <= 258047)
        else {
            if (g6_string.size() < 4) {
                throw std::runtime_error("This input seems to be too short for a 4-byte expansion :(");
            } else {
                vertexCount = 0;
                adjIdx = 4;

                for (int i = 0; i < 3; i++) {
                    vertexCount = vertexCount | ((g6_string[1 + i] - 63) << ((2 - i) * 6));
                }
            }
        }
    }

    // Only start working when we actually have a vertex count
    if (vertexCount > -1) {
        _ncount = vertexCount;
        int edgeCount = 0;
        // Extract adjacency list from graph6,
        for (int n = 1; n < vertexCount; n++) {
            // Natural numbers sum n*(n+1)/2, but only need previous n so (n-1)*(n+1-1)/2
            int nPos = n * (n - 1) / 2;
            for (int m = 0; m < n; m++) {
                int bitPos = nPos + m;
                int charIndex = adjIdx + int(std::floor(bitPos / 6));
                int charVal = g6_string[charIndex] - 63;

                int bitIndex = bitPos % 6;
                int bitShift = (5 - bitIndex);
                int bitMask = (1 << bitShift);

                int bitVal = ((charVal & bitMask) >> bitShift);

                // If the bit for this particular position is 1 then {m, n} is an edge
                if (bitVal == 1) {
                    _elist.push_back(n + 1);
                    _elist.push_back(m + 1);
                    edgeCount++;
                }
            }
        }
        _ecount = edgeCount;
    } else {
        throw std::runtime_error("Was not able to get vertex count from graph6 file.");
    }
}



Graph Graph::perm_graph(const Permutation &perm) const {
    if(perm.size() != _ncount){
        throw std::runtime_error("Size of graph and permutation do not match.");
    }
    Graph perm_graph(_ncount, _ecount, {});
    perm_graph._elist.resize(2*_ecount);
    for(size_t i = 0; i < 2*_ecount; i++){
        perm_graph._elist[i] = perm[_elist[i]-1];
    }
    return perm_graph;
}

Coloring Graph::dsatur(Permutation &ordering) const {
    NeighborList neighbors = get_neighbor_list();
    return dsatur(ordering, neighbors);
}

Coloring Graph::dsatur(Permutation &ordering, const NeighborList& neighbors) const {
    ordering.clear();
    ordering.reserve(_ncount);
    Coloring coloring;
    //initialise all vertices to be uncolored
    std::vector<int> vertex_color(_ncount, -1);
    std::set<Vertex> unselected_vertices;
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++) {
        unselected_vertices.insert(unselected_vertices.begin(), vertex);
    }

    //find max degree vertex and select it as first colored vertex. ties are broken randomly
    Vertex max_degree_vertex = 0;
    std::set<Vertex> tied_max_degree;
    int max_degree = -1;
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++){
        if(int(neighbors[vertex-1].size()) > max_degree){
            max_degree = int(neighbors[vertex-1].size());
            max_degree_vertex = vertex;
            tied_max_degree.clear();
            tied_max_degree.insert(tied_max_degree.end(), max_degree_vertex);
        }
        else if(int(neighbors[vertex-1].size()) == max_degree){
            tied_max_degree.insert(tied_max_degree.end(), vertex);
        }
    }
    //if more than one vertex with max degree, choose one at random
    if(tied_max_degree.size() > 1 and opt!= NULL and opt->random_tie_breaks){
        max_degree_vertex = random_element(tied_max_degree);
    }

    //color first selected vertex
    ordering.push_back(max_degree_vertex);
    unselected_vertices.erase(max_degree_vertex);
    coloring.push_back({}); //new color class
    coloring[0].insert(max_degree_vertex); //assign first color to max degree vertex
    vertex_color[max_degree_vertex-1] = 0;

    //record for an uncolored vertex to how many colored vertices it is connected
    std::vector<int> saturation_level(_ncount, 0);
    //don't want to select already selected vertex
    saturation_level[max_degree_vertex - 1] = std::numeric_limits<int>::min();
    //update saturation level for neighbors of first selected vertex
    for(Vertex neighbor : neighbors[max_degree_vertex-1]){
        saturation_level[neighbor - 1]++;
    }

    //main loop looking for the next vertex according to its saturation level
    while(not unselected_vertices.empty()){
        Vertex max_saturated_vertex = 0;
        int max_saturated = -1;
        //keep set of tied vertices to later choose one at random
        std::set<Vertex> saturated_tied_max_degrees;
        for(Vertex vertex : unselected_vertices){
            if(saturation_level[vertex - 1] > max_saturated){
                //found vertex with higher saturation, update
                max_saturated = saturation_level[vertex - 1];
                max_saturated_vertex = vertex;

                saturated_tied_max_degrees.clear();
                saturated_tied_max_degrees.insert(vertex);
            } else if(saturation_level[vertex - 1] == max_saturated){
                if(neighbors[vertex - 1].size() > neighbors[max_saturated_vertex - 1].size()) {
                    //found vertex with equal saturation but higher degree, update
                    max_saturated_vertex = vertex;

                    saturated_tied_max_degrees.clear();
                    saturated_tied_max_degrees.insert(vertex);
                }
                else if(neighbors[vertex - 1].size() == neighbors[max_saturated_vertex - 1].size()){
                    //same saturation level, same degree, add to set of possible candidates
                    saturated_tied_max_degrees.insert(vertex);
                }
            }
        }
        //select max saturated vertex, at random if there exist multiple
        if(saturated_tied_max_degrees.size() > 1 and opt!= NULL and opt->random_tie_breaks) {
            max_saturated_vertex = int(random_element(saturated_tied_max_degrees));
        }
        ordering.push_back(max_saturated_vertex);
        unselected_vertices.erase(max_saturated_vertex);

        //get all the colors of colored neighbors from the chosen vertex, ignore uncolored neighbors
        std::set<int> saturation_colors;
        saturation_colors.clear();
        for(Vertex neighbor : neighbors[max_saturated_vertex-1]){
            if(vertex_color[neighbor - 1] != -1)
                saturation_colors.insert(vertex_color[neighbor-1]);
        }

        //find lowest color not used by any neighbors
        int lowest_color = 0;
        for(int color = 0; color <= int(saturation_colors.size()); color++){
            if(saturation_colors.count(color) == 0){
                lowest_color = color;
                break;
            }
        }

        if(lowest_color < int(coloring.size())){
            //we can add chosen vertex to some existing color class
            coloring[lowest_color].insert(max_saturated_vertex);
            vertex_color[max_saturated_vertex-1] = lowest_color;
        } else if(lowest_color == int(coloring.size())){
            //new color class is created for max_saturated_vertex
            //before that we try to recolor two vertices to not have to use a new color
            bool recolored = try_color_swap(max_saturated_vertex, neighbors, coloring, vertex_color);
            if(not recolored){
                //not able to recolor, add new color class
                coloring.push_back({});
                coloring[lowest_color].insert(max_saturated_vertex);
                vertex_color[max_saturated_vertex-1] = lowest_color;
            }
        } else {
            std::cout << "Error: lowest color is too high, something must have gone wrong." << std::endl;
        }

        //new vertex has been selected so update the saturation levels of all its neighbors if it's a different color
        int msv_color = vertex_color[max_saturated_vertex - 1];
        for(Vertex vertex : neighbors[max_saturated_vertex - 1]){
            bool is_new_color = true;
            //possibly check that neighbor has been selected, i.e. it's level is equal to numeric_limits<int>::min()
            if(saturation_level[vertex - 1] != std::numeric_limits<int>::min()){
                for(Vertex adj_vertex : neighbors[vertex - 1]){
                    if (adj_vertex != max_saturated_vertex and vertex_color[adj_vertex - 1] == msv_color){
                        //found neighbor of neighbor that uses msv color so we don't increase the saturation level
                        is_new_color = false;
                        break;
                    }
                }
                if(is_new_color)
                    saturation_level[vertex - 1]++;
            }
        }
        //remove selected vertex from any further consideration
        saturation_level[max_saturated_vertex - 1] = std::numeric_limits<int>::min();
    }
    //colored all vertices
    return coloring;
}


Coloring Graph::dsatur_original(Permutation &ordering) const {
    NeighborList neighbors = get_neighbor_list();
    return dsatur_original(ordering, neighbors);
}


Coloring Graph::dsatur_original(Permutation &ordering, const NeighborList& neighbors) const {
    ordering.clear();
    ordering.reserve(_ncount);
    Coloring coloring;
    //initialise all vertices to be uncolored
    std::vector<int> vertex_color(_ncount, -1);
    std::set<Vertex> unselected_vertices;
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++) {
        unselected_vertices.insert(unselected_vertices.begin(), vertex);
    }

    //find max degree vertex and select it as first colored vertex. ties are broken randomly
    Vertex max_degree_vertex = 0;
    std::set<Vertex> tied_max_degree;
    int max_degree = -1;
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++){
        if(int(neighbors[vertex-1].size()) > max_degree){
            max_degree = int(neighbors[vertex-1].size());
            max_degree_vertex = vertex;
            tied_max_degree.clear();
            tied_max_degree.insert(tied_max_degree.end(), max_degree_vertex);
        }
        else if(int(neighbors[vertex-1].size()) == max_degree){
            tied_max_degree.insert(tied_max_degree.end(), vertex);
        }
    }
    //if more than one vertex with max degree, choose one at random
    if(tied_max_degree.size() > 1 and opt!= NULL and opt->random_tie_breaks){
        max_degree_vertex = random_element(tied_max_degree);
    }

    //color first selected vertex
    ordering.push_back(max_degree_vertex);
    unselected_vertices.erase(max_degree_vertex);
    coloring.push_back({}); //new color class
    coloring[0].insert(max_degree_vertex); //assign first color to max degree vertex
    vertex_color[max_degree_vertex-1] = 0;

    //record for an uncolored vertex to how many colored vertices it is connected
    //also track degree of degree in uncolored subgraph
    std::vector<int> saturation_level(_ncount, 0);
    std::vector<int> uncolored_subgraph_degree(_ncount);
    for(Vertex v = 1; v <= _ncount; v++)
        uncolored_subgraph_degree[v-1] = int(neighbors[v-1].size());
    //don't want to select already selected vertex, remove from consideration
    saturation_level[max_degree_vertex - 1] = std::numeric_limits<int>::min();
    uncolored_subgraph_degree[max_degree_vertex - 1] = std::numeric_limits<int>::min();
    //update saturation level for neighbors of first selected vertex
    for(Vertex neighbor : neighbors[max_degree_vertex-1]){
        saturation_level[neighbor - 1]++;
        uncolored_subgraph_degree[neighbor - 1]--;
    }

    //main loop looking for the next vertex according to its saturation level
    while(not unselected_vertices.empty()){
        Vertex max_saturated_vertex = 0;
        int max_saturated = -1;
        //keep set of tied vertices to later choose one at random
        std::set<Vertex> saturated_tied_max_degrees;
        for(Vertex vertex : unselected_vertices){
            if(saturation_level[vertex - 1] > max_saturated){
                //found vertex with higher saturation, update
                max_saturated = saturation_level[vertex - 1];
                max_saturated_vertex = vertex;

                saturated_tied_max_degrees.clear();
                saturated_tied_max_degrees.insert(vertex);
            } else if(saturation_level[vertex - 1] == max_saturated){
                if(uncolored_subgraph_degree[vertex - 1] > uncolored_subgraph_degree[max_saturated_vertex - 1]) {
                    //found vertex with equal saturation but higher degree, update
                    max_saturated_vertex = vertex;

                    saturated_tied_max_degrees.clear();
                    saturated_tied_max_degrees.insert(vertex);
                }
                else if(neighbors[vertex - 1].size() == neighbors[max_saturated_vertex - 1].size()){
                    //same saturation level, same degree, add to set of possible candidates
                    saturated_tied_max_degrees.insert(vertex);
                }
            }
        }
        //select max saturated vertex, at random if there exist multiple
        if(saturated_tied_max_degrees.size() > 1 and opt!= NULL and opt->random_tie_breaks) {
            max_saturated_vertex = int(random_element(saturated_tied_max_degrees));
        }
        ordering.push_back(max_saturated_vertex);
        unselected_vertices.erase(max_saturated_vertex);

        //get all the colors of colored neighbors from the chosen vertex, ignore uncolored neighbors
        std::set<int> saturation_colors;
        saturation_colors.clear();
        for(Vertex neighbor : neighbors[max_saturated_vertex-1]){
            if(vertex_color[neighbor - 1] != -1)
                saturation_colors.insert(vertex_color[neighbor-1]);
        }

        //find lowest color not used by any neighbors
        int lowest_color = 0;
        for(int color = 0; color <= int(saturation_colors.size()); color++){
            if(saturation_colors.count(color) == 0){
                lowest_color = color;
                break;
            }
        }

        if(lowest_color < int(coloring.size())){
            //we can add chosen vertex to some existing color class
            coloring[lowest_color].insert(max_saturated_vertex);
            vertex_color[max_saturated_vertex-1] = lowest_color;
        } else if(lowest_color == int(coloring.size())){
            //new color class is created for max_saturated_vertex
            //before that we try to recolor two vertices to not have to use a new color
            bool recolored = try_color_swap(max_saturated_vertex, neighbors, coloring, vertex_color);
            if(not recolored){
                //not able to recolor, add new color class
                coloring.push_back({});
                coloring[lowest_color].insert(max_saturated_vertex);
                vertex_color[max_saturated_vertex-1] = lowest_color;
            }
        } else {
            std::cout << "Error: lowest color is too high, something must have gone wrong." << std::endl;
        }

        //new vertex has been selected so update the saturation levels of all its neighbors if it's a different color
        int msv_color = vertex_color[max_saturated_vertex - 1];
        for(Vertex vertex : neighbors[max_saturated_vertex - 1]){
            bool is_new_color = true;
            //possibly check that neighbor has been selected, i.e. it's level is equal to numeric_limits<int>::min()
            if(saturation_level[vertex - 1] != std::numeric_limits<int>::min()){
                for(Vertex adj_vertex : neighbors[vertex - 1]){
                    if (adj_vertex != max_saturated_vertex and vertex_color[adj_vertex - 1] == msv_color){
                        //found neighbor of neighbor that uses msv color so we don't increase the saturation level
                        is_new_color = false;
                        break;
                    }
                }
                if(is_new_color){
                    saturation_level[vertex - 1]++;
                }
                //also update degree of neighbor in uncolored subgraph
                uncolored_subgraph_degree[vertex -1]--;
            }
        }
        //remove selected vertex from any further consideration
        saturation_level[max_saturated_vertex - 1] = std::numeric_limits<int>::min();
        uncolored_subgraph_degree[max_saturated_vertex - 1] = std::numeric_limits<int>::min();
    }
    //colored all vertices
    return coloring;
}


Coloring Graph::max_connected_degree_coloring(Permutation &ordering) const {
    NeighborList neighbors = get_neighbor_list();
    return max_connected_degree_coloring(ordering, neighbors);
}

Coloring Graph::max_connected_degree_coloring(Permutation &ordering, const NeighborList& neighbors) const {
    ordering.clear();
    ordering.reserve(_ncount);
    Coloring coloring;
    //initialise all vertices to be uncolored
    std::vector<int> vertex_color(_ncount, -1);
    std::set<Vertex> unselected_vertices;
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++) {
        unselected_vertices.insert(unselected_vertices.begin(), vertex);
    }

    //find max degree vertex and select it as first colored vertex. ties are broken randomly
    Vertex max_degree_vertex = 0;
    std::set<Vertex> tied_max_degree;
    int max_degree = -1;
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++){
        if(int(neighbors[vertex-1].size()) > max_degree){
            max_degree = int(neighbors[vertex-1].size());
            max_degree_vertex = vertex;
            tied_max_degree.clear();
            tied_max_degree.insert(tied_max_degree.end(), max_degree_vertex);
        }
        else if(int(neighbors[vertex-1].size()) == max_degree){
            tied_max_degree.insert(tied_max_degree.end(), vertex);
        }
    }
    //if more than one vertex with max degree, choose one at random
    if(tied_max_degree.size() > 1 and opt!= NULL and opt->random_tie_breaks){
        max_degree_vertex = random_element(tied_max_degree);
    }

    //color first selected vertex
    ordering.push_back(max_degree_vertex);
    unselected_vertices.erase(max_degree_vertex);
    coloring.push_back({}); //new color class
    coloring[0].insert(max_degree_vertex); //assign first color to max degree vertex
    vertex_color[max_degree_vertex-1] = 0;

    //record for an uncolored vertex to how many colored vertices it is connected
    std::vector<int> num_selected_neighbors(_ncount, 0);
    //don't want to select already selected vertex
    num_selected_neighbors[max_degree_vertex - 1] = std::numeric_limits<int>::min();
    //update saturation level for neighbors of first selected vertex
    for(Vertex neighbor : neighbors[max_degree_vertex-1]){
        num_selected_neighbors[neighbor - 1]++;
    }

    //main loop looking for the next vertex according to its saturation level
    while(not unselected_vertices.empty()){
        Vertex max_saturated_vertex = 0;
        int max_saturated = -1;
        //keep set of tied vertices to later choose one at random
        std::set<Vertex> saturated_tied_max_degrees;
        for(Vertex vertex : unselected_vertices){
            if(num_selected_neighbors[vertex - 1] > max_saturated){
                //found vertex with higher saturation, update
                max_saturated = num_selected_neighbors[vertex - 1];
                max_saturated_vertex = vertex;

                saturated_tied_max_degrees.clear();
                saturated_tied_max_degrees.insert(vertex);
            } else if(num_selected_neighbors[vertex - 1] == max_saturated){
                if(neighbors[vertex - 1].size() > neighbors[max_saturated_vertex - 1].size()) {
                    //found vertex with equal saturation but higher degree, update
                    max_saturated_vertex = vertex;

                    saturated_tied_max_degrees.clear();
                    saturated_tied_max_degrees.insert(vertex);
                }
                else if(neighbors[vertex - 1].size() == neighbors[max_saturated_vertex - 1].size()){
                    //same saturation level, same degree, add to set of possible candidates
                    saturated_tied_max_degrees.insert(vertex);
                }
            }
        }
        //select max saturated vertex, at random if there exist multiple
        if(saturated_tied_max_degrees.size() > 1 and opt!= NULL and opt->random_tie_breaks) {
            max_saturated_vertex = int(random_element(saturated_tied_max_degrees));
        }
        ordering.push_back(max_saturated_vertex);
        unselected_vertices.erase(max_saturated_vertex);

        //get all the colors of colored neighbors from the chosen vertex, ignore uncolored neighbors
        std::set<int> saturation_colors;
        saturation_colors.clear();
        for(Vertex neighbor : neighbors[max_saturated_vertex-1]){
            if(vertex_color[neighbor - 1] != -1)
                saturation_colors.insert(vertex_color[neighbor-1]);
        }

        //find lowest color not used by any neighbors
        int lowest_color = 0;
        for(int color = 0; color <= int(saturation_colors.size()); color++){
            if(saturation_colors.count(color) == 0){
                lowest_color = color;
                break;
            }
        }

        if(lowest_color < int(coloring.size())){
            //we can add chosen vertex to some existing color class
            coloring[lowest_color].insert(max_saturated_vertex);
            vertex_color[max_saturated_vertex-1] = lowest_color;
        } else if(lowest_color == int(coloring.size())){
            //new color class is created for max_saturated_vertex
            //before that we try to recolor two vertices to not have to use a new color
            bool recolored = try_color_swap(max_saturated_vertex, neighbors, coloring, vertex_color);
            if(not recolored){
                //not able to recolor, add new color class
                coloring.push_back({});
                coloring[lowest_color].insert(max_saturated_vertex);
                vertex_color[max_saturated_vertex-1] = lowest_color;
            }
        } else {
            std::cout << "Error: lowest color is too high, something must have gone wrong." << std::endl;
        }

        //new vertex has been selected so update the the number of colored neighbors for that vertex's neighbors
        int msv_color = vertex_color[max_saturated_vertex-1];
        for(Vertex vertex : neighbors[max_saturated_vertex - 1]){
            //possibly check that neighbor has been selected, i.e. it's level is equal to numeric_limits<int>::min()
            if(num_selected_neighbors[vertex - 1] != std::numeric_limits<int>::min()){
                num_selected_neighbors[vertex - 1]++;
            }
        }
        //remove selected vertex from any further consideration
        num_selected_neighbors[max_saturated_vertex - 1] = std::numeric_limits<int>::min();
    }
    return coloring;
}


Permutation Graph::max_connected_degree_ordering() const {
    NeighborList neighbors = get_neighbor_list();
    return max_connected_degree_ordering(neighbors);
}


Permutation Graph::max_connected_degree_ordering(const NeighborList& neighbors) const {
    Permutation ordering;
    ordering.reserve(_ncount);
    std::set<Vertex> unselected_vertices;
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++) {
        unselected_vertices.insert(unselected_vertices.begin(), vertex);
    }

    Vertex max_degree_vertex = 0;
    int max_degree = -1;
    //find max degree vertex and select it first
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++){
        if(int(neighbors[vertex-1].size()) > max_degree){
            max_degree = int(neighbors[vertex-1].size());
            max_degree_vertex = vertex;
        }
    }
    ordering.push_back(max_degree_vertex);
    unselected_vertices.erase(max_degree_vertex);


    //record for an unselected vertex to how many selected vertices it is connected
    std::vector<int> num_connected(_ncount);
    //don't want to select already selected vertex
    num_connected[max_degree_vertex - 1] = std::numeric_limits<int>::min();
    //update connection count for neighbors of first selected
    for(Vertex neighbor : neighbors[max_degree_vertex-1]){
        num_connected[neighbor - 1]++;
    }

    while(not unselected_vertices.empty()){
        Vertex max_connections_vertex = 0;
        int max_connections = -1;
        for(Vertex vertex : unselected_vertices){
            if(num_connected[vertex - 1] > max_connections){
                max_connections = num_connected[vertex - 1];
                max_connections_vertex = vertex;
            } else if(num_connected[vertex - 1] == max_connections and (neighbors[vertex - 1].size() > neighbors[max_connections_vertex - 1].size())) {
                max_connections_vertex = vertex; //saturation is the same but vertex has a higher degree
            }
        }
        //select max conneceted vertex
        ordering.push_back(max_connections_vertex);
        unselected_vertices.erase(max_connections_vertex);
        for(Vertex neighbor : neighbors[max_connections_vertex - 1]){
            if(num_connected[neighbor - 1] != std::numeric_limits<int>::min()){
               num_connected[neighbor - 1]++;
            }
        }
        //remove selected vertex from any further consideration
        num_connected[max_connections_vertex - 1] = std::numeric_limits<int>::min();
    }
    return ordering;
}



Permutation Graph::max_degree_ordering() const {
    NeighborList neighbors = get_neighbor_list();
    return max_degree_ordering(neighbors);
}


Permutation Graph::max_degree_ordering(const NeighborList& neighbors) const {
    Permutation ordering;
    ordering.reserve(_ncount);
    for (Vertex i = 1; i <= _ncount; i++) {
        ordering.push_back(i);
    }
    std::sort(ordering.begin(), ordering.end(), [&neighbors](Vertex v, Vertex w){return neighbors[v-1].size() < neighbors[w-1].size();});
    std::reverse(ordering.begin(), ordering.end());

    return ordering;
}


void Graph::remove_last_vertices(int num_to_be_removed) {
    std::set<Vertex> remove_vertices;
    for(int i = 0; i < num_to_be_removed; i++){
        remove_vertices.insert(remove_vertices.begin(), int(_ncount) - i);
    }
    for(int edge_index = int(_ecount)-1; edge_index >= 0; edge_index--){
        if(remove_vertices.count(_elist[2 * edge_index]) or remove_vertices.count(_elist[2 * edge_index + 1])){
            _elist.erase(std::next(_elist.begin(), 2*edge_index+1));
            _elist.erase(std::next(_elist.begin(), 2*edge_index));
            _ecount--;
        }
    }
    _ncount -= num_to_be_removed;
}


void Graph::peel_graph_ordered(int peeling_degree) {
    NeighborList neighbors = get_neighbor_list();
    peel_graph_ordered(peeling_degree, neighbors);
}

void Graph::peel_graph_ordered(int peeling_degree, const NeighborList& neighbors) {
    std::set<Vertex> small_degree_vertices;
    for(Vertex i = 1; i <= _ncount; i++){
        if(neighbors[i-1].size() < peeling_degree){
            small_degree_vertices.insert(small_degree_vertices.end(), i);
        }
    }
    remove_vertices(small_degree_vertices);
    std::cout << "peeled: " << ncount() << ", " << ecount() << std::endl;
}

void Graph::remove_vertex(Vertex remove) {
    for(int edge_index = int(_ecount)-1; edge_index >= 0; edge_index--){
        if(remove == _elist[2 * edge_index] or remove ==_elist[2 * edge_index + 1]){
            _elist.erase(std::next(_elist.begin(), 2*edge_index+1));
            _elist.erase(std::next(_elist.begin(), 2*edge_index));
            _ecount--;
        }
        else{
            if( _elist[2 * edge_index + 1] > remove) _elist[2 * edge_index + 1]--;
            if( _elist[2 * edge_index] > remove) _elist[2 * edge_index]--;
        }
        bool a = false;
        if(_elist[2 * edge_index] <= 1 or _elist[2 * edge_index + 1] <= 1){
            a = true;
        }
    }
    _ncount--;
}


void Graph::remove_vertices(const std::set<Vertex>& to_remove) {
    //order of removing vertices is important here, be removing later vertices first, we keep the numbering intact
    for( auto v_it = to_remove.rbegin(); v_it != to_remove.rend(); v_it++){
        remove_vertex(*v_it);
    }
}


void Graph::remove_dominated_vertices() {
    NeighborList neighbors = get_neighbor_list();
    remove_dominated_vertices(neighbors);
}


void Graph::remove_dominated_vertices(const NeighborList& neighbors) {

    std::set<Vertex> dominated_vertices;
    //TDOD only check for dominated if one set of neighbors is larger than the other
    for(Vertex v = 1; v <= _ncount; v++){
        for(Vertex w = v+1; w <= _ncount; w++){
            if(std::includes(neighbors[v-1].begin(), neighbors[v-1].end(),
                             neighbors[w-1].begin(), neighbors[w-1].end())){
                dominated_vertices.insert(dominated_vertices.end(), w);
            }else if(std::includes(neighbors[w-1].begin(), neighbors[w-1].end(),
                               neighbors[v-1].begin(), neighbors[v-1].end())){
                dominated_vertices.insert(dominated_vertices.end(), v);
                break;
            }
        }
    }
    remove_vertices(dominated_vertices);
    std::cout << "dominated: " << ncount() << ", " << ecount() << std::endl;
}


void Graph::vertex_fusion(const std::set<Vertex>& clique) {
   NeighborList neighbors = get_neighbor_list();
    vertex_fusion(clique, neighbors);
    ;
}


void Graph::vertex_fusion(const std::set<Vertex>& clique, const NeighborList& neighbors) {
    for(Vertex v = 1; v <= _ncount; v++) {
        if(clique.count(v))
            continue;
        for(Vertex w : clique){
            if(neighbors[v-1].count(w))
                continue;

            bool adjacent_to_all = true;
            for(Vertex u : clique) {
                if (w == u)
                    continue;

            }
        }
    }
}

void Graph::edge_addition(const std::set<Vertex>& clique) {
    NeighborList neighbors = get_neighbor_list();
    edge_addition(clique, neighbors);
}

void Graph::edge_addition(const std::set<Vertex>& clique, const NeighborList& neighbors) {
    for(Vertex v = 1; v <= _ncount; v++){
        if(clique.count(v))
            continue;
        for(Vertex w = v+1; w <= _ncount; w++){
            if(clique.count(w) or neighbors[v-1].count(w))
                continue;

            bool edge_always_exists = true;
            for(Vertex u : clique){
                if(v == u or w == u)
                    continue;

                if((not neighbors[v-1].count(u)) and (not neighbors[w-1].count(u))){
                    edge_always_exists = false;
                    break;
                }
            }
            if(edge_always_exists){
                //add edge {v,w}
                _elist.push_back(v);
                _elist.push_back(w);
                _ecount++;
            }
        }
    }
    std::cout << "edge addition : " << ncount() << ", " << ecount() << std::endl;
}


std::set<Vertex> Graph::find_clique(int nrbranches) const{
    //find size of some (optimally maximal) clique
    COLORset* cliques = (COLORset *)NULL;
    int ncliques = 0;
    int pval = 0;
    if(nrbranches == -1)
        nrbranches = int(ncount())/10 + 10;
    std::vector<int> one_weights(ncount(), 1);
    int rval = COLORclique_ostergard(&cliques, &ncliques, int(ncount()), int(ecount()), shifted_elist().data(),
                                     one_weights.data(), COLOR_MAXINT, &pval, nrbranches);
    if(rval) std::cout <<"Failed in COLORclique_ostergard." << std::endl;

    std::set<Vertex> clique;
    std::transform(cliques[0].members, cliques[0].members+cliques[0].count, std::inserter(clique, clique.end()), [](int v){return (Vertex)v+1;});

    std::cout << "Found clique of size " << pval << std::endl;
    return clique;
}




bool try_color_swap(Vertex max_saturated_vertex, const NeighborList &neighbors, Coloring &coloring,
                    std::vector<int> &vertex_color) {
    for(unsigned int j = 0; j < coloring.size(); j++){
        for(unsigned int k = j+1; k < coloring.size(); k++){
            std::set<Vertex> neighbors_v_intersect_color_j;
            std::set_intersection(neighbors[max_saturated_vertex-1].begin(), neighbors[max_saturated_vertex-1].end(),
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
                    coloring[j].insert(max_saturated_vertex);
                    coloring[k].insert(u);
                    vertex_color[max_saturated_vertex-1] = int(j);
                    vertex_color[u-1] = int(k);
                    k = coloring.size();
                    j = coloring.size();
                    return true;//was able to swap colors
                }
            }
        }
    }
    return false; //unable to swap two colors
}













