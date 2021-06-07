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


Permutation perm_inverse(const Permutation& perm) {
    Permutation inverse(perm.size());
    for(size_t i=0, max = perm.size(); i<max; i++){
        inverse[perm[i]-1] = i+1;
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
    std::vector<int> vertex_color(_ncount, -1);
    std::set<Vertex> unselected_vertices;

    unsigned int max_degree_vertex = 0;
    int max_degree = 0;
    //find max degree vertex and select it first
    for(size_t vertex = 1; vertex <= _ncount; vertex++){
        if(int(neighbors[vertex-1].size()) > max_degree){
            max_degree = int(neighbors[vertex-1].size());
            max_degree_vertex = vertex;
        }
        unselected_vertices.insert(unselected_vertices.begin(), vertex);
    }
    ordering.push_back(max_degree_vertex);
    unselected_vertices.erase(max_degree_vertex);
    coloring.push_back({}); //new color class
    coloring[0].insert(max_degree_vertex); //assign first color to max degree vertex
    vertex_color[max_degree_vertex-1] = 0;

    //record for an unselected vertex to how many selected vertices it is connected
    std::vector<int> saturation_level(_ncount);
    //don't want to select already selected vertex
    saturation_level[max_degree_vertex - 1] = std::numeric_limits<int>::min();
    //update saturation level for neighbors of first selected
    for(Vertex neighbor : neighbors[max_degree_vertex-1]){
        saturation_level[neighbor - 1]++;
    }

    while(not unselected_vertices.empty()){
        int max_saturated_vertex = -1;
        int max_saturated = -1;
        std::set<int> saturation_colors;
        for(Vertex vertex : unselected_vertices){
            if(saturation_level[vertex - 1] > max_saturated){
                max_saturated = saturation_level[vertex - 1];
                max_saturated_vertex = vertex;

                saturation_colors.clear();
                for(Vertex neighbor : neighbors[vertex-1]){
                    saturation_colors.insert(vertex_color[neighbor-1]);
                }
            } else if(saturation_level[vertex - 1] == max_saturated
                    and (neighbors[vertex - 1].size() > neighbors[max_saturated_vertex - 1].size())) { //TODO use degree in uncolored subgraph
                max_saturated_vertex = vertex; //saturation is the same but vertex has a higher degree

                saturation_colors.clear();
                for(Vertex neighbor : neighbors[vertex-1]){
                    saturation_colors.insert(vertex_color[neighbor-1]);
                }
            }
        }
        //select max saturated vertex
        ordering.push_back(max_saturated_vertex);
        unselected_vertices.erase(max_saturated_vertex);

        //find lowest color not used by any neighbors
        int lowest_color = 0;
        for(int color = 0; color <= int(saturation_colors.size()); color++){
            if(saturation_colors.count(color) == 0){
                lowest_color = color;
                break;
            }
        }

        vertex_color[max_saturated_vertex-1] = lowest_color;
        if(lowest_color < int(coloring.size())){
            coloring[lowest_color].insert(max_saturated_vertex);
        } else if(lowest_color == int(coloring.size())){
            coloring.push_back({});
            coloring[lowest_color].insert(max_saturated_vertex);
        } else {
            std::cout << "Error: lowest color is too high, something must have gone wrong." << std::endl;
        }

        //new vertex has been selected so update the saturation levels of all nodes
        for(Vertex neighbor : neighbors[max_saturated_vertex-1]){
            //possibly check that neighbor has been selected, i.e. it's level is equal to numeric_limits<int>::min()
            saturation_level[neighbor-1]++;
        }
        saturation_level[max_saturated_vertex-1] = std::numeric_limits<int>::min();
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

    int max_degree_vertex = -1;
    int max_degree = -1;
    //find max degree vertex and select it first
    for(size_t vertex = 1; vertex <= int(_ncount); vertex++){
        if(int(neighbors[vertex-1].size()) > max_degree){
            max_degree = int(neighbors[vertex-1].size());
            max_degree_vertex = vertex;
        }
        unselected_vertices.insert(unselected_vertices.begin(), vertex);
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
        unsigned int max_connections_vertex = 0;
        unsigned int max_connections = 0;
        for(Vertex vertex : unselected_vertices){
            if(num_connected[vertex - 1] > max_connections){
                max_connections = num_connected[vertex - 1];
                max_connections_vertex = vertex;
            } else if(num_connected[vertex - 1] == max_connections and (neighbors[vertex - 1].size() > neighbors[max_connections_vertex - 1].size())) {
                max_connections_vertex = vertex; //saturation is the same but vertex has a higher degree
            }
        }
        //select max saturated vertex
        ordering.push_back(max_connections_vertex);
        unselected_vertices.erase(max_connections_vertex);
        num_connected[max_connections_vertex - 1] = std::numeric_limits<int>::min();
        for(Vertex neighbor : neighbors[max_connections_vertex - 1]){
            num_connected[neighbor - 1]++;
        }
    }
    return ordering;
}

Graph Graph::peel_graph(int peeling_degree) const {
    NeighborList neighbors = get_neighbor_list();
    return peel_graph(peeling_degree, neighbors);
}

Graph Graph::peel_graph(int peeling_degree, const NeighborList& neighbors) const {

    std::set<Vertex> small_degree_vertices;
    for(Vertex i = 1; i <= _ncount; i++){
        if(neighbors[i-1].size() < peeling_degree){
            small_degree_vertices.insert(small_degree_vertices.end(), i);
        }
    }

    Permutation perm(_ncount);
    for (int i = 1; i <= _ncount; ++i){
        perm[i-1] = i;
    }
    int index = 0;
    for(auto remove_it = small_degree_vertices.rbegin(); remove_it != small_degree_vertices.rend(); remove_it++, index++){
        std::swap(perm[*remove_it - 1], perm[_ncount-1 - index]);
    }
    Graph permutated_graph = perm_graph(perm_inverse(perm));
    permutated_graph.remove_last_vertices(int(small_degree_vertices.size()));
    return permutated_graph;
}


void Graph::remove_last_vertices(int num_to_be_removed) {
    std::set<Vertex> remove_vertices;
    for(int i = 0; i < num_to_be_removed; i++){
        remove_vertices.insert(remove_vertices.begin(), int(_ncount) - i);
    }
    NeighborList neighbors = get_neighbor_list();
    for(int edge_index = int(_ecount)-1; edge_index >= 0; edge_index--){
        if(remove_vertices.count(_elist[2 * edge_index]) or remove_vertices.count(_elist[2 * edge_index + 1])){
            _elist.erase(std::next(_elist.begin(), 2*edge_index+1));
            _elist.erase(std::next(_elist.begin(), 2*edge_index));
            _ecount--;
        }
    }
    _ncount -= num_to_be_removed;
}

void Graph::print() const {
    std::cout << "Graph has " << _ncount << " vertices and " << _ecount << " edges." << std::endl;
    NeighborList neighbors = get_neighbor_list();
    for(std::set<Vertex>& adj : neighbors){
            std::cout << adj << std::endl;
    }
}



















