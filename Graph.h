#ifndef DDCOLORS_GRAPH_H
#define DDCOLORS_GRAPH_H

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>

extern "C" {
#include "color.h"
}

typedef unsigned int Vertex;
typedef std::vector<std::set<Vertex> > NeighborList;
typedef std::vector<unsigned int> Permutation;
typedef unsigned int Color;
typedef std::set<Color> ColorClass;
typedef std::vector<ColorClass> Coloring;


template<typename T>
std::ostream &operator<<(std::ostream &s, const std::set<T> &set);

template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec);

template<typename T>
T random_element(std::set<T> const &v);

Permutation identity(size_t size);

Permutation perm_inverse(const Permutation &perm);


class Graph{
public:
    Graph(unsigned int ncount, unsigned int ecount, std::vector<Vertex> elist);

    Graph(const char *filename);

    Graph(std::string g6_string);//only used this for testing

    ~Graph() = default;

    //get information about the graph
    NeighborList get_neighbor_list() const;

    unsigned int ncount() const;

    unsigned int ecount() const;

    std::vector<Vertex> elist() const;

    std::vector<int> shifted_elist() const; //used as helper function for clique algorithm from exactcolors
    void print() const;

    void use_random_tiebreaks() { random_tiebreaks = true; }

private:
    //data fields to store the graph
    unsigned int _ncount;
    unsigned int _ecount;
    std::vector<Vertex> _elist;

    //helper functions to build the graph
    void read_dimacs(const char *filename);

    void read_graph6(const char *filename);

    //this function of reading in the graph6 format was taken from treedecomposition.com
    void read_graph6_string(std::string g6_string);

    //setting use randomness during the coloring heuristic when breaking a tie
    bool random_tiebreaks = false;
public:
    enum OrderType {Lexicographic, Dsatur, Dsatur_original, Max_Connected_degree, MinWidth};

    Graph perm_graph(const Permutation &perm) const;

    //a clique can be passed to the coloring heuristic, oftentimes yielding better bound and ordering
    //coloring of the clique will definitely be optimal in the sense they can't but all get different colors

    //dsatur with max saturated degrees but first ties are broken by their degree in the whole graph
    Coloring dsatur(Permutation &ordering, const std::set<Vertex> &clique = {}) const;

    Coloring dsatur(Permutation &ordering, const NeighborList &neighbors, const std::set<Vertex> &clique = {}) const;

    //dsatur with max saturated degrees but first ties are broken by their degree in the uncolored subgraph
    Coloring dsatur_original(Permutation &ordering, const std::set<Vertex> &clique = {}) const;

    Coloring
    dsatur_original(Permutation &ordering, const NeighborList &neighbors, const std::set<Vertex> &clique = {}) const;

    //max connected degree, either only gets that ordering or also gets a greedy coloring using that heuristic/ordering
    Coloring max_connected_degree_coloring(Permutation &ordering, const std::set<Vertex> &clique = {}) const;

    Coloring max_connected_degree_coloring(Permutation &ordering, const NeighborList &neighbors,
                                           const std::set<Vertex> &clique = {}) const;

    Permutation max_connected_degree_ordering(const std::set<Vertex> &clique = {}) const;

    Permutation max_connected_degree_ordering(const NeighborList &neighbors, const std::set<Vertex> &clique = {}) const;

    //simply sort the vertices by descending degrees
    Permutation max_degree_ordering() const;

    Permutation max_degree_ordering(const NeighborList &neighbors) const;

    //removes single or multiple given vertices from the graph, relabels remaining vertices accordingly
    void remove_vertex(Vertex remove);

    void remove_vertices(const std::set<Vertex> &to_remove);

    void remove_vertices_together(const std::set<Vertex> &to_remove);

    //removes the k last vertices, no relabeling of vertices is necessary
    void remove_last_vertices(int num_to_be_removed);

    //removes all vertices with degree strictly less that the peeling_degree
    void peel_graph_ordered(int peeling_degree);

    void peel_graph_ordered(int peeling_degree, const NeighborList &neighbors);

    //removes all dominated vertices from the graph
    void remove_dominated_vertices();

    void remove_dominated_vertices(const NeighborList &neighbors);

    std::set<Vertex> find_clique(int nrbranches = -1) const;


    //only applicable for deciding k-colorability for given k = size of some clique
    //can not help at all if the chromatic number is greater than the largest clique
    //TODO remove for hand-in
    void vertex_fusion(const std::set<Vertex> &clique);

    void vertex_fusion(const std::set<Vertex> &clique, const NeighborList &neighbors);

    void edge_addition(const std::set<Vertex> &clique);

    void edge_addition(const std::set<Vertex> &clique, const NeighborList &neighbors);


    int constraint_graph_width();

    Permutation constraint_graph_ordering();

};

//helper functions to possibly swap to colors and improve coloring found through heuristics
//also adjusts the vertex_color vector if one is given
bool try_color_swap(Vertex max_saturated_vertex, const NeighborList &neighbors, Coloring &coloring,
                    std::vector<int> &vertex_color);


#endif //DDCOLORS_GRAPH_H
