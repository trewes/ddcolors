#ifndef DDCOLORS_GRAPH_H
#define DDCOLORS_GRAPH_H

/*
 * Graph.h
 * Purpose: Implementation of a class to maintain the data of a graph and calling several functions on a graph.
 * Includes reading in and building a graph as well functions to get a permuted graph, several different
 * coloring heuristics or orderings, finding a clique and functions to remove (certain) vertices from the graph.
 */

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
#include "color.h"  //needed for finding a clique, using the COLORclique_ostergard method from exactcolors
}

typedef unsigned int Vertex;
//type to store for a vertex j the set of neighboring vertices as neighbors[j]
typedef std::vector<std::set<Vertex> > NeighborList;
//note that permutations are given as positive integers, there is no i s.t. perm[i] = 0, it's always positive
typedef std::vector<unsigned int> Permutation;
typedef unsigned int Color;
typedef std::set<Color> ColorClass;
typedef std::vector<ColorClass> Coloring;

//allows to print sets and vectors in a nice format via std::cout
template<typename T>
std::ostream &operator<<(std::ostream &s, const std::set<T> &set);
template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec);

//helper function to chose a random element from a given set of elements
template<typename T>
T random_element(std::set<T> const &v);

//returns the identity permutation, i.e. a vector with vec[i] = i+1 of the given size
Permutation identity(size_t size);

//returns the inverse permutation s.t. inverse[perm[i]] = i+1
Permutation perm_inverse(const Permutation &perm);



/*
 * Graph
 * Purpose: store the basic data of a graph, i.e. node count and edges and allowing to perform functions on a graph
 */
class Graph{
public:

    /*
     * Member variables
     * _ncount : the number of vertices
     * _ecount : the number of edges
     * _elist : the edges of the graph, for i = 0, ..., _ecount we have that {_elist[2i], _elist[2i+1]} is an edge
     * random_tiebreaks : if set to true, allows for randomness during choices of the coloring heuristics/orderings
     *
     * Graph : builds a graph from given vertex/edge counts and an edge list
     *         or  from a given file in either dimacs or Brendan McKay's sparse graph format ".g6" (only the first line)
     *         or  reads in a string containing a graph in Brendan McKay's sparse graph format ".g6"
     * get_neighbor_list : Puts the edge list in a more useful format, stores to a vertex
     *                     the set of adjacent vertices as neighbors[j]
     * ncount, ecount, elist : return the corresponding data field of the class
     * shifted_elist : a helper function to ease communication with the clique algorithm from exactcolors,
     *                 shifts the vertices to be labeled beginning from 0 instead of 1
     * perm_graph : returns the permuted graph for a given permutation
     * use_random_tiebreaks : enables the use of randomness during choices in the coloring heuristics and orderings
     * print : outputs the graph in an adjacency list format
     *
     */

    Graph(unsigned int ncount, unsigned int ecount, std::vector<Vertex> elist);
    Graph(const char *filename);
    Graph(std::string g6_string);

    NeighborList get_neighbor_list() const;

    unsigned int ncount() const;
    unsigned int ecount() const;
    std::vector<Vertex> elist() const;

    std::vector<int> shifted_elist() const;

    Graph perm_graph(const Permutation &perm) const;

    void use_random_tiebreaks() { random_tiebreaks = true; }

    void print() const;

private:
    unsigned int _ncount;
    unsigned int _ecount;
    std::vector<Vertex> _elist;
    bool random_tiebreaks = false;

    //helper functions to build the graph, either reading in dimacs or .g6 format
    void read_dimacs(const char *filename);
    void read_graph6(const char *filename);
    void read_graph6_string(std::string g6_string);

public:

    /*
     * Several coloring heuristics are implemented:
     *
     * dsatur: select the next vertex as the one with the highest saturation degree, ties
     *         are broken by the highest degree and randomly after that
     * dsatur_original: dsatur how it was originally introduced by Brelaz, select the next vertex
     *                  as the one with the highest saturation degree, ties are broken by the
     *                  highest degree in the uncolored subgraph and randomly after that
     * max_connected_degree: select the next vertex as the one with the highest amount of colored neighbors,
     *                       not necessarily colored differently as for dsatur.
     *                       Ties are broken by the highest degree and randomly after that.
     *                       _coloring yields a coloring, _ordering only the vertex ordering
     * constraint_graph_ordering : a ordering only algorithm, computes the ordering minimizing the
     *                             number of conflicts in the conflict graph
     *
     * Additionally, a clique can be passed to each algorithm whose vertices colors will be fixed. This can lead
     * to fewer mistakes in the heuristic since we already know these vertices will have to be colored distinctly.
     *
     * try_color_swap : implements the recolor technique for the coloring heuristics used.
     *                  if a new color class would be created, check if we can swap the colors of two vertices
     *                  to avoid having to use a new color class. This generally improves the heuristics
     *
     * Also defines OrderType, according to the algorithm by which the ordering was obtained
     */
    enum OrderType {Lexicographic, Dsatur, DsaturOriginal, MaxConnectedDegree, MinWidth};


    Coloring dsatur(Permutation &ordering, const std::set<Vertex> &clique = {}) const;
    Coloring dsatur(Permutation &ordering, const NeighborList &neighbors, const std::set<Vertex> &clique = {}) const;

    Coloring dsatur_original(Permutation &ordering, const std::set<Vertex> &clique = {}) const;
    Coloring
    dsatur_original(Permutation &ordering, const NeighborList &neighbors, const std::set<Vertex> &clique = {}) const;


    Coloring max_connected_degree_coloring(Permutation &ordering, const std::set<Vertex> &clique = {}) const;
    Coloring max_connected_degree_coloring(Permutation &ordering, const NeighborList &neighbors,
                                           const std::set<Vertex> &clique = {}) const;

    Permutation max_connected_degree_ordering(const std::set<Vertex> &clique = {}) const;
    Permutation max_connected_degree_ordering(const NeighborList &neighbors, const std::set<Vertex> &clique = {}) const;

    int constraint_graph_width();
    Permutation constraint_graph_ordering();

    static bool try_color_swap(Vertex max_saturated_vertex, const NeighborList &neighbors, Coloring &coloring,
                        std::vector<int> &vertex_color);

    //this is never used and not useful as a vertex ordering for this problem
    Permutation max_degree_ordering() const;
    Permutation max_degree_ordering(const NeighborList &neighbors) const;

    /*
     *Several functions to remove (certain) vertices from the graph.
     * note that these function are done on the graph itself and not a new one
     *
     * remove_vertex : removes a single vertex from the graph and relabels the vertices to
     *                 be labeled from 1 to n-1
     * remove_vertices : removes the given set of vertices from the graph one by one, taking care
     *                   that the right vertex is removed even after the labels change
     * remove_vertices_together : same as previous but removes the vertices all at once
     *
     * peel_graph : removes all vertices with degree strictly smaller than the specified number
     * remove_dominated_vertices : removes all dominated vertices from the graph. A vertex is dominated by another
     *                             vertex if all his neighbors are also adjacent to that other vertex
     *
     * find_clique : uses COLORclique_ostergard from exactcolors to find a clique, the argument specifies the
     *               number of branches to explore before stopping (if large enough, finds a largest clique)
     */

    void remove_vertex(Vertex remove);
    void remove_vertices(const std::set<Vertex> &to_remove);
    void remove_vertices_together(const std::set<Vertex> &to_remove);

    void peel_graph(int peeling_degree);
    void peel_graph(int peeling_degree, const NeighborList &neighbors);

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

};

#endif //DDCOLORS_GRAPH_H
