#ifndef DDCOLORS_GRAPH_H
#define DDCOLORS_GRAPH_H

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

typedef unsigned int Vertex;
typedef std::vector<std::set<Vertex> > NeighborList;
typedef std::vector<unsigned int> Permutation;
typedef std::set<unsigned int> ColorClass;
typedef std::vector<ColorClass> Coloring;


template<typename T> std::ostream& operator<<(std::ostream& s, const std::set<T> & set);
template<typename T> std::ostream& operator<<(std::ostream& s, const std::vector<T> & vec);

Permutation perm_inverse(const Permutation& perm);

class Graph{
public:
    Graph(unsigned int ncount, unsigned int ecount, std::vector<Vertex> elist);
    Graph(const char* filename);
    Graph(std::string g6_string);
    NeighborList get_neighbor_list() const;
    unsigned int ncount() const;
    unsigned int ecount() const;
private:
    unsigned int _ncount;
    unsigned int _ecount;
    std::vector<Vertex> _elist;
    void read_dimacs(const char *filename);
    void read_graph6(const char *filename);
    //this function of reading in the graph6 format was taken from treedecomposition.com
    void read_graph6_string(std::string g6_string);
public:
    enum OrderType {Lexicographic, Dsatur, Max_Connected_degree};
    Graph perm_graph(const Permutation &perm) const;
    Coloring dsatur(Permutation &ordering) const;
    Coloring dsatur(Permutation &ordering, const NeighborList& neighbors) const;
    Permutation max_connected_degree_ordering() const;
    Permutation max_connected_degree_ordering(const NeighborList& neighbors) const;

    Graph peel_graph(int peeling_degree) const;
    Graph peel_graph(int peeling_degree, const NeighborList& neighbors) const;

    void remove_last_vertices(int num_to_be_removed);

    void print() const;
};

#endif //DDCOLORS_GRAPH_H
