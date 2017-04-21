#ifndef __MYGRAPH__
#define __MYGRAPH__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

using namespace boost;

typedef adjacency_list_traits<vecS, vecS, directedS> Traits;

typedef adjacency_list<setS, vecS, directedS, 
		property<vertex_index_t, int>,
		property<edge_capacity_t, double,
		property<edge_residual_capacity_t, double,
		property<edge_reverse_t, Traits::edge_descriptor> > > > DiGraph;

typedef graph_traits<DiGraph>::vertex_iterator vertex_iterator;
typedef graph_traits<DiGraph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<DiGraph>::edge_iterator edge_iterator;
typedef graph_traits<DiGraph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<DiGraph>::edge_descriptor edge_descriptor;

typedef property_map<DiGraph, edge_capacity_t>::type edge_capacity_map;
typedef property_map<DiGraph, edge_reverse_t>::type edge_reverse_map;

static vertex_descriptor VNULL = graph_traits<DiGraph>::null_vertex();

#endif
