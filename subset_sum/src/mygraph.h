#ifndef __MYGRAPH__
#define __MYGRAPH__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

using namespace boost;

typedef adjacency_list<setS, vecS, directedS> DiGraph;

typedef graph_traits<DiGraph>::vertex_iterator vertex_iterator;
typedef graph_traits<DiGraph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<DiGraph>::edge_iterator edge_iterator;
typedef graph_traits<DiGraph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<DiGraph>::edge_descriptor edge_descriptor;

static vertex_descriptor VNULL = graph_traits<DiGraph>::null_vertex();

#endif
