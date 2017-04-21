#include "subset_sum_maxflow.h"

#include <boost/graph/push_relabel_max_flow.hpp>
#include <cfloat>

SubsetSumMaxflow::SubsetSumMaxflow(const SubsetSum &ss)
	:sss(ss)
{
	alpha = sss.k;
}

int SubsetSumMaxflow::solve()
{
	build_graph();
	compute_maxflow();
	return 0;
}

int SubsetSumMaxflow::build_graph()
{
	edge_capacity_map e2w = get(edge_capacity, gr);
	edge_reverse_map rev = get(edge_reverse, gr);

	// 0: source 
	gr.clear();
	add_vertex(gr);

	// [1, sum_array.size() + 1) -> sss.vertices
	for(int i = 0; i < sss.sum_array.size(); i++) 
	{
		add_vertex(gr);
		PEB p1 = add_edge(0, i + 1, gr);
		PEB p2 = add_edge(i + 1, 0, gr);
		e2w[p1.first] = alpha;
		e2w[p2.first] = 0;
		rev[p1.first] = p2.first;
		rev[p2.first] = p1.first;
	}

	// [sum_array.size() + 1, X) -> sss.edges
	for(int i = 0; i < sss.sum_array.size(); i++)
	{
		for(int j = 0; j < sss.sum_array[i].size(); j++)
		{
			int k = sss.sum_array[i][j].next_index;
			int n = num_vertices(gr);
			add_vertex(gr);
			PEB p1 = add_edge(i + 1, n, gr);
			PEB p2 = add_edge(n, i + 1, gr);
			PEB p3 = add_edge(k + 1, n, gr);
			PEB p4 = add_edge(n, k + 1, gr);
			e2w[p1.first] = alpha * sss.sum_array.size() + 1;
			e2w[p3.first] = alpha * sss.sum_array.size() + 1;
			e2w[p2.first] = 0;
			e2w[p4.first] = 0;
			rev[p1.first] = p2.first;
			rev[p2.first] = p1.first;
			rev[p3.first] = p4.first;
			rev[p4.first] = p3.first;
		}
	}

	// last element: sink
	int n = num_vertices(gr);
	add_vertex(gr);
	for(int i = sss.sum_array.size() + 1; i < n; i++)
	{
		PEB p1 = add_edge(i, n, gr);
		PEB p2 = add_edge(n, i, gr);
		e2w[p1.first] = 1;
		e2w[p2.first] = 0;
		rev[p1.first] = p2.first;
		rev[p2.first] = p1.first;
	}

	return 0;
}

int SubsetSumMaxflow::build_test_graph()
{
	edge_capacity_map e2w = get(edge_capacity, gr);
	edge_reverse_map rev = get(edge_reverse, gr);

	// 0: source 
	gr.clear();
	add_vertex(gr);

	for(int i = 0; i < 3; i++) 
	{
		add_vertex(gr);
		PEB p1 = add_edge(0, i + 1, gr);
		PEB p2 = add_edge(i + 1, 0, gr);
		e2w[p1.first] = 4;
		e2w[p2.first] = 0;
		rev[p1.first] = p2.first;
		rev[p2.first] = p1.first;
	}

	for(int j = 0; j < 2; j++)
	{
		for(int i = 0; i < 3; i++)
		{
			int n = num_vertices(gr);
			add_vertex(gr);
			PEB p1 = add_edge(i + 1, n, gr);
			PEB p2 = add_edge(n, i + 1, gr);
			e2w[p1.first] = 100;
			e2w[p2.first] = 0;
			rev[p1.first] = p2.first;
			rev[p2.first] = p1.first;
		}
	}

	// last element: sink
	int n = num_vertices(gr);
	add_vertex(gr);
	for(int i = 4; i <= 5; i++)
	{
		PEB p1 = add_edge(i, n, gr);
		PEB p2 = add_edge(n, i, gr);
		e2w[p1.first] = 1;
		e2w[p2.first] = 0;
		rev[p1.first] = p2.first;
		rev[p2.first] = p1.first;
	}

	return 0;
}

int SubsetSumMaxflow::compute_maxflow()
{
	double flow = push_relabel_max_flow(gr, 0, num_vertices(gr) - 1);
	printf("max-flow = %.3lf\n", flow);
	return 0;
}

int SubsetSumMaxflow::print()
{
	printf("alpha = %.2lf\n", alpha);
	printf("number of input spectrum = %lu\n", sss.spectrum.size());
	printf("number of vertices in subsetsum = %lu\n", sss.sum_array.size());
	printf("number of edges in subsetsum = %lu\n", num_vertices(gr) - 2 - sss.sum_array.size());
	printf("number of edges in network = %lu\n", num_edges(gr));
	//sss.analyze_degree_distribution();
	return 0;
}

