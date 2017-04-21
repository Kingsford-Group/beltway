#include "subset_sum_maxflow.h"

SubsetSumMaxflow::SubsetSumMaxflow(const SubsetSum &ss)
	:sss(ss)
{}

int SubsetSumMaxflow::solve()
{
	return 0;
}

int SubsetSumMaxflow::print()
{
	return 0;
}

int SubsetSumMaxflow::build_graph()
{
	// 0: source 
	gr.clear();
	add_vertex(gr);

	// [1, sum_array.size() + 1) -> sss.vertices
	for(int i = 0; i < sss.sum_array.size(); i++) 
	{
		add_vertex(gr);
		add_edge(0, i + 1, gr);
	}

	// [sum_array.size() + 1, X) -> sss.edges
	for(int i = 0; i < sss.sum_array.size(); i++)
	{
		for(int j = 0; j < sss.sum_array[i].size(); j++)
		{
			int k = sss.sum_array[i][j].next_index;
			int n = num_vertices(gr);
			add_vertex(gr);
			add_edge(i + 1, n, gr);
			add_edge(k + 1, n, gr);
		}
	}

	// last element: sink
	int n = num_vertices(gr);
	add_vertex(gr);
	for(int i = sss.sum_array.size() + 1; i < n; i++)
	{
		add_edge(i, n, gr);
	}

	return 0;
}
