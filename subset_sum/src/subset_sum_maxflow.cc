#include "subset_sum_maxflow.h"
#include <cfloat>

SubsetSumMaxflow::SubsetSumMaxflow(const SubsetSum &ss)
	:sss(ss)
{
	alpha = sss.k;
}

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
		PEB p = add_edge(0, i + 1, gr);
		assert(p.second == true);
		e2w.insert(PED(p.first, alpha));
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
			PEB p2 = add_edge(k + 1, n, gr);
			assert(p1.second == true);
			assert(p2.second == true);
			e2w.insert(PED(p1.first, DBL_MAX));
			e2w.insert(PED(p2.first, DBL_MAX));
		}
	}

	// last element: sink
	int n = num_vertices(gr);
	add_vertex(gr);
	for(int i = sss.sum_array.size() + 1; i < n; i++)
	{
		PEB p = add_edge(i, n, gr);
		assert(p.second == true);
		e2w.insert(PED(p.first, 1.0));
	}

	return 0;
}
