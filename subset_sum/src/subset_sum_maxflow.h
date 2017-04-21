#ifndef __SUBSET_SUM_MAXFLOW_H__
#define __SUBSET_SUM_MAXFLOW_H__

#include "subset_sum.h"
#include "mygraph.h"
#include <map>

using namespace std;

typedef pair<edge_descriptor, bool> PEB;
typedef pair<edge_descriptor, double> PED;
typedef map<edge_descriptor, double> MED;


class SubsetSumMaxflow
{
public:
	SubsetSumMaxflow(const SubsetSum &);
	
public:
	const SubsetSum &sss;
	double alpha;
	DiGraph gr;
	MED e2w;

public:
	int solve();
	int print();

public:
	int build_graph();
};

#endif
