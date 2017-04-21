#ifndef __SUBSET_SUM_MAXFLOW_H__
#define __SUBSET_SUM_MAXFLOW_H__

#include "subset_sum.h"
#include "mygraph.h"

using namespace std;

class SubsetSumMaxflow
{
public:
	SubsetSumMaxflow(const SubsetSum &);
	
public:
	const SubsetSum &sss;
	DiGraph gr;

public:
	int solve();
	int print();

public:
	int build_graph();
};

#endif
