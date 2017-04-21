#ifndef __SUBSET_SUM_MAXFLOW_H__
#define __SUBSET_SUM_MAXFLOW_H__

#include <subset_sum.h>

using namespace std;

class SubsetSumMaxflow
{
public:
	SubsetSumMaxflow(const SubsetSum &);
	
public:
	const SubsetSum &sss;

public:
	int solve();
	int print();
};

#endif
