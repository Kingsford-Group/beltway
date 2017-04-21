#include "subset_sum.h"
#include "subset_sum_maxflow.h"

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		printf("usage: %s <num-amino-acid> <spectrum-file>\n", argv[0]);
		return 0;
	}

    ifstream f(argv[2]);

    SubsetSum s = SubsetSum(atoi(argv[1]), f);
    s.clean();
    
    SubsetSumMaxflow maxflow(s);
    
    maxflow.solve();
    maxflow.print();

	return 0;
}
