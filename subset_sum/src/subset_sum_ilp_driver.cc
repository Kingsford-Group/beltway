#include "subset_sum.h"
#include "subset_sum_ilp.h"

int main(int argc, char* argv[]){
	//Kpath::best_final_num_bits = -1;

    ifstream f(argv[2]);

    SubsetSum s = SubsetSum(atoi(argv[1]),f);
    s.clean();
    
    SubsetSumILP ilp = SubsetSumILP(&s);
    
    ilp.solve();
    ilp.print();
}
