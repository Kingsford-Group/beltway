#ifndef __SUB_SUM_H__
#define __SUB_SUM_H__

#include "pointer.h"

#include <vector>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <queue>

using namespace std;


class Kpath{
    public:
        vector<int> path;
        bool* theoretical_spectrum_coverage;
        int theoretical_spectrum_length;
        vector<int> spectrum;
        vector<int> spectrum_ends_with;

        Kpath(Kpath*);
        Kpath(int, int, int);
        Kpath(int);
        void add_step(int, int*);

        static bool compare_coverage(Kpath* a, Kpath* b);
	static bool same_path(Kpath*, Kpath*);
	static int best_final_num_bits;

        int num_bits();
        void set_num_bits(int);

        int num_unmatched;
        static int min_num_unmatched;
        
    private:
        int private_num_bits;
        void count_bits();
};

class Recursive_call{
        public:
                Recursive_call( Kpath* a, int b, int c){ processed = false; so_far = a; set_index = b; size_index = c; };
		Kpath* so_far;
                int set_index;
                int size_index;
		bool processed;

                static bool Compare(Recursive_call* a, Recursive_call* b){ return Kpath::compare_coverage(a->so_far,b->so_far); };
                static bool Equal(Recursive_call* a, Recursive_call* b){ return Kpath::same_path(a->so_far,b->so_far); };
};

class SubsetSum{

    public:
        SubsetSum(int, istream&);
        
    	int k;
    	vector<int> spectrum;
        vector< vector< pointer > > sum_array;
        int* reverse_spectrum;
        bool*** paths;
        vector< Kpath* > kcycles;
        int max_value;
        bool** valid_to_include;
	    vector< Recursive_call* >  call_queue;
	    
	    void recover_k_paths();

        vector< Kpath* > recursive_set_recover(Recursive_call*);
        void recover_sets();
        void find_paths();
        void subset_sum();
        bool match_by_rotation_reversal(int, int);
        int verify_with_theoretical_spectrum(int a);
        
        void clean();
        void prune_edges();

};

#endif
