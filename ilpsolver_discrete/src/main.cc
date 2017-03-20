#include <iostream>
#include "ilpsolver.h"
#include "config.h"

using namespace std;


vector<double> slice(const vector<double>& v, int start=0, int end=-1) {
    int oldlen = v.size();
    int newlen;

    if (end == -1 or end >= oldlen){
        newlen = oldlen-start;
    } else {
        newlen = end-start;
    }

    vector<double> nv(newlen);

    for (int i=0; i<newlen; i++) {
        nv[i] = v[start+i];
    }
    return nv;
}

int main(int argc, const char ** argv)
{
	if(argc == 1)
	{
		cout << "usage: " << argv[0] << " <alphabet-file> <weights-file> <output-file> [-t ilp-time-limit] [-max_error maximum-allowed-greedy-search-error] [-no_infinity] " << endl;
		return 0;
	}

	parse_arguments(argc, argv);

	try
	{
		ilpsolver sv(argv[1], argv[2]);
		
		/*vector<double> spectrum_save = sv.spectrum;
		
		//sv.spectrum = slice(sv.spectrum,0,sv.spectrum.size()/2);
		
		for(int i=sv.slots/2;i<spectrum_save.size();i+=sv.slots/2){
            sv.spectrum = slice(spectrum_save,0,i);
        
            sv.solve();
            sv.write(argv[3]);
            sv.print();
            sv.reset();
		}
		sv.spectrum = spectrum_save;*/
		
		//sv.greedy_warm_start();
		double best_error = -1;
		if(!no_priming_contraints) best_error = sv.graph_greedy_warm_start();
		printf("Best Error At End: %.3lf\n",best_error);
		if(best_error == -1 || best_error > 0.001) sv.solve();
        sv.write(argv[3]);
        sv.print();
        sv.reset();
	}
	catch(GRBException &e)
	{
		printf("%s\n", e.getMessage().c_str());
	}

	return 0;
}

