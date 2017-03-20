#include "config.h"
#include <string>
#include <cstdlib>

double ilp_time_limit = 600;
bool no_infinity_contraints = false;
bool no_priming_contraints = false;
double max_error_allowed = -1;
bool lp_relax = false;
bool use_mvars = true;

int parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-t")
		{
			ilp_time_limit = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-no_infinity")
        {
            no_infinity_contraints = true;
        }		
        else if(string(argv[i]) == "-no_priming")
        {
            no_priming_contraints = true;
        }
        else if(string(argv[i]) == "-max_error")
        {
            max_error_allowed = atof(argv[i + 1]);
            i++;
        }
        else if(string(argv[i]) == "-lp_relaxation"){
            lp_relax = true;
        }
        else if(string(argv[i]) == "-ulvars"){
            use_mvars = false;
        }
        else if(string(argv[i]) == "-mvars"){
            use_mvars = true;
        }
	}
	return 0;
}
