#include "config.h"
#include <string>
#include <cstdlib>

double ilp_time_limit = 600;
bool no_infinity_contraints = false;
double max_error_allowed = -1;

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
        else if(string(argv[i]) == "-max_error")
        {
            max_error_allowed = atof(argv[i + 1]);
            i++;
        }
	}
	return 0;
}
