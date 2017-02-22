#include "config.h"
#include <string>
#include <cstdlib>

double ilp_time_limit = 600;
int mip_focus = 0;

int parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-t")
		{
			ilp_time_limit = atof(argv[i + 1]);
			i++;
		}
		if(string(argv[i]) == "-f")
		{
			mip_focus = atoi(argv[i + 1]);
			i++;
		}
	}
	return 0;
}
