#ifndef __CONFIG_H__
#define __CONFIG_H__

using namespace std;

extern double ilp_time_limit;
extern bool no_infinity_contraints;
extern double max_error_allowed;
extern bool no_priming_contraints;
extern bool lp_relax;
extern bool use_mvars;
int parse_arguments(int argc, const char ** argv);

#endif
