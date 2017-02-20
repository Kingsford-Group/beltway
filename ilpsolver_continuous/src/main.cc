#include <iostream>
#include "ilpsolver.h"
#include "config.h"

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc == 1)
	{
		cout << "usage: " << argv[0] << " <weights-file> <output-file> [-t ilp-time-limit]" << endl;
		return 0;
	}

	parse_arguments(argc, argv);

	try
	{
		ilpsolver sv(argv[1]);
		sv.solve();
		sv.write(argv[2]);
		sv.print();
	}
	catch(GRBException &e)
	{
		printf("%s\n", e.getMessage().c_str());
	}

	return 0;
}

