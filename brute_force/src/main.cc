#include <iostream>
#include "brute_force.h"
#include "config.h"

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc == 1)
	{
		cout << "usage: " << argv[0] << " <alphabet-file> <weights-file> <output-file> [-t ilp-time-limit]" << endl;
		return 0;
	}

	parse_arguments(argc, argv);

	brute_force sv(argv[1], argv[2]);
	sv.solve();
	sv.write(argv[3]);
	sv.print();

	return 0;
}

