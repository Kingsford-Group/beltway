#include <iostream>
#include "ilpsolver.h"
#include "config.h"

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc != 4)
	{
		cout << "usage: " << argv[0] << " <alphabet-file> <weights-file> <output-file>" << endl;
		return 0;
	}

	parse_arguments(argc, argv);

	ilpsolver sv(argv[1], argv[2]);
	sv.solve();
	sv.write(argv[3]);
	sv.print();

	return 0;
}
