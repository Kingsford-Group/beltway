#include <iostream>
#include "alphabet.h"
#include "ilpsolver.h"

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc != 3)
	{
		cout << "usage: " << argv[0] << " <alphabet-file> <weights-file>" << endl;
		return 0;
	}

	ilpsolver sv(argv[1], argv[2]);
	sv.print();

	return 0;
}
