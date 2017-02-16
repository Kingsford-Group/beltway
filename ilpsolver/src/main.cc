#include <iostream>
#include "alphabet.h"

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc != 3)
	{
		cout << "usage: " << argv[0] << " <alphabet-file> <weights-file>" << endl;
		return 0;
	}

	alphabet ab(argv[1]);
	ab.print();

	return 0;
}
