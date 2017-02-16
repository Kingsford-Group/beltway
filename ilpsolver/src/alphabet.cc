#include "alphabet.h"
#include <fstream>
#include <sstream>
#include <iostream>

alphabet::alphabet(const string &file)
{
	read(file);
}

int alphabet::read(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) 
	{
		cout << "open " << file.c_str() << " error" << endl;
		return 0;
	}

	string line;
	char a[1024];
	double w;
	while(getline(fin, line))
	{
		stringstream sstr(line);
		sstr >> a >> w;
		string s(a);
		if(a2w.find(s) == a2w.end()) a2w.insert(PSD(s, w));
		else a2w[s] = w;
	}
	return 0;
}

int alphabet::print()
{
	for(MSD::iterator it = a2w.begin(); it != a2w.end(); it++)
	{
		string s = it->first;
		double w = it->second;
		printf("alphabet %s -> %.2lf\n", s.c_str(), w);
	}
	return 0;
}
