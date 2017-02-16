#ifndef __ALPHABET_H__
#define __ALPHABET_H__

#include <map>
#include <string>

using namespace std;

typedef map<string, double> MSD;
typedef pair<string, double> PSD;

class alphabet
{
public:
	alphabet(const string &file);

public:
	MSD a2w;

public:
	int read(const string &file);
	int print();
};

#endif
