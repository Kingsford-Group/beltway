#ifndef __ILP_H__
#define __ILP_H__

#include <vector>
#include <set>
#include <map>
#include <string>

using namespace std;

typedef map<string, double> MSD;
typedef pair<string, double> PSD;

class brute_force
{
public:
	brute_force(const string &alphabet_file, const string &spectrum_file);
	~brute_force();

public:
	// instance
	int slots;
	double ubound;
	double current_best;
	vector<double> spectrum;
	vector<string> aa_list;
	vector<double> aa_mass;

	// testing assignments
	vector<int> xassign_test;
        vector<int> lassign_test;
        vector<int> uassign_test;

	// results
	vector<int> xassign;
	vector<int> lassign;
	vector<int> uassign;

public:
	int solve();
	int print();
	int write(const string &file);

private:
	// read input files
	int read_alphabet(const string &file);
	int read_spectrum(const string &file);
	int compute_upper_bound();

	bool increment_xassign_test(int);
	bool increment_lassign_test(int);
	bool increment_uassign_test(int);
	double calculate_objective();

};

#endif
