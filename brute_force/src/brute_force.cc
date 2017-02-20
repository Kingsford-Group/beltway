#include "brute_force.h"
#include "config.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h> 

brute_force::brute_force(const string &alphabet_file, const string &spectrum_file)
{
	read_alphabet(alphabet_file);
	read_spectrum(spectrum_file);
	compute_upper_bound();

}

brute_force::~brute_force()
{
}

int brute_force::solve()
{

	bool first = true;
	xassign_test.clear();
	xassign_test.resize(slots);
	do{
		lassign_test.clear();
		lassign_test.resize(spectrum.size());
		do{
	                printf("Current best objective: %lf\n", current_best);
			uassign_test.clear();
        	        uassign_test.resize(spectrum.size());
	                do{
				double test_objective = calculate_objective();
				if(first || test_objective < current_best){
					first = false;
					current_best = test_objective;
					xassign = xassign_test;
					lassign = lassign_test;
					uassign = uassign_test;
				}
			} while(increment_uassign_test(0));
		}while(increment_lassign_test(0));
	}while(increment_xassign_test(0));
	
	return 0;
}

bool brute_force::increment_xassign_test(int index){
	if(index >= slots) return false;
	if(xassign_test[index] < aa_list.size() - 1){
		xassign_test[index]++;
		return true;	
	}
	xassign_test[index] = 0;
	return increment_xassign_test(index + 1);
}

bool brute_force::increment_lassign_test(int index){
        if(index >= spectrum.size()) return false;
        if(lassign_test[index] < slots - 1){
                lassign_test[index]++;
                return true;
        }
        lassign_test[index] = 0;
        return increment_lassign_test(index + 1);
}

bool brute_force::increment_uassign_test(int index){
        if(index >= spectrum.size()) return false;
        if(uassign_test[index] < slots - 1){
                uassign_test[index]++;
                return true;
        }
        uassign_test[index] = 0;
        return increment_uassign_test(index + 1);
}


double brute_force::calculate_objective(){
	double error = 0;


	for(int i = 0; i < spectrum.size(); i++)
	{
		double w = 0;
		//printf("i: %d\tl: %d\tu: %d\tx: %d\tmass: %lf\n",i,lassign_test[i],uassign_test[i],xassign_test[lassign_test[i]], aa_mass[xassign_test[lassign_test[i]]]);
		if(lassign_test[i]<=uassign_test[i]){
			for(int l=lassign_test[i]; l<=uassign_test[i]; l++){
				w += aa_mass[xassign_test[l]];
			}
		}else{	
			for(int l=lassign_test[i]; l<=slots; l++){
                                w += aa_mass[xassign_test[l]];
                        }
			for(int l=0; l<=uassign_test[i]; l++){
                                w += aa_mass[xassign_test[l]];
                        }
		}

		//printf("Error: %lf\tSpectrum: %lf\tw: %lf\t",error,spectrum[i],w);
		error += (abs(spectrum[i] - w));
	}
	return error;
}




int brute_force::read_alphabet(const string &file)
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
		aa_list.push_back(a);
		aa_mass.push_back(w);
		/*
		string s(a);
		if(aa2m.find(s) == aa2m.end()) aa2m.insert(PSD(s, w));
		else aa2m[s] = w;
		*/
	}

	fin.close();
	return 0;
}

int brute_force::read_spectrum(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) 
	{
		cout << "open " << file.c_str() << " error" << endl;
		return 0;
	}

	string line;
	char buf[10240];
	while(getline(fin, line))
	{
		if(line == "") continue;
		stringstream sstr(line);
		sstr >> buf >> buf >> slots;
		break;
	}

	spectrum.clear();
	double m;
	while(getline(fin, line))
	{
		if(line == "") continue;
		stringstream sstr(line);
		sstr >> m;
		spectrum.push_back(m);
	}

	fin.close();
	return 0;
}

int brute_force::compute_upper_bound()
{
	ubound = 0;
	for(int i = 0; i < spectrum.size(); i++)
	{
		if(spectrum[i] > ubound) ubound = spectrum[i];
	}
	ubound = ubound * slots;
	return 0;
}

int brute_force::print()
{
	/*
	for(MSD::iterator it = aa2m.begin(); it != aa2m.end(); it++)
	{
		string s = it->first;
		double m = it->second;
		printf("amino acid %s -> %.3lf\n", s.c_str(), m);
	}
	*/
	for(int i = 0; i < aa_list.size(); i++)
	{
		printf("amino acid %d: %s -> %.3lf\n", i, aa_list[i].c_str(), aa_mass[i]);
	}
	printf("number of amino acid = %d\n", slots);
	for(int i = 0; i < spectrum.size(); i++)
	{
		printf("spectrum: %.3lf\n", spectrum[i]);
	}
	printf("upper bound = %.3lf\n", ubound);

	for(int i = 0; i < xassign.size(); i++)
	{
		int k = xassign[i];
		printf("slot %d is assigned amino acid %d: %s -> %.3lf\n", i, k, aa_list[k].c_str(), aa_mass[k]);
	}
	for(int i = 0; i < spectrum.size(); i++)
	{
		int l = lassign[i];
		int u = uassign[i];
		double w = 0;//wassign[i];
		double e = 0;//eassign[i];
		printf("spectrum %d with mass %.3lf is assigned to interval [%d, %d], with mass %.3lf and error %.3lf\n", i, spectrum[i], l, u, w, e);
	}
	return 0;
}

int brute_force::write(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail()) return 0;

	for(int i = 0; i < xassign.size(); i++)
	{
		int k = xassign[i];
		fout << aa_list[k].c_str() << endl;
	}

	fout.close();
	return 0;
}
