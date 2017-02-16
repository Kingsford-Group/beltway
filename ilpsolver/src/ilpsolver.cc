#include "ilpsolver.h"
#include <iostream>
#include <fstream>
#include <sstream>

ilpsolver::ilpsolver(const string &alphabet_file, const string &spectrum_file)
{
	read_alphabet(alphabet_file);
	read_spectrum(spectrum_file);

	env = new GRBEnv();
	model = new GRBModel(*env);
}

ilpsolver::~ilpsolver()
{
	if(model != NULL) delete model;
	if(env != NULL) delete env;
}

int ilpsolver::solve()
{
	add_amino_acid_variables();
	add_lower_endpoints_variables();
	add_upper_endpoints_variables();
	add_range_variables();
	add_error_variables();
	return 0;
}

int ilpsolver::read_alphabet(const string &file)
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
	return 0;
}

int ilpsolver::read_spectrum(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) 
	{
		cout << "open " << file.c_str() << " error" << endl;
		return 0;
	}

	string line;
	while(getline(fin, line))
	{
		if(line == "") continue;
		stringstream sstr(line);
		sstr >> slots;
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
	return 0;
}

int ilpsolver::add_amino_acid_variables()
{
	xvars.clear();
	xvars.resize(slots);
	for(int i = 0; i < slots; i++)
	{
		for(int k = 0; k < aa_list.size(); k++)
		{
			GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
			xvars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_lower_endpoints_variables()
{
	lvars.clear();
	lvars.resize(spectrum.size());
	for(int i = 0; i < spectrum.size(); i++)
	{
		for(int k = 0; k < aa_list.size(); k++)
		{
			GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
			lvars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_upper_endpoints_variables()
{
	uvars.clear();
	uvars.resize(spectrum.size());
	for(int i = 0; i < spectrum.size(); i++)
	{
		for(int k = 0; k < slots; k++)
		{
			GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
			uvars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_range_variables()
{
	rvars.clear();
	rvars.resize(slots);
	for(int i = 0; i < slots; i++)
	{
		for(int k = 0; k < slots; k++)
		{
			GRBVar var = model->addVar(0, ubound, 0, GRB_CONTINUOUS);
			rvars[i].push_back(var);
		}
	}
	model->update();
	return 0;
}

int ilpsolver::add_error_variables()
{
	evars.clear();
	for(int i = 0; i < spectrum.size(); i++)
	{
		GRBVar var = model->addVar(0, ubound, 0, GRB_CONTINUOUS);
		evars.push_back(var);
	}
	model->update();
	return 0;
}

int ilpsolver::print()
{
	/*
	for(MSD::iterator it = aa2m.begin(); it != aa2m.end(); it++)
	{
		string s = it->first;
		double m = it->second;
		printf("amino acid %s -> %.2lf\n", s.c_str(), m);
	}
	*/
	assert(aa_list.size() == aa_mass.size());
	for(int i = 0; i < aa_list.size(); i++)
	{
		printf("amino acid %d: %s -> %.2lf\n", i, aa_list[i].c_str(), aa_mass[i]);
	}
	printf("number of amino acid = %d\n", slots);
	for(int i = 0; i < spectrum.size(); i++)
	{
		printf("spectrum: %.2lf\n", spectrum[i]);
	}
	return 0;
}
