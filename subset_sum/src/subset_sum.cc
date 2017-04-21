#include "subset_sum.h"
#include <map>

using namespace std;

int Kpath::best_final_num_bits = -1;
int Kpath::min_num_unmatched = -1;


void Kpath::count_bits(){
    if(private_num_bits != -1) return;
    private_num_bits = 0;
    for(int i=0;i<theoretical_spectrum_length;i++) private_num_bits += (theoretical_spectrum_coverage[i])?1:0;
}

int Kpath::num_bits(){
    if(private_num_bits == -1) count_bits();
    return private_num_bits;
}

void Kpath::set_num_bits(int in){
    private_num_bits = in;
}

std::ostream &operator<<(std::ostream &os, Kpath const &m) { 
    stringstream s; 
    for(int i=0;i<m.path.size();i++){
        if(i != 0) s << "-";
        s << m.path[i];
    }

    return os << s.str();
}

bool Kpath::compare_coverage(Kpath* a, Kpath* b){
        //cerr << "Compare " << *a << " and " << *b << endl;
        if(a->num_unmatched==b->num_unmatched){
                if(a->num_bits()==b->num_bits() && a->path.size() != b->path.size()){
                        return a->path.size()>b->path.size();
                }else if(a->num_bits()==b->num_bits() && a->path.size()>0 && b->path.size()>0){
                        int min = (a->path.size()>b->path.size())?b->path.size():a->path.size();
			for(int i=0;i<min;i++){
				if(a->path[i]!=b->path[i]) return a->path[i]<b->path[i];
			}
                }else{
                        return a->num_bits()>b->num_bits();
                }
        }else{
                return a->num_unmatched<b->num_unmatched;
        }
	return false;
}

bool Kpath::same_path(Kpath* a, Kpath* b){
	if(a->num_bits()==b->num_bits() &&
		a->path.size() == b->path.size() && 
		a->num_unmatched==b->num_unmatched){
		for(int i=0; i<a->path.size(); i++){
			if(a->path[i] != b->path[i]) return false;
		}
		return true;
	}
	return false;
}

Kpath::Kpath(Kpath* from){
    //path = from->path;
    path.clear();
    for(int i=0;i<from->path.size();i++) path.push_back(from->path[i]);
    //spectrum = from->spectrum;
    spectrum.clear();
    for(int i=0;i<from->spectrum.size();i++) spectrum.push_back(from->spectrum[i]);
    theoretical_spectrum_length = from->theoretical_spectrum_length;
    theoretical_spectrum_coverage = new bool[theoretical_spectrum_length];
    for(int i=0;i<theoretical_spectrum_length;i++) theoretical_spectrum_coverage[i] = from->theoretical_spectrum_coverage[i];
    private_num_bits = from->private_num_bits;
    spectrum_ends_with = from->spectrum_ends_with;
    num_unmatched = from->num_unmatched;
}

Kpath::Kpath(int spectrum_length){
	path.clear();
	spectrum.clear();
	spectrum_ends_with.clear();
	private_num_bits = -1;
	theoretical_spectrum_coverage = new bool[spectrum_length];
    for(int i=0;i<spectrum_length;i++) theoretical_spectrum_coverage[i] = false;
    theoretical_spectrum_length = spectrum_length;
    num_unmatched = 0;
}

Kpath::Kpath(int s, int spectrum_index, int spectrum_length){
    path.clear();
    path.push_back(s);
    spectrum.clear();
    spectrum.push_back(s);
    spectrum_ends_with.clear();
    spectrum_ends_with.push_back(0);
    private_num_bits = -1;
    theoretical_spectrum_coverage = new bool[spectrum_length];
    for(int i=0;i<spectrum_length;i++) theoretical_spectrum_coverage[i] = false;
    theoretical_spectrum_coverage[spectrum_index] = true;
    theoretical_spectrum_length = spectrum_length;
    num_unmatched = 0;
}

void Kpath::add_step(int s, int* reverse_spectrum){
    int spectrum_length = spectrum.size();
    path.push_back(s);
    int my_index = path.size()-1;
    spectrum.push_back(s);
    spectrum_ends_with.push_back(my_index);
    int input_index = reverse_spectrum[s];
    if(input_index > -1) theoretical_spectrum_coverage[input_index] = true;

    for(int i=0;i<spectrum_length;i++){
        if(spectrum_ends_with[i] == my_index-1){
            spectrum.push_back(spectrum[i] + s);
            spectrum_ends_with.push_back(my_index);
            int input_index = reverse_spectrum[spectrum[i] + s];
            if(input_index > -1) theoretical_spectrum_coverage[input_index] = true;
            else num_unmatched++;
        }
    }
    private_num_bits = -1;
}

string to_binary(int length, unsigned long long num){
    if(length == 1 && num == 1) return "1";
    if(length == 1 && num == 0) return "0";
    if(length == 1) return "ERROR";
    stringstream s;
    s << to_binary(length-1, num/2);
    if(num % 2 == 1) s << "1";
    else s << "0";
    return s.str();
}

void SubsetSum::subset_sum(){
        sum_array.resize(max_value+1);
        for(int i=1;i<=max_value;i++){
                //bit_masks[i] = 0;
                sum_array[i].clear();
                for(int j=0;j<spectrum.size();j++){
                        if(i==spectrum[j]){
                                sum_array[i].push_back(pointer(0,j));
                                //bit_masks[i] |= (unsigned long long)pow(2,j);
                        }else if(i > spectrum[j]){
                                if(sum_array[i - spectrum[j]].size() > 0){
                                        sum_array[i].push_back(pointer((i-spectrum[j]),j));
                                        //bit_masks[i] |= (unsigned long long)pow(2,j);
                                }
                        }
                }
		//cerr << "Sum Array: " << i << " is sized " << sum_array[i].size() << " and has bit mask " << to_binary(spectrum.size(),bit_masks[i]) << endl;
        }
}

void SubsetSum::find_paths(){
	paths = new bool**[2];
	paths[0] = new bool*[max_value];
	paths[1] = new bool*[max_value];
	for(int i=0;i<max_value;i++){
		paths[0][i] = new bool[k+1];
		paths[1][i] = new bool[k+1];
		for(int j=0;j<=k;j++) paths[0][i][j] = false;
		for(int j=0;j<=k;j++) paths[1][i][j] = false;
		for(int j=0;j<sum_array[i+1].size();j++){
			if(sum_array[i+1][j].next_index == 0){
				paths[0][i][1] = true;
			}else{
				for(int l=2;l<=k;l++){
					paths[0][i][l] |= paths[0][sum_array[i+1][j].next_index-1][l-1];
				}
			}
		}
		/*if(sum_array[i+1].size() > 0){
			cerr << i << ":";
			for(int j=0;j<=k;j++) cerr << "\t" << paths[0][i][j];
			cerr << endl;
		}*/
	}

	for(int i=0;i<max_value;i++){
                int inverse = max_value - i - 1;

		if(i==0){
			for(int j=0;j<sum_array[inverse+1].size();j++){
				//cerr << "Looking at " << sum_array[inverse+1][j].next_index << endl;
				if(sum_array[inverse+1][j].next_index != 0){
					paths[1][sum_array[inverse+1][j].next_index-1][1] = true;
				}else{
				    paths[1][max_value-1][1] = true;
				}
			}
		}else{
			for(int j=0;j<sum_array[inverse+1].size();j++){
				for(int l=1;l<k;l++){
					if(sum_array[inverse+1][j].next_index != 0){
						paths[1][sum_array[inverse+1][j].next_index-1][l+1] |= paths[1][inverse][l];
					}else{
					    paths[1][max_value-1][l+1] |= paths[1][inverse][l];
					}
				}
			}
		}

	}

	valid_to_include = new bool*[max_value];
	for(int i=0;i<max_value;i++){
	    valid_to_include[i] = new bool[k+1];
	    for(int set_index=0;set_index<=k;set_index++){
            valid_to_include[i][set_index] = true;
            valid_to_include[i][set_index] &= (paths[0][i][1]);
            valid_to_include[i][set_index] &= (paths[1][i][1]);
            for(int j=1;j<set_index;j++) valid_to_include[i][set_index] &= (paths[0][i][j]);
            for(int j=1;j<k - set_index;j++) valid_to_include[i][set_index] &= (paths[1][i][j]);
        }
    }


	/*
	cerr << "------" << endl;
	for(int i=0;i<max_value;i++){
		bool sum = false;
		for(int j=0;j<=k;j++) sum |= paths[1][i][j];
		if(sum){
			cerr << i << ":";
                        for(int j=0;j<=k;j++) cerr << "\t" << paths[1][i][j];
                        cerr << endl;
		}
	}*/
}

vector< Kpath* > SubsetSum::recursive_set_recover(Recursive_call* call){
	//cerr << "Zero" << endl;
	Kpath* so_far = call->so_far;
	//cerr << "One" << endl;
	int set_index = call->set_index;
	//cerr << "Two" << endl;
	int size_index = call->size_index;
	call->processed = true;
	//cerr << "A " << call->so_far->path.size() << endl;
	//cerr << "Recursive Call "<< call << " " << *so_far << "/" << so_far<< "(" << so_far->num_unmatched << ")," << set_index << "," << size_index << "\tsum_array[" << size_index << "].size = " << sum_array[size_index].size() << "\tBest: " << Kpath::best_final_num_bits << endl;
	vector< Kpath* > r;
	r.clear();
	 if(Kpath::min_num_unmatched != -1 && Kpath::min_num_unmatched < so_far->num_unmatched) return r;
	vector< Kpath* > new_paths;
	new_paths.clear();
	//if(size_index % (max_value/2) != 0) new_bit_vector |= bit_masks[size_index];
	for(int i=0;i<sum_array[size_index].size();i++){
		if(sum_array[size_index][i].next_index == 0){
			if(set_index == 1){// && new_bit_vector == (unsigned long long)(pow(2,spectrum.size())-1)){
				Kpath* k = new Kpath(so_far);
				k->add_step(spectrum[sum_array[size_index][i].spectrum_index],reverse_spectrum);
				if(k->num_bits() >= Kpath::best_final_num_bits && (Kpath::min_num_unmatched == -1 || Kpath::min_num_unmatched>k->num_unmatched)){
				    r.push_back(k);
				    Kpath::best_final_num_bits = k->num_bits();
				    Kpath::min_num_unmatched = k->num_unmatched;
				    cerr << "New Best: " << Kpath::best_final_num_bits << "\tUnmatched: " << Kpath::min_num_unmatched << "\tPath: " << *k <<  endl;
				}else{
				    //delete k;
				}
				//r.push_back(new Kpath(spectrum[sum_array[size_index][i].spectrum_index],sum_array[size_index][i].spectrum_index,spectrum.size()));
				return r;
			}
		}else{
			bool okay = valid_to_include[size_index-1][set_index];
			if(okay){
				Kpath* k = new Kpath(so_far);
				k->add_step(spectrum[sum_array[size_index][i].spectrum_index],reverse_spectrum);
				if(Kpath::min_num_unmatched == -1 || Kpath::min_num_unmatched >= k->num_unmatched){
					new_paths.push_back(k);
					call_queue.push_back(new Recursive_call(k, set_index-1, sum_array[size_index][i].next_index));
				}
				//else delete k;
			}
		}
	}
	//delete so_far;
	return r;
	//cerr << "Before sort " << endl;
	sort(new_paths.begin(),new_paths.end(),Kpath::compare_coverage);
	for(int i=0;i<new_paths.size();i++){
		call_queue.push_back(new Recursive_call(new_paths[i], set_index-1, size_index - new_paths[i]->path[new_paths[i]->path.size()-1]));
		//vector< Kpath* > t = recursive_set_recover(new_paths[i], set_index-1, size_index - new_paths[i]->path[new_paths[i]->path.size()-1]);
		//r.insert(r.end(),t.begin(),t.end());
	}
	return r;
}

void SubsetSum::recover_sets(){
	call_queue.clear();
	//kcycles = recursive_set_recover(new Kpath(spectrum.size()), k, max_value);
	call_queue.push_back(new Recursive_call(new Kpath(spectrum.size()), k, max_value));
	bool keep_going = true;
	while(keep_going){
		sort(call_queue.begin(),call_queue.end(),Recursive_call::Compare);
		int i=0;
		while(call_queue[i]->processed) i++;
		keep_going = (Kpath::min_num_unmatched == -1 || Kpath::min_num_unmatched > call_queue[i]->so_far->num_unmatched);
		//cerr << "Call Queue Size: " << call_queue.size() << endl;
		if(keep_going){
			vector< Kpath* > t = recursive_set_recover(call_queue[i]);	
		
			//cerr << "Call Queue Size: " << call_queue.size() << endl;
			kcycles.insert(kcycles.end(),t.begin(),t.end());
		}
	}
	/*for(int i=0;i<kcycles.size();i++){
		vector<int> temp = kcycles[i];
		int min = 0;
		for(int j=1; j<kcycles[i].size(); j++){
			//cerr << kcycles[i][j] << " < " << kcycles[i][min] << endl;
			if(kcycles[i][j] < kcycles[i][min]) min = j;
		}
		for(int j=0; j<kcycles[i].size(); j++){
			cerr << "mapping " << ((j+min)%kcycles[i].size()) << "->" << j << "(" << kcycles[i][j] << "->" << temp[(j+min)%kcycles[i].size()] << ")" << endl;
			kcycles[i][j] = temp[(j+min)%kcycles[i].size()];
		}
	}*/
}

bool SubsetSum::match_by_rotation_reversal(int a, int b){

    for(int i=0; i<kcycles[a]->path.size(); i++){
        if(kcycles[a]->path[i] == kcycles[b]->path[0]){
            bool full_match = true;
            bool full_match_reverse = true;
            int jr = 0;
            for(int j=0; j<kcycles[b]->path.size(); j++){
                int ip = (i + j) % (kcycles[a]->path.size());
                full_match &= (kcycles[a]->path[ip] == kcycles[b]->path[j]);

                full_match_reverse &= (kcycles[a]->path[ip] == kcycles[b]->path[jr]);
                jr--;
                if(jr == -1) jr += kcycles[b]->path.size();
            }
            if(full_match || full_match_reverse) return true;
        }
    }
    return false;
}

int SubsetSum::verify_with_theoretical_spectrum(int a){
    vector<int> theoretical;
    theoretical.clear();
    for(int i=0;i<kcycles[a]->path.size();i++){
        theoretical.push_back(kcycles[a]->path[i]);
        int sum = kcycles[a]->path[i];
        for(int j=1;j<kcycles[a]->path.size();j++){
            int l = (i+j)%kcycles[a]->path.size();
            sum += kcycles[a]->path[l];
            theoretical.push_back(sum);
        }
    }
    sort(theoretical.begin(),theoretical.end());
    vector<int>::iterator it = std::unique (theoretical.begin(), theoretical.end());
    spectrum.resize( std::distance(theoretical.begin(),it) );

	sort(theoretical.begin(), theoretical.end());

    /*int** matching = new int*[spectrum.size()];
    int** matching_dir = new int*[spectrum.size()];
    for(int i=0;i<spectrum.size();i++){
        matching[i] = new int[theoretical.size()];
        matching_dir[i] = new int[theoretical.size()];
        for(int j=0;j<theoretical.size();j++){
            if(i==0 || j==0){
                matching[i][j] = (spectrum[i] == theoretical[j])?1:0;
                matching_dir[i][j] = 0;
                if(j!=0){
                    if(matching[i][j-1] > matching[i][j]) matching[i][j] = matching[i][j-1];
                    matching_dir[i][j] = 1;
                }
                if(i!=0){
                    if(matching[i-1][j] > matching[i][j]) matching[i][j] = matching[i-1][j];
                    matching_dir[i][j] = 2;
                }
            }else{
                matching[i][j] = matching[i-1][j-1] + ((spectrum[i] == theoretical[j])?1:0);
                matching_dir[i][j] = 0;
                if(matching[i-1][j] > matching[i][j]){
                    matching[i][j] = matching[i-1][j];
                    matching_dir[i][j] = 2;
                }
                if(matching[i][j-1] > matching[i][j]){
                    matching[i][j] = matching[i][j-1];
                    matching_dir[i][j] = 1;
                }
            }
            //cerr << matching[i][j] << "(" << (matching_dir[i][j]==0?"M":(matching_dir[i][j]==1)?"<":"^") << ")" << " ";
        }
       // cerr << endl;
    }
    */
    /*for(int i=spectrum.size()-1;i>=0;i--) cerr << (kcycles[a]->theoretical_spectrum_coverage[i]?"1":"0");
    cerr << endl;
    int i = spectrum.size()-1;
    int j = theoretical.size()-1;
    while(i >= 0 || j >= 0){
        if(matching_dir[i][j] == 0){
            cerr << (spectrum[i] == theoretical[j])?"1":"0";
            i--;
            j--;
        }else if(matching_dir[i][j] = 1){
            j--;
        }else if(matching_dir[i][j] = 2){
            i--;
            cerr << "0";
        }
    }
    cerr << endl;
    exit(0);*/
    bool* matched = new bool[spectrum.size()];
    for(int i=0;i<spectrum.size();i++) matched[i] = false;
    for(int i=0;i<theoretical.size();i++){
        if(reverse_spectrum[theoretical[i]] != -1){
            matched[reverse_spectrum[theoretical[i]]] = true;
        }
    }
    int num_bits = 0;
    for(int i=0;i<spectrum.size();i++) num_bits += (matched[i]?1:0);
    return num_bits;
    //kcycles[a]->set_num_bits(num_bits);
    //return matching[spectrum.size()-1][theoretical.size()-1];
    /*if(spectrum.size()==theoretical.size()){
        for(int i=0;i<spectrum.size();i++){
            if(spectrum[i] != theoretical[i]) return false;
        }
        return true;
    }
    return false;*/
}

SubsetSum::SubsetSum(int kin, istream& f){
   k = kin;
   spectrum.clear();

   while(!f.eof()){
        string line;
        getline(f,line);
        if(line[0] != '#' && line.length()>=1){
            spectrum.push_back(atoi(line.c_str()));
        }
    }

	sort(spectrum.begin(), spectrum.end());

    vector<int>::iterator it = std::unique (spectrum.begin(), spectrum.end());
    spectrum.resize( std::distance(spectrum.begin(),it) );

	sort(spectrum.begin(), spectrum.end());
	max_value = spectrum[spectrum.size()-1];

	reverse_spectrum = new int[max_value+1];
	for(int i=0;i<=max_value;i++) reverse_spectrum[i] = -1;
	for(int i=0;i<spectrum.size();i++) reverse_spectrum[spectrum[i]] = i;


    subset_sum();

}

void SubsetSum::prune_edges(){
    int locations_removed = 0;
    for(int i=0;i<max_value;i++){
        bool keep = paths[0][i][k] || paths[1][i][k];
        
	    for(int j=1;j<k;j++){
	        keep |= (paths[0][i][j] && paths[0][i][k-j]);
	    }
	    
	    if(!keep){
	        sum_array[i].clear();
	        locations_removed++;
	    }
    }
    cout << "Locations pruned: " << locations_removed << endl;
}

void SubsetSum::clean(){
    cerr << "Between Sum and Find Paths..." << endl;
	find_paths();
	
	prune_edges();
}

void SubsetSum::recover_k_paths(){

    cerr << "Between Find Paths and Recover..." << endl;
	recover_sets();
	cerr << "Between Recover and Distinct..." << endl;

	sort(kcycles.begin(),kcycles.end(),Kpath::compare_coverage);

	cerr << "K Cycle Length: " << kcycles.size() << endl;
	bool* distinct = new bool[kcycles.size()];
	int max_verified = 0;
	int* verified = new int[kcycles.size()];

	distinct[0] = true;
	for(int j=0;j<kcycles.size();j++){
	    cerr << "On cycle " << j << "/" << kcycles.size() << " ("<< ((j*10000)/(kcycles.size())/100.0) << "%)\r";
	    distinct[j] = true;
	    for(int l=0;l<j;l++){
	        distinct[j] &= !(match_by_rotation_reversal(j,l));
	    }
	    if(distinct[j]){
	        verified[j] = verify_with_theoretical_spectrum(j);
	        if(verified[j]>max_verified){
	            max_verified = verified[j];
	            cerr << j << "::";
	            for(int i=0;i<kcycles[j]->path.size();i++){
                    if(i!=0) cerr << "-";
                    cerr << kcycles[j]->path[i];
                }
                cerr << ": " << verified[j] << "/" << spectrum.size() << " bits covered: " << kcycles[j]->num_bits() << endl;
                if(verified[j] == spectrum.size()){
                    for(int i=0;i<kcycles[j]->path.size();i++){
                        if(i!=0) cout << "-";
                        cout << kcycles[j]->path[i];
                    }
                    cout << endl;
                    exit(0);

                }
	        }
	    }else{
	        verified[j] = 0;
	    }
	}


	cerr << "Between Distinct and Verify..." << endl;


	for(int i=0;i<kcycles.size();i++){
	    if(verified[i] == max_verified){
            for(int j=0;j<kcycles[i]->path.size();j++){
                if(j!=0) cout << "-";
                cout << kcycles[i]->path[j];
            }
            cout << ": " << verified[i] << endl;
		}
	}
	
}

void SubsetSum::analyze_degree_distribution() const
{
	map<int, int> m;
	for(int i = 0; i < sum_array.size(); i++)
	{
		int k = sum_array[i].size();
		if(m.find(k) == m.end()) m.insert(pair<int, int>(k, 1));
		else m[k]++;
	}
	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		int k = it->first;
		int c = it->second;
		printf("sumsetsum graph %d vertices has degree of %d\n", c, k);
	}
}
