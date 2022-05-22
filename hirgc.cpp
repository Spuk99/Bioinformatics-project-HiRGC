#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <sys/time.h>

using namespace std;

//Target sequence
string id_tg;           // identifier from target FASTA file
vector<int> t_seq_len;  // vector of sequence lengths for each line
string t_seq_L = "";    // all lines of sequences concatenated into one string to form one sequence
vector<int> t_low_pos;  // vector of positions of lowercase letters in target sequence
vector<int> t_low_len;  // vector of lengths of lowercase letters in target sequence
string t_seq_L1 = "";   // sequence with lowercase letters converted to uppercase
vector<int> t_N_pos;    // vector of positions of N letters in target sequence
vector<int> t_N_len;    // vector of lengths of N letters in target sequence
vector<int> t_oth_pos;  // vector of positions of other letters in target sequence
vector<char> t_oth_ch;  // vector of characters of other letters in target sequence
vector<int> t_oth_len;  // vector of lengths of other letters in target sequence
string t_seq_L3 = "";   // sequence with other letters than A, C, G and T removed from target sequence
string t_final="";      // final sequence encoded - A=0, C=1, G=2, T=3 

// Reference sequence
string id_r;            // identifier from reference FASTA file
vector<int> r_seq_len;  // vector of sequence lengths for each line
string r_seq_L = "";    // all lines of sequence concatenated into one string to form one sequence
string r_seq_L1 = "";   // sequence with lowercase letters converted to uppercase
string r_seq_L3 = "";   // sequence with other letters than A, C, G and T removed from reference sequence
string r_final="";      // final reference sequence encoded in the same way as the target sequence

const int k = 20;  // This is said to be the best value for k in the article
const int max_hash_size = 1 << 20; // Maximum length of the hash table
const int max_char_num = 1 << 29; // Maximum length of a chromosome

int *point = new int[max_hash_size]; // An array of entries 
int *loc = new int[max_char_num]; // An array of header pointers



// Written by Katarina Misura
// This function reads the target sequence from the FASTA file and stores the final encoded sequence in the variable t_final
void target_preprocess(string file_name){
    //reading target FASTA file
    string line;
    ifstream myfile(file_name);
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            if (line[0] == '>')
            {
                id_tg = line;
            }
            else
            {
				if(line[line.size() - 1] == '\r'){
					line.erase(line.size() - 1);
				}
                t_seq_L += line;
                t_seq_L +=')'; //adding the sing for end of the line (we used the right bracket as the end of the line)
                t_seq_len.push_back(line.length());
            }

        }
        myfile.close();
    }
    else
    {
        cout << "Unable to open file";
        throw runtime_error("Closing program!");
    }  
    //preprocesing target sequence
    int interval; // interval 
    int last = -1; //flag for last position
    //turning lowercase letters to uppercase
    for (int i = 0; i < t_seq_L.length(); i++)
    {
        if (t_seq_L[i] >= 'a' && t_seq_L[i] <= 'z' && last == -1)
        {
            last = i;
            t_low_pos.push_back(last);
        }
        else if (last != -1 && !(t_seq_L[i] >= 'a' && t_seq_L[i] <= 'z'))
        {
            interval = i - last;
            last = -1;
            t_low_len.push_back(interval);
        }
        t_seq_L1 += toupper(t_seq_L[i]);
    }
    if (last != -1)
    {
        interval = t_seq_L.length() - last;
        t_low_len.push_back(interval);
    }

    //collecting information about letters N
    last = -1;
    for (int i = 0; i < t_seq_L1.length(); i++)
    {
        if (t_seq_L1[i] == 'N' && last == -1)
        {
            last = i;
            t_N_pos.push_back(last);
        }
        else if (last != -1 && t_seq_L1[i] != 'N')
        {
            interval = i - last;
            last = -1;
            t_N_len.push_back(interval);
        }
    }
    if (last != -1)
    {
        interval = t_seq_L1.length() - last;
        t_N_len.push_back(interval);
    }

    //collecting information about other letters than A, C, G, T and N and removing them from the sequence
    last = -1;
    for (int i = 0; i < t_seq_L1.length(); i++)
    {
        if (last == -1 &&!(t_seq_L1[i] == 'A' || t_seq_L1[i] == 'C' || t_seq_L1[i] == 'G' || t_seq_L1[i] == 'T'|| t_seq_L1[i] == 'N') )
        {
            last = i;
            t_oth_pos.push_back(last);
            t_oth_ch.push_back(t_seq_L1[i]);
        }else if(last != -1 && !(t_seq_L1[i] == 'A' || t_seq_L1[i] == 'C' || t_seq_L1[i] == 'G' || t_seq_L1[i] == 'T' || t_seq_L1[i] == 'N')){
            t_oth_ch.push_back(t_seq_L1[i]);
        }
        else if (last != -1 && (t_seq_L1[i] == 'A' || t_seq_L1[i] == 'C' || t_seq_L1[i] == 'G' || t_seq_L1[i] == 'T' || t_seq_L1[i] == 'N'))
        {
            interval = i - last;
            last = -1;
            t_oth_len.push_back(interval);
        }
        if ((t_seq_L1[i] == 'A' || t_seq_L1[i] == 'C' || t_seq_L1[i] == 'G' || t_seq_L1[i] == 'T'))
        {
            t_seq_L3 += t_seq_L1[i];
        }
    }
    if (last != -1)
    {
        interval = t_seq_L1.length() - last;
        t_oth_len.push_back(interval);
    }

    //encoding target sequence
    for(int i=0; i<t_seq_L3.length();i++){
        if(t_seq_L3[i]=='A'){
            t_final+="0";
        }
        else if(t_seq_L3[i]=='C'){
            t_final+="1";
        }
        else if(t_seq_L3[i]=='G'){
            t_final+="2";
        }
        else if(t_seq_L3[i]=='T'){
            t_final+="3";
        }
    }

}

// Written by Katarina Misura
// This function reads the reference sequence from the FASTA file and stores the final encoded sequence in the variable r_final
void refrence_preprocess(string file_name){
    //reading reference FASTA file
    string line;
    ifstream myfile(file_name);
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            if (line[0] == '>')
            {
                id_r = line;
            }
            else
            {
				if(line[line.size() - 1] == '\r'){
					line.erase(line.size() - 1);
				}
                r_seq_L += line;
                r_seq_len.push_back(line.length());
            }
        }
        myfile.close();
    }
    else
    {
        cout << "Unable to open file";
        throw runtime_error("Closing program!");
    }
    //preprocesing reference sequence
    //turning lowercase letters to uppercase
    for (int i = 0; i < r_seq_L.length(); i++)
    {
        r_seq_L1 += toupper(r_seq_L[i]);
    }
    //removing all letters than A, C, G, T from reference sequence
    for (int i = 0; i < r_seq_L1.length(); i++)
    {
        if ((r_seq_L1[i] == 'A' || r_seq_L1[i] == 'C' || r_seq_L1[i] == 'G' || r_seq_L1[i] == 'T'))
        {
            r_seq_L3 += r_seq_L1[i];
        }
    }
    //encoding reference sequence
    for(int i=0; i<r_seq_L3.length();i++){
        if(r_seq_L3[i]=='A'){
            r_final+="0";
        }
        else if(r_seq_L3[i]=='C'){
            r_final+="1";
        }
        else if(r_seq_L3[i]=='G'){
            r_final+="2";
        }
        else if(r_seq_L3[i]=='T'){
            r_final+="3";
        }
    }
}

// Written by Katarina Misura
//function for writing auxiliary information to file
void saveDataToFile(ofstream &myfile){
    //write size of target sequence
    myfile << t_seq_L.length() << endl;

    //write positions and lengths of intervals of lowercase letters to file
    for(int i=0; i<t_low_pos.size();i++){
        myfile << t_low_pos[i] << "-"<<t_low_len[i] << " ";
    }
    myfile << endl;

    //write positions and lengths of intervals of letters N to file
    for(int i=0; i<t_N_pos.size();i++){
        myfile << t_N_pos[i] << "-"<<t_N_len[i] << " ";
    }
    myfile << endl;
    int j=0;
    int k=0;
    //write lengths of intervals of other letters and the interval to file
    for(int i=0; i<t_oth_pos.size();i++){
        myfile << t_oth_pos[i] << "-";
        for(; j<t_oth_len[i]+k;){
            myfile << t_oth_ch[j];
            j++;
        }
        k=j;
        myfile << " ";
    }
    myfile << endl;    
}



// Written by Marko Marfat
// This function constructs a hash table from the tuple values
// of a reference sequence. 
void initHT(){
	// Initialize entries
	for (int i = 0; i < max_hash_size; i++){
        	point[i] = -1;
    }

	uint64_t value = 0;
	
	// Only k values is needed, skip rest
	for (int i = 0; i < k - 1; i++){
		value = value << 2; // Not 100% sure about this one, I think it shifts the value so new char can be added (not really explained in the article)
		value = value + r_final[i]; // Add char to tuple
	}
    
    uint64_t remain = ((uint64_t) 1 << (2 * k)) - 1; // Used as a mask for which bits will remain (only those below 2 * k)
	for (int i = k - 1; i < r_final.length(); i++){
		value = value << 2;
		value = value + r_final[i];
		value = value & remain; // Only bits before 2 * k remain
		int hash = value % max_hash_size; // By formula in the article [V mod s(size of hash table)]
		loc[i - k + 1] =  point[hash]; // Tuple index has to be adjusted by k - 1
		point[hash] = i - k + 1; 
	}	
}

// Written by Marko Marfat
// Run Length encoding
// Turns a sequence like: 0 0 0 0 0 1 2 2 1 1 1 into a format 0-5 1-1 2-2 1-3
// Very useful for when theres a lot of repeating in a sequence
string RLE(){
	string output;
	int count = 1;
	output += to_string(t_seq_len[0]);
	int val_before = t_seq_len[0];
	for(int i = 1; i < t_seq_len.size(); i++){
		if (t_seq_len[i] == val_before){
			count++;
		}
		else{
			output += "-";
			output += to_string(count);
			output += " ";
			output += to_string(t_seq_len[i]);
			val_before = t_seq_len[i];
			count = 1;
		}
	}
	output += "-";
	output += to_string(count);
	output += "\n";

	return output;
}


// Written by Marko Marfat
// This function finds the best match for a given interval of the target sequence
string greedyMatching(){
	// Indexes for target tuple and mismatched subsequence
	int i = 0; 
	int p = 0; 
	string output;
	
	while(i < (t_final.length() - k + 1)){
		
		// The code up until the // line is doing the exact same thing as the initHT() function
		// but its doing it for the target sequence
		uint64_t value = 0; 
		for (int x = 0; x < k; x++){
			value = value << 2; 
			value = value + t_final[i + x]; 
		}
		
		int hash = value % max_hash_size;
		/////////////////////////////////////////////////////////////////////////////////////////////////
		
		int j = point[hash]; // Find the reference tuple with same hash as target tuple
		int pmax = -1; // Longest match in sequences
		int lmax = 0; // Length of longest match
		
		// Loop through all tuples in a hash map bucket
		while(j != -1){
			int l = 0; // Length of sequence
			
			// Check for matches
			// These are the conditions under which a match can be considered, they are pulled out for clarity
			// This isn't the "right" way to code, but it makes it clearer to see whats happening
			
			bool cond_1 = (j + l) < t_final.length() ; // Check if index is in bounds
			bool cond_2 = (j + l) < r_final.length(); // Check if index is in bounds
			bool cond_3 = r_final[j + l] == t_final[i + l]; // Actual char match checking
			
			// Loop and match
			while(cond_1 && cond_2 && cond_3){
				l++;  // Increase match length
				
				// Refresh conditions
				cond_1 = (j + l) < t_final.length() ;
				cond_2 = (j + l) < r_final.length(); 
				cond_3 = r_final[j + l] == t_final[i + l];
			}
			
			// Check for length of match and if its longer than last longest match, replace it
			if (l >= k && l > lmax){
				pmax = j;
				lmax = l;
			}
			
			// Next iteration (next tuple index)
			j = loc[j];
		}
		
		if (lmax > 0){
			int p_ = p;
			while (p_ <= i - 1){
				output += t_final[p_];
				
				p_++;
				if(p_ > (i - 1)){
					output += "\n";
				}
			}
			
			// This is the matched information
			output += to_string(pmax); // Match index
			output += ","; // Delimiter
			output += to_string(lmax); // Length of match
			output += "\n"; // Delimiter (\n since it will be stored in a file)
			
			p = i + lmax; // Update p
		}
		
		// Next iteration (sequence)
		i = i + lmax + 1;
	}
	
	// The remaining subsequence is a mismatched subsequence
	while (p < t_final.length()){
		output += t_final[p];
		p++;
	}
	
	return output;
}


// Written by Katarina Misura
//This function prints out instructions for the user
void user_interface(){
    cout<< "HiRGC v1.0" << endl;
    cout<< "Use: ./hirgc -r <reference> -t <target>"<< endl;
    cout<< "-r is the reference FASTA file -- required input" << endl;
    cout<< "-t is the target FASTA file -- required input" << endl;
    cout<< "Examples:\n";
    cout<< "hirgc -r testref.fa -t testtar.fa\n";
}

int main(int argc, char *argv[])
{	
	// Written by Katarina Misura
    	char *ref_file = NULL, *tar_file = NULL;
    	if(argc != 5){
        	cout<< "Wrong number of arguments" << endl;
        	user_interface();
        	return 0;
    	}
    	else{
        	if(strcmp(argv[1], "-r") == 0){
            		ref_file = argv[2];
        	}
        	else{
            		cout << "Input must contain -r in front of reference FASTA file." << argv[1] << endl;
            		user_interface();
            		return 0;
        	}
        	if(strcmp(argv[3], "-t") == 0){
            		tar_file = argv[4];
        	}
        	else{
            		cout << "Input must contain -t in front of target FASTA file. Your input was: "<<argv[3] << endl;
            		user_interface();
            		return 0;
        	}
    	}

	target_preprocess(tar_file);
	refrence_preprocess(ref_file);
		
	// Written by Marko Marfat
	struct  timeval  start;
	struct  timeval  end;
	unsigned long timer;
	gettimeofday(&start,NULL);
	
	ofstream myfile;
	myfile.open("output.fa");

	initHT(); // Initialize the Hash Table
	string rle = RLE(); // Conduct run length encoding
	string gm = greedyMatching(); // Conduct greedy matching
	
	myfile << id_tg << endl; // Save identifier to file
	myfile << rle;  	 // Save RLE encoded input to file
	saveDataToFile(myfile);  // Save auxilliary data to file
	myfile << gm;            // Save result of greedy matching to file
	myfile.close();
	
	// Compress output file using PPMd method from 7za library
	system("7za a -m0=PPMd output.7z output.fa");

	// Delete the output file (it's now saved in a compressed form)
	system("rm output.fa");

	gettimeofday(&end, NULL);
	timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	printf("total compression timer = %lf ms; %lf min\n", timer/1000.0, timer/1000.0/1000.0/60.0);
    
	return 0;
}
