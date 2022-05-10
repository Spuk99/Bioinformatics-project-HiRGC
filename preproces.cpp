#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

//target sequence
string id_tg;           // identifier from target FASTA file
vector<int> t_seq_len; // vector of sequence lengths for each line
string t_seq_L = "";   // all sequences concatenated

vector<int> t_low_pos; // vector of positions of lowercase letters
vector<int> t_low_len; // vector of lengths of lowercase letters
string t_seq_L1 = "";  // all sequences concatenated with lowercase letters converted to uppercase

vector<int> t_N_pos;  // vector of positions of N letters in t_seq_L1
vector<int> t_N_len;  // vector of lengths of N letters in t_seq_L1
string t_seq_L2 = ""; // all sequences concatenated with N letters removed from t_seq_L1

vector<int> t_oth_pos; // vector of positions of other letters in t_seq_L2
vector<int> t_oth_len; // vector of lengths of other letters in t_seq_L2
string t_seq_L3 = "";  // all sequences concatenated with other letters removed from t_seq_L2

string t_final=""; // final sequence

// Reference sequence
string id_r;           // identifier from reference FASTA file
vector<int> r_seq_len; // vector of sequence lengths for each line
string r_seq_L = "";   // all sequences concatenated
string r_seq_L1 = "";  // all sequences concatenated with lowercase letters converted to uppercase
string r_seq_L3 = "";  // all sequences concatenated with other letters removed from t_seq_L2
string r_final=""; // final sequence

const int k = 20; // This is said to be the best value for k in the article
const int max_hash_size = 1 << 20; // Maximum length of the hash table
const int max_char_num = 1 << 28; // Maximum length of a chromosome

int *point = new int[max_hash_size]; // An array of entries 
int *loc = new int[max_char_num]; // An array of header pointers

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
                t_seq_L += line;
                t_seq_len.push_back(line.length());
            }

        }
        myfile.close();
    }
    else
    {
        cout << "Unable to open file";
    }  
    //preprocesing target sequence
    int interval; // interval 
    int last = -1; //flag for last position
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
        if (t_seq_L1[i] != 'N')
        {
            t_seq_L2 += t_seq_L1[i];
        }
    }
    if (last != -1)
    {
        interval = t_seq_L1.length() - last;
        t_N_len.push_back(interval);
    }

    last = -1;
    for (int i = 0; i < t_seq_L2.length(); i++)
    {
        if (!(t_seq_L2[i] == 'A' || t_seq_L2[i] == 'C' || t_seq_L2[i] == 'G' || t_seq_L2[i] == 'T') && last == -1)
        {
            last = i;
            t_oth_pos.push_back(last);
        }
        else if (last != -1 && (t_seq_L2[i] == 'A' || t_seq_L2[i] == 'C' || t_seq_L2[i] == 'G' || t_seq_L2[i] == 'T'))
        {
            interval = i - last;
            last = -1;
            t_oth_len.push_back(interval);
        }
        if ((t_seq_L2[i] == 'A' || t_seq_L2[i] == 'C' || t_seq_L2[i] == 'G' || t_seq_L2[i] == 'T'))
        {
            t_seq_L3 += t_seq_L2[i];
        }
    }
    if (last != -1)
    {
        interval = t_seq_L2.length() - last;
        t_oth_len.push_back(interval);
    }

    //encoding
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
                r_seq_L += line;
                r_seq_len.push_back(line.length());
            }
        }
        myfile.close();
    }
    else
    {
        cout << "Unable to open file";
    }
    //preprocesing reference sequence
    for (int i = 0; i < r_seq_L.length(); i++)
    {
        r_seq_L1 += toupper(r_seq_L[i]);
    }

    for (int i = 0; i < r_seq_L1.length(); i++)
    {
        if ((r_seq_L1[i] == 'A' || r_seq_L1[i] == 'C' || r_seq_L1[i] == 'G' || r_seq_L1[i] == 'T'))
        {
            r_seq_L3 += r_seq_L1[i];
        }
    }
    //encoding
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


/* This function constructs a hash table from the tuple values
of a reference sequence.  */
void initHT(){
	// Initialize entries
    for (int i = 0; i < max_hash_size; i++){
        point[i] = -1;
    }

	uint64_t value = 0;
	
	// Only k values is needed, skip rest
	for (int i = 0; i < k - 1; i++){
		value = value << 2; // Not 100% sure about this one, I think it shifts the value so new char can be added (not really explained in the article)
		// cout << "Value (1): " << value << endl; debugging
		value = value + r_final[i]; // Add char to tuple
		// cout << "Value (2): " << value << endl; debugging
	}
	
	uint64_t remain = ((uint64_t) 1 << (2 * k)) - 1; // Used as a mask for which bits will remain (only those below 2 * k)
	for (int i = k - 1; i < r_final.length(); i++){
		value = value << 2;
		// cout << "Value (1): " << value << endl; debugging
		value = value + r_final[i];
		// cout << "Value (2): " << value << endl; debugging
		value = value & remain; // Only bits before 2 * k remain
		int hash = value % max_hash_size; // By formula in the article [V mod s(size of hash table)]
		loc[i - k + 1] =  point[hash]; // Tuple index has to be adjusted by k - 1
		point[hash] = i - k + 1; 
	}	

}


int main()
{
    //target file preprocessing
    target_preprocess("testtar.txt");
    cout << "Encoded target sequence: " << t_final << endl;

    //reference file preprocessing
    refrence_preprocess("testref.txt");
    cout << "Encoded reference sequence: " << r_final << endl;
	
	initHT();
    
    
    return 0;
}