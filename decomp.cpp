#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <tuple>
#include <sys/time.h>

using namespace std;

// Reference sequence
string id_r;           // identifier from reference FASTA file
vector<int> r_seq_len; // vector of sequence lengths for each line
string r_seq_L = "";   // all sequences concatenated
string r_seq_L1 = "";  // all sequences concatenated with lowercase letters converted to uppercase
string r_final=""; // final reference sequence

//Target sequence
string id_tg; // identifier from target FASTA file
string t_rle =""; // RLE of target sequence
vector<string> t_rle_seq; // vector of RLE sequences
vector<int> t_rle_len; // vector of sequence lengths for each line
vector<int> t_rle_pos;  // all positions of RLE
string t_low = "";  // lower case auxiliary information
vector<string> t_low_seq; // vector of lowercase letters sequences
vector<int> t_low_pos; // vector of positions of lowercase letters
vector<int> t_low_len; // vector of lengths of lowercase letters
string t_N=""; // N's auxiliary information
vector<string> t_N_seq; // vector of N sequences
vector<int> t_N_pos;  // vector of positions of N letters in t_seq_L1
vector<int> t_N_len;  // vector of lengths of N letters in t_seq_L1
string t_oth=""; // other auxiliary information
vector<string> t_oth_seq; // vector of other sequences
vector<int> t_oth_pos; // vector of positions of other letters in t_seq_L2
vector<string> t_oth_ch; // vector of characters of other letters in t_seq_L2
vector<int> t_oth_len; // vector of lengths of other letters in t_seq_L2

vector<int> match_pos; // vector of positions of matches
vector<int> match_len; // vector of lengths of matches
vector<string> miss_ch; // vector of characters of missmatches

string t_final=""; // final sequence


// Written by Katarina Misura
// Get sequence from reference FASTA file
void refrence_get_seq(string file_name){
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
        cout << "Unable to open file"<<endl;
        //add exit from program
        throw runtime_error("Closing program!");
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
            r_final += r_seq_L1[i];
        }
    }
}

// Written by Katarina Misura
//splits strings into vectors of strings - first split that is required to get seperate values
vector<string> split_string(string str, char delimiter){
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while(getline(ss, tok, delimiter)) { //Extracts characters from ss and stores them into tok until the delimitation character delimiter is found 
        internal.push_back(tok);
    }

    return internal;
}

// Written by Katarina Misura
//Splits string into two integers, one represents positions, other represents length of interval (or length of line - RLE)
tuple<vector<int>,vector<int>> split_extract_pos_len(vector<string> input, char delimiter){
    vector<int> pos;
    vector<int> len;
    for (int i = 0; i < input.size(); i++)
    {
        vector<string> temp = split_string(input[i], delimiter);
        pos.push_back(stoi(temp[0]));
        len.push_back(stoi(temp[1]));
    }
    return make_tuple(pos,len);
}

// Written by Katarina Misura
//splits string into integer and string (example: 21-PSDDKL => 21, PSDDKL)
tuple<vector<int>,vector<string>> split_extract_pos_str(vector<string> input, char delimiter){
    vector<int> pos;
    vector<string> str;
    for (int i = 0; i < input.size(); i++)
    {
        vector<string> temp = split_string(input[i], delimiter);
        pos.push_back(stoi(temp[0]));
        str.push_back(temp[1]);
    }
    return make_tuple(pos,str);
}

// Written by Katarina Misura
//read target compressed file and extract RLE, lowercase, N, other information, matches and missmatches
void read_target_comp(string file_name){
    string line;
    ifstream myfile(file_name);
    if (myfile.is_open())
    {
        getline(myfile, line);
        id_tg = line; // identifier from target FASTA file
        getline(myfile, line);
        t_rle = line; // RLE of target sequence
        getline(myfile, line);
        t_low = line; // lower case auxiliary information
        getline(myfile, line);
        t_N = line; // N's auxiliary information
        getline(myfile, line);
        t_oth = line; // other auxiliary information

        //read matches and missmatches and store them in vectors
        while (getline(myfile, line))
        {
            vector<string> temp = split_string(line, ',');
            if(temp.size()==2){
                match_pos.push_back(stoi(temp[0]));
                match_len.push_back(stoi(temp[1]));
            }else{
                miss_ch.push_back(temp[0]);
            }
        }
        myfile.close();
    }
    else
    {
        cout << "Unable to open file"<<endl;
        //add exit from program 
        throw runtime_error("Closing program!");
    }

    //split RLE
    t_rle_seq = split_string(t_rle, ' ');
    //split lower case
    t_low_seq = split_string(t_low, ' ');
    //split N's
    t_N_seq = split_string(t_N, ' ');
    //split other
    t_oth_seq = split_string(t_oth, ' ');
    
    //split RLE sequence and extract position and length of RLE
    tie(t_rle_pos,t_rle_len) = split_extract_pos_len(t_rle_seq, '-');
    //split lower case sequence and extract position and length of lower case
    tie(t_low_pos,t_low_len) = split_extract_pos_len(t_low_seq, '-');
    //split N's sequence and extract position and length of N's
    tie(t_N_pos,t_N_len) = split_extract_pos_len(t_N_seq, '-');
    //split other sequence and extract position and character of other
    tie(t_oth_pos,t_oth_ch) = split_extract_pos_str(t_oth_seq, '-');
}




// Written by Katarina Misura
void user_interface(){
    cout<< "HiRGC decompression v1.0" << endl;
    cout<< "Use: ./decomp -r <reference> -t <target>"<< endl;
    cout<< "-r is the reference FASTA file -- required input" << endl;
    cout<< "-t is the target sequence .7z format file -- required input" << endl;
    cout<< "Examples:\n";
    cout<< "decomp -r testref.fa -t output.fa.7z\n";
}

int main(int argc, char *argv[]){
    //Written by Katarina Misura
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
            cout << "Input must contain -t in front of target .7z file. Your input was: "<<argv[3] << endl;
            user_interface();
            return 0;
        }
    }
    
    struct  timeval  start;
	struct  timeval  end;
	unsigned long timer;
	gettimeofday(&start,NULL);
    
    refrence_get_seq(ref_file);
    /*
    //decompress the target file -- nez jel ovo uopce dobro???
    string command = "7z e " + string(tar_file);
    system(command.c_str());
    //get the name of the target file
    string tar_file_name = tar_file;
    tar_file_name.erase(tar_file_name.length()-3);
    */

    //proces the target file
    read_target_comp(tar_file);

    //testing of currnet code:
    cout << "RLE pos and length: " << endl;
    for(int i = 0; i < t_rle_pos.size(); i++){
        cout << t_rle_pos[i] << " " << t_rle_len[i] << "   ";
    }
    cout << endl;
    cout << "Low pos and length: " << endl;
    for(int i = 0; i < t_low_pos.size(); i++){
        cout << t_low_pos[i] << " " << t_low_len[i] << "   ";
    }
    cout << endl;
    cout << "N pos and length: " << endl;
    for(int i = 0; i < t_N_pos.size(); i++){
        cout << t_N_pos[i] << " " << t_N_len[i] << "   ";
    }
    cout << endl;
    cout << "Other pos and character: " << endl;
    for(int i = 0; i < t_oth_pos.size(); i++){
        cout << t_oth_pos[i] << " " << t_oth_ch[i] << "   ";
    }
    cout << endl;
    cout << "Match pos and length: " << endl;
    for(int i = 0; i < match_pos.size(); i++){
        cout << match_pos[i] << " " << match_len[i] << "   ";
    }
    cout << endl;
    cout << "Miss match: " << endl;
    for(int i = 0; i < miss_ch.size(); i++){
        cout << miss_ch[i] << "   ";
    }
    cout << endl;

    
    gettimeofday(&end, NULL);
	timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	printf("total compression timer = %lf ms; %lf min\n", timer/1000.0, timer/1000.0/1000.0/60.0);
    
    return 0;
}