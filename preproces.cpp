#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

//target sequence
string id_t;           // identifier from target FASTA file
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
                id_t = line;
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

int main()
{
    //target file preprocessing
    target_preprocess("testtar.txt");
    cout << "Encoded target sequence: " << t_final << endl;

    //reference file preprocessing
    refrence_preprocess("testref.txt");
    cout << "Encoded reference sequence: " << r_final << endl;
    
    return 0;
}