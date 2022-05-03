#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

string id;           // identifier from FASTA file
vector<int> seq_len; // vector of sequence lengths for each line
string seq_L = "";   // all sequences concatenated

vector<int> low_pos; // vector of positions of lowercase letters
vector<int> low_len; // vector of lengths of lowercase letters
string seq_L1 = "";  // all sequences concatenated with lowercase letters converted to uppercase

vector<int> N_pos;  // vector of positions of N letters in seq_L1
vector<int> N_len;  // vector of lengths of N letters in seq_L1
string seq_L2 = ""; // all sequences concatenated with N letters removed from seq_L1

vector<int> oth_pos; // vector of positions of other letters in seq_L2
vector<int> oth_len; // vector of lengths of other letters in seq_L2
string seq_L3 = "";  // all sequences concatenated with other letters removed from seq_L2

int main()
{

    for (string line; getline(cin, line);)
    {
        if (line[0] == '>')
        {
            id = line;
        }
        else
        {
            seq_L += line;
            seq_len.push_back(line.length());
        }
    }
    int interval; // interval lowercase letters
    int last = -1;
    for (int i = 0; i < seq_L.length(); i++)
    {
        if (seq_L[i] >= 'a' && seq_L[i] <= 'z' && last == -1)
        {
            last = i;
            low_pos.push_back(last);
        }
        else if (last != -1 && !(seq_L[i] >= 'a' && seq_L[i] <= 'z'))
        {
            interval = i - last;
            last = -1;
            low_len.push_back(interval);
        }
        seq_L1 += toupper(seq_L[i]);
    }
    if (last != -1)
    {
        interval = seq_L.length() - last;
        low_len.push_back(interval);
    }

    last = -1;
    for (int i = 0; i < seq_L1.length(); i++)
    {
        if (seq_L1[i] == 'N' && last == -1)
        {
            last = i;
            N_pos.push_back(last);
        }
        else if (last != -1 && seq_L1[i] != 'N')
        {
            interval = i - last;
            last = -1;
            N_len.push_back(interval);
        }
        if (seq_L1[i] != 'N')
        {
            seq_L2 += seq_L1[i];
        }
    }
    if (last != -1)
    {
        interval = seq_L1.length() - last;
        N_len.push_back(interval);
    }

    last = -1;
    for (int i = 0; i < seq_L2.length(); i++)
    {
        if (!(seq_L2[i] == 'A' || seq_L2[i] == 'C' || seq_L2[i] == 'G' || seq_L2[i] == 'T') && last == -1)
        {
            last = i;
            oth_pos.push_back(last);
        }
        else if (last != -1 && (seq_L2[i] == 'A' || seq_L2[i] == 'C' || seq_L2[i] == 'G' || seq_L2[i] == 'T'))
        {
            interval = i - last;
            last = -1;
            oth_len.push_back(interval);
        }
        if ((seq_L2[i] == 'A' || seq_L2[i] == 'C' || seq_L2[i] == 'G' || seq_L2[i] == 'T'))
        {
            seq_L3 += seq_L2[i];
        }
    }
    if (last != -1)
    {
        interval = seq_L2.length() - last;
        oth_len.push_back(interval);
    }

    cout << id << endl;
    cout << "Sequence length: " << seq_L3.length() << endl;
    cout << "Sequemce: " << seq_L3 << endl;
    return 0;
}