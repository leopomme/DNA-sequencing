#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

/*-------------------------------------------------------------

    struct DNA_reg
    
    DNA region

-------------------------------------------------------------*/

struct DNA_reg
{
    char label;        // Region label ('N' unknown, 'C' coded)
    int num;        // Region number
    int i0;            // Starting basepair in region
    int iN;            // Final basepair in region
};

/*-------------------------------------------------------------

    class DNA
    
    DNA sequence

-------------------------------------------------------------*/

/*    DNA definition */

class DNA
{
private:
    vector<char> seq;            // Base pairs represented by a character sequence
    
    int gid_num;                // Identifiers for the DNA sequence
    string ref_str;
    string name_str;

    vector<DNA_reg> N_regions;
    vector<DNA_reg> C_regions;
    
    double N_reg_count;            // Gap region count
    double C_reg_count;            // Coded region count
    double T_reg_count;            // Total region count
    
    int G_count;                // G character count
    int A_count;                // N character count
    int T_count;                // ...
    int C_count;
    int R_count;
    int Y_count;
    int M_count;
    int K_count;
    int S_count;
    int W_count;
    int H_count;
    int B_count;
    int V_count;
    int D_count;
    int N_count;
    int X_count;                // unexpected character
    
public:
    DNA();
    ~DNA();
    
    void set_gid(int gid_in) {gid_num = gid_in;}
    void set_ref(string ref_in) {ref_str = ref_in;}
    void set_name(string name_in) {name_str = name_in;}
    int get_gid() {return gid_num;}
    string get_ref() {return ref_str;}
    string get_name() {return name_str;}

    void push_back(char seq_ch);
    void update_stats();
    void find_seq_man(string seq_str);
    void find_seq_file(string seq_file_str);

    void print_seq(int i0, int iN);
    void print_gap_region(int gnum);
    void print_coded_region(int cnum);
    void print_bp_range(int i0, int iN);
    void print_stats();
    void print_help();
    int size();

};

/*    DNA constructors and destructors */

DNA::DNA()
{
    gid_num = 0;

    N_reg_count = 0;
    C_reg_count = 0;
    T_reg_count = 0;

    G_count = 0;
    A_count = 0;
    T_count = 0;
    C_count = 0;
    R_count = 0;
    Y_count = 0;
    M_count = 0;
    K_count = 0;
    S_count = 0;
    W_count = 0;
    H_count = 0;
    B_count = 0;
    V_count = 0;
    D_count = 0;
    N_count = 0;
    X_count = 0;
}

DNA::~DNA()
{

}

/*    DNA member functions */

void DNA::push_back(char seq_ch)
{
    seq.push_back(seq_ch);
}

void DNA::print_seq(int i0, int iN)
{
    for (int i = i0; i <= iN; i++)
        cout << seq[i];
    cout << endl;
}

void DNA::print_gap_region(int gnum)
{
    cout << "Selected sequence:" << endl;
    cout << "Base pair range: ("
         << N_regions[gnum].i0+1 << ","
         << N_regions[gnum].iN+1 << ")" << endl;
    cout << "Gap region number: " << gnum+1 << endl;
    cout << endl;
    cout << "Sequence:" << endl;
    print_seq(N_regions[gnum].i0, N_regions[gnum].iN);
    cout << endl;
}

void DNA::print_coded_region(int cnum)
{
    cout << "Selected sequence:" << endl;
    cout << "Base pair range: ("
         << C_regions[cnum].i0+1 << ","
         << C_regions[cnum].iN+1 << ")" << endl;
    cout << "Coded region number: " << cnum+1 << endl;
    cout << endl;
    cout << "Sequence:" << endl;
    print_seq(C_regions[cnum].i0, C_regions[cnum].iN);
    cout << endl;
}

void DNA::print_bp_range(int i0, int iN)
{
    cout << "Selected sequence:" << endl;
    cout << "Base pair range: ("
         << i0+1 << ","
         << iN+1 << ")" << endl;
    cout << endl;
    cout << "Sequence:" << endl;
    print_seq(i0, iN);
    cout << endl;
}

void DNA::update_stats()
{
    DNA_reg temp_region;
    
    N_reg_count = 0;
    C_reg_count = 0;
    T_reg_count = 0;
    
    G_count = 0;
    A_count = 0;
    T_count = 0;
    C_count = 0;
    R_count = 0;
    Y_count = 0;
    M_count = 0;
    K_count = 0;
    S_count = 0;
    W_count = 0;
    H_count = 0;
    B_count = 0;
    V_count = 0;
    D_count = 0;
    N_count = 0;
    X_count = 0;
    
    for (int i = 0; i < seq.size(); i++)        // run through entire sequence
    {

        if (seq[i] == 'N')                        // read non-coding region
        {
            N_reg_count++;
            N_count++;
            
            temp_region.label = 'N';
            temp_region.num = N_reg_count;
            temp_region.i0 = i;
            while (i+1 < seq.size() && seq[i+1] == 'N')
            {
                i++;
                N_count++;
            }
            temp_region.iN = i;
            
            N_regions.push_back(temp_region);
        }
        else                                    // read coding region
        {
            C_reg_count++;
            if (seq[i] == 'G')
                G_count++;
            else if (seq[i] == 'A')
                A_count++;
            else if (seq[i] == 'T')
                T_count++;
            else if (seq[i] == 'C')
                C_count++;
            else if (seq[i] == 'R')
                R_count++;
            else if (seq[i] == 'Y')
                Y_count++;
            else if (seq[i] == 'M')
                M_count++;
            else if (seq[i] == 'K')
                K_count++;
            else if (seq[i] == 'S')
                S_count++;
            else if (seq[i] == 'W')
                W_count++;
            else if (seq[i] == 'H')
                H_count++;
            else if (seq[i] == 'B')
                B_count++;
            else if (seq[i] == 'V')
                V_count++;
            else if (seq[i] == 'D')
                D_count++;
            else
                X_count++;
            
            temp_region.label = 'C';
            temp_region.num = C_reg_count;
            temp_region.i0 = i;
            while (i+1 < seq.size() && seq[i+1] != 'N')
            {
                i++;
                
                if (seq[i] == 'G')
                    G_count++;
                else if (seq[i] == 'A')
                    A_count++;
                else if (seq[i] == 'T')
                    T_count++;
                else if (seq[i] == 'C')
                    C_count++;
                else if (seq[i] == 'R')
                    R_count++;
                else if (seq[i] == 'Y')
                    Y_count++;
                else if (seq[i] == 'M')
                    M_count++;
                else if (seq[i] == 'K')
                    K_count++;
                else if (seq[i] == 'S')
                    S_count++;
                else if (seq[i] == 'W')
                    W_count++;
                else if (seq[i] == 'H')
                    H_count++;
                else if (seq[i] == 'B')
                    B_count++;
                else if (seq[i] == 'V')
                    V_count++;
                else if (seq[i] == 'D')
                    D_count++;
                else
                    X_count++;
            }
            temp_region.iN = i;
            
            C_regions.push_back(temp_region);
        }
    }
    T_reg_count = N_reg_count + C_reg_count;

    
}

void DNA::find_seq_man(string seq_str)
{
    vector<int> matches;
    
    if (seq_str.size() >= 10)
        {

        for (int i=0; i < seq.size(); i++)
        {
            bool match = true;
            for (int j=0; j < seq_str.length(); j++)
            {
                if (seq[i+j] == seq_str[j])
                {
                    match = match && true;
                }
                else
                {
                    match = false;
                    break;
                }
            }
        
            if (match)
                matches.push_back(i);
        }
    
        if (matches.size() > 0)
        {
            for (int i=0; i<matches.size(); i++)
            {
                cout << i+1 << ". Base pair range: ("
                     << matches[i]+1 << "," << matches[i]+1 + seq_str.size()-1 << ")"
                     << endl;
                print_seq(matches[i]-20, matches[i]-1);
                print_seq(matches[i], matches[i]+ seq_str.size()-1);
                print_seq(matches[i]+ seq_str.size()-1+1, matches[i]+ seq_str.size()-1+20);
                cout << endl;
            }
        }
        else
        {
            cout << "No matching sequence found." << endl;
            cout << endl;
        }
    }
    else
    {
        cout << "Sequences must have 10 or more nucleotides." << endl;
        cout << endl;
    }
    
}

void DNA::find_seq_file(string file_name)
{
    vector<char> dna_seq;

    ifstream fin;
    string line_str;
    char temp_ch;
    
    fin.open(file_name,ios_base::in);
    
    if (fin.is_open())
    {
        vector<int> matches;
        
        cout << "Loading " << file_name << "..." << endl;
        getline(fin, line_str);
    
        while(!fin.eof())
        {
            fin.get(temp_ch);
        
            if (temp_ch == '\n')            // skip new line characters
                ;
            else
                dna_seq.push_back(temp_ch);
        }
        fin.close();
        
        cout << "Successful loading of " << file_name << endl;
        cout << endl;

        if (dna_seq.size() >= 10)
        {
            for (int i=0; i < seq.size(); i++)
            {
                bool match = true;
                for (int j=0; j < dna_seq.size(); j++)
                {
                    if (seq[i+j] == dna_seq[j])
                    {
                        match = match && true;
                    }
                    else
                    {
                        match = false;
                        break;
                    }
                }
        
                if (match)
                    matches.push_back(i);
            }
        
            if (matches.size() > 0)
            {
                cout << matches.size() << " matching sequence(s) found." << endl;
                cout << endl;
        
                for (int i=0; i<matches.size(); i++)
                {
                    cout << i+1 << ". Base pair range: ("
                         << matches[i]+1 << "," << matches[i]+1 + dna_seq.size()-1 << ")"
                         << endl;
                    print_seq(matches[i]-20, matches[i]-1);
                    print_seq(matches[i], matches[i]+ dna_seq.size()-1);
                    print_seq(matches[i]+ dna_seq.size()-1+1, matches[i]+ dna_seq.size()-1+20);
                    cout << endl;
                }
            }
            else
            {
                cout << "No matching sequence found." << endl;
                cout << endl;
            }
        }
        else
        {
            cout << "Sequences must have 10 or more nucleotides." << endl;
            cout << endl;
        }
        
    }
    else
    {
        cout << "Could not open " << file_name << endl;
    }
}

void DNA::print_stats()
{
    cout << "Sequence identifiers:" << endl;
    cout << "Name: \t\t" << name_str << endl;
    cout << "GID:  \t\t" << gid_num << endl;
    cout << "REF:  \t\t" << ref_str << endl;
    cout << endl;
    cout << "Region characteristics:" << endl;
    cout << "# regions: \t" << T_reg_count << endl;
    cout << "# N regions: \t" << N_reg_count << endl;
    cout << "# C regions: \t" << C_reg_count << endl;
    cout << endl;
    cout << "Base pair characteristics:" << endl;
    cout << "# base pairs: \t" << seq.size() << endl;
    cout << "G:    \t\t" << G_count << endl;
    cout << "A:    \t\t" << A_count << endl;
    cout << "T:    \t\t" << T_count << endl;
    cout << "C:    \t\t" << C_count << endl;
    cout << "R:    \t\t" << R_count << endl;
    cout << "Y:    \t\t" << Y_count << endl;
    cout << "M:    \t\t" << M_count << endl;
    cout << "K:    \t\t" << K_count << endl;
    cout << "S:    \t\t" << S_count << endl;
    cout << "W:    \t\t" << W_count << endl;
    cout << "H:    \t\t" << H_count << endl;
    cout << "B:    \t\t" << B_count << endl;
    cout << "V:    \t\t" << V_count << endl;
    cout << "D:    \t\t" << D_count << endl;
    cout << "N:    \t\t" << N_count << endl;
    cout << "Unknown:\t" << X_count << endl;
    cout << endl;
}

void DNA::print_help()
{
    cout << "Code Base Description" << endl;
    cout << "G    Guanine" << endl;
    cout << "A    Adenine" << endl;
    cout << "T    Thymine (Uracil in RNA)" << endl;
    cout << "C    Cytosine" << endl;
    cout << "R    Purine (A or G)" << endl;
    cout << "Y    Pyrimidine (C or T or U)" << endl;
    cout << "M    Amino (A or C)" << endl;
    cout << "K    Ketone (G or T)" << endl;
    cout << "S    Strong interaction (C or G)" << endl;
    cout << "W    Weak interaction (A or T)" << endl;
    cout << "H    Not-G (A or C or T) H follows G in the alphabet" << endl;
    cout << "B    Not-A (C or G or T) B follows A in the alphabet" << endl;
    cout << "V    Not-T (not-U) (A or C or G) V follows U in the alphabet" << endl;
    cout << "D    Not-C (A or G or T) D follows C in the alphabet" << endl;
    cout << "N    Any (A or C or G or T)" << endl;
    cout << endl;
}

int DNA::size()
{
    return seq.size();
}

/*-------------------------------------------------------------

    class DNA DB
    
    DNA Database

-------------------------------------------------------------*/

class DNA_DB
{
private:
    vector<DNA> dna_db;                    // store DNA sequences - by gaps
    
public:
    DNA_DB();
    ~DNA_DB();
    
    bool load_db(string file_name);
    void clear();
    int size();
    void print_stats();
    
    DNA& operator[] (int x) { return dna_db[x]; }
};

DNA_DB::DNA_DB()
{

}

DNA_DB::~DNA_DB()
{

}


bool DNA_DB::load_db(string file_name)
{
    DNA dna_seq;

    ifstream fin;
    
    string line_str;
    string trash;
    string gid_str;
    int gid_num;
    string ref_str;
    string name_str;
    
    char temp_ch;
    
    fin.open(file_name,ios_base::in);
    
    if (fin.is_open())
    {
        cout << "Loading " << file_name << "..." << endl;
    
        getline(fin, line_str);
    
        stringstream sin(line_str);
        getline(sin,trash,'|');
        getline(sin,gid_str,'|');
        getline(sin,trash,'|');
        getline(sin,ref_str,'|');
        getline(sin,name_str,'\n');
    
        gid_num = stoi(gid_str);
    
        dna_seq.set_gid(gid_num);
        dna_seq.set_ref(ref_str);
        dna_seq.set_name(name_str);
    
        while(!fin.eof())
        {
            fin.get(temp_ch);
        
            if (temp_ch == '\n')            // skip new line characters
                ;
            else
                dna_seq.push_back(temp_ch);
        }
    
        dna_seq.update_stats();
        //dna_seq.print_seq();
        //dna_seq.print_stats();
        
        dna_db.push_back(dna_seq);
        
        fin.close();
        
        cout << "Successful loading of " << file_name << endl;

        return true;
    }
    else
    {
        cout << "Could not open " << file_name << endl;
    
        return false;
    }
}

void DNA_DB::clear()
{
    dna_db.clear();
}

int DNA_DB::size()
{
    return dna_db.size();
}

void DNA_DB::print_stats()
{

    cout << "The DNA Sequence Database holds " << dna_db.size() << " sequences." << endl;
    cout << endl;
    for (int i = 0; i < dna_db.size(); i++)
    {
        cout << "Sequence " << i+1 << ":" << endl;
        cout << "Name: \t\t" << dna_db[i].get_name() << endl;
        cout << "GID:  \t\t" << dna_db[i].get_gid() << endl;
        cout << "REF:  \t\t" << dna_db[i].get_ref() << endl;
        cout << "# base pairs:\t" << dna_db[i].size() << endl;
        cout << endl;
    }

}



/*-------------------------------------------------------------

    Functions

-------------------------------------------------------------*/

vector<string> menu_filenames();
char menu_db(vector<string> file_names);
char menu_dna();

int request_gap_region();
int request_coded_region();
pair<int, int> request_bp_range();
string request_seq_str();
string request_seq_file_str();

/*-------------------------------------------------------------

    main()

-------------------------------------------------------------*/

int main()
{
    DNA_DB dna_db;
    int db_i = 0;
    
    DNA seq1;
    
    string line_str;
    string trash;

    vector<string> file_names;
    
    bool m1_fns_valid = false;
    char option_ch = ' ';
    
    
    // OPENING SCREEN //
    
    cout << endl;
    cout << "DNA Sequence Database Software" << endl;
    cout << endl;
    
    // MENU 1 // ask user for data files
    
    while (!m1_fns_valid)
    {
        file_names = menu_filenames();
        
        if (file_names.size() > 0)
            m1_fns_valid = true;
            
        for (int i = 0; i < file_names.size(); i++)
        {
            bool state;
            state = dna_db.load_db(file_names[i]);
            m1_fns_valid = m1_fns_valid && state;
        }
        
        if (!m1_fns_valid)
        {
            file_names.clear();
            dna_db.clear();
        }

        cout << endl;
    }
    
    
    // MENU 2 //
    
    while (!(option_ch == 'Q' || option_ch == 'q'))
    {
        option_ch = menu_db(file_names);
            
        if (option_ch == 'S' || option_ch == 's')
            dna_db.print_stats();
        else
        {
            for (int i = 0; i < dna_db.size(); i++)
            {
                if (option_ch == (i+'1'))
                {
                    db_i = option_ch - '1';
                    while (!(option_ch == 'Q' || option_ch == 'q' || option_ch == 'R' || option_ch == 'r'))
                    {
                        option_ch = menu_dna();

                        if (option_ch == 'H' || option_ch == 'h')
                            dna_db[db_i].print_help();
                        else if (option_ch == 'S' || option_ch == 's')
                            dna_db[db_i].print_stats();
                        
                        //
                        else if (option_ch == '1')
                        {
                            // analyse gap region
                            int gap_num;
                            gap_num = request_gap_region();
                            dna_db[db_i].print_gap_region(gap_num-1);
                        }
                        else if (option_ch == '2')
                        {
                            // analyse coded region
                            int coded_num;
                            coded_num = request_coded_region();
                            dna_db[db_i].print_coded_region(coded_num-1);
                        }
                        else if (option_ch == '3')
                        {
                            // analyse base pair range
                            pair<int, int> rng;
                            rng = request_bp_range();
                            dna_db[db_i].print_bp_range(rng.first-1, rng.second-1);
                         }
                        else if (option_ch == '4')
                        {
                            // find DNA sequence by manual input
                            string seq_str;
                            seq_str = request_seq_str();
                            dna_db[db_i].find_seq_man(seq_str);
                            
                        }
                        else if (option_ch == '5')
                        {
                            // find DNA sequence by file input
                            string seq_file_str;
                            seq_file_str = request_seq_file_str();
                            dna_db[db_i].find_seq_file(seq_file_str);
                        }
                    }
                }
            }
        }
    }
    
    cout << "Program ended." << endl;

    return 0;
}

vector<string> menu_filenames()
{
    vector<string> file_names;
    string line_str;
    string temp_str;

    // display menu text
    cout << "Specify the name of DNA sequence file names you would like to load. For multiple files, add a ',' between each file name." << endl;
    cout << ">";
            
    // receive file names
    getline(cin, line_str);
    cout << endl;
    
    // parse file names
    stringstream sin(line_str);
    while (sin.good())
    {
        getline(sin, temp_str, ',');
        
        stringstream temp_sin(temp_str);
        temp_sin >> temp_str;
        
        file_names.push_back(temp_str);
    }
    
    return file_names;
}

char menu_db(vector<string> file_names)
{
    string line_str;
    char option_ch;

    cout << "Select one of the following options:" << endl;
    cout << "(S) Summary statistics of the DNA database" << endl;
    for (int i = 0; i < file_names.size(); i++)
        cout << "(" << i+1 << ") Analyse " << file_names[i] << endl;
    cout << "(Q) Quit" << endl;
    cout << ">";

    getline(cin, line_str);
    
    stringstream sin(line_str);
    
    sin >> option_ch;

    cout << endl;
    
    return option_ch;
}

char menu_dna()
{
    string line_str;
    char option_ch;

    cout << "Select one of the following options" << endl;
    cout << "(H) Help" << endl;
    cout << "(S) Summary statistics of the DNA sequence" << endl;
    cout << "(1) Analyse gap region" << endl;
    cout << "(2) Analyse coded region" << endl;
    cout << "(3) Analyse base pair range" << endl;
    cout << "(4) Find DNA sequence by manual input" << endl;
    cout << "(5) Find DNA sequence by file input" << endl;
    cout << "(R) Return to the previous menu" << endl;
    cout << "(Q) Quit" << endl;
    cout << ">";
    getline(cin, line_str);
    
    stringstream sin(line_str);
    sin >> option_ch;
    cout << endl;
    
    return option_ch;
}

int request_gap_region()
{
    string line_str;
    int option_int;

    cout << "Enter gap region number:" << endl;
    cout << ">";
    getline(cin, line_str);
    
    stringstream sin(line_str);
    sin >> option_int;
    cout << endl;
    
    return option_int;
}

int request_coded_region()
{
    string line_str;
    int option_int;

    cout << "Enter coded region number:" << endl;
    cout << ">";
    getline(cin, line_str);
    
    stringstream sin(line_str);
    sin >> option_int;
    cout << endl;
    
    return option_int;
}

pair<int, int> request_bp_range()
{
    string line_str;
    string temp_str;
    pair<int,int> rng_int;

    cout << "Enter a comma ',' separated base pair range:" << endl;
    cout << ">";
    getline(cin, line_str);
    
    stringstream sin(line_str);
    getline(sin, temp_str, ',');
            
    stringstream temp1_sin(temp_str);
    temp1_sin >> rng_int.first;
    
    sin >> rng_int.second;
    cout << endl;
    
    return rng_int;
}


string request_seq_str()
{
    string line_str;
    string temp_str;

    // display menu text
    cout << "Specify the DNA sequence nucleotides you would like to find:" << endl;
    cout << ">";
            
    // receive file names
    getline(cin, line_str);
    cout << endl;
    
    return line_str;
}

string request_seq_file_str()
{
    string line_str;
    string temp_str;

    // display menu text
    cout << "Specify the DNA sequence file you would like to find:" << endl;
    cout << ">";
            
    // receive file names
    getline(cin, line_str);
    cout << endl;
    
    return line_str;
}



