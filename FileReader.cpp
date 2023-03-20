/*
   BGEN Extraction Tool (bgex)
   FileReader.cpp
   Functions related to reading human readable files
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#include "FileReader.h"
#include "StringOperations.h"

using namespace std;

void ReadBGENListFile(string& f, map<string,string>& b) {
    /*
    Function to read list of bgen and samples files to maps: chr->bgen, chr->sample
    Parameters:
      f : ref to string variable holding chromosome bgen list filename
      b : ref to map to fill with chromosome and corresponding bgen filenames
    Returns:
      na
    */
    string l, chr, bgen;
    ifstream bgenListFile(f.c_str());
    if (bgenListFile.is_open()) {
        while (bgenListFile.good()) {
            getline (bgenListFile, l);
            if (!l.empty()) {
                istringstream iss;
                iss.str(l);
                iss >> chr;
                iss >> bgen;
                b[chr] = bgen;
            }
        }
        bgenListFile.close();
    }
    else {
        cout << " error" << endl;
        cout << "Unable to open " << f << " for reading" << endl;
        exit(EXIT_FAILURE);
    }
}


void ReadBGENSampleFile(string& f, map<int, SAMPLE>& m) {
    /*
    Function to read IDs from bgen sample list to vector
    Parameters:
      f : ref to string variable holding bgen sample filename
      v : ref to vector to fill with sample ids - fid:iid format
    Returns:
      na
    */
    string l, iid, fid;
    int index = 0;
    ifstream sampleFile(f.c_str());
    if (sampleFile.is_open()) {
        // skip the first 2 lines
        getline (sampleFile, l);
        getline (sampleFile, l);
        // read remaining lines
        while (sampleFile.good()) {
            getline (sampleFile, l);
            if (!l.empty()) {
                istringstream iss;
                iss.str(l);
                iss >> fid;
                iss >> iid;
                SAMPLE thisSample;
                thisSample.id = fid+":"+iid;
                m[index] = thisSample;
                index++;
            }
        }
        sampleFile.close();
    }
    else {
        cout << " error" << endl;
        cout << "Unable to open " << f << " for reading" << endl;
        exit(EXIT_FAILURE);
    }
}


void ReadVariantListFile(string& f, map<string,VARIANT>& m, vector<string>& c) {
    /*
    Function to read variants to extract/use into VARIANT struct vector held within a map of chromosomes
    Parameters:
      f : ref to string variable holding variant list filename. Assumes no header
      m : ref to map of user defined chr:pos:a1:a2 -> VARIANT struct
      c : ref to vector of strings to populate with chromosome cods
    Returns:
      NA
    */
    string l, chr, pos, a1, a2;
    vector<string> d;
    ifstream variantFile(f.c_str());
    if (variantFile.is_open()) {
        while (variantFile.good()) {
            bool pos_num = true;
            getline (variantFile, l);
            if (!l.empty()) {
                split(l, "\t", d);
                // check pos is a number
                for (string::iterator it = d[1].begin(); it != d[1].end(); ++it) {
                    if (!isdigit(*it)) {
                        pos_num = false;
                        break;
                    }
                }
                // if number carry on
                if (pos_num) {
                    VARIANT thisVar;
                    thisVar.chr = d[0];
                    thisVar.pos = d[1];
                    thisVar.a1  = d[2];
                    thisVar.a2  = d[3];
                    thisVar.a1_weight = 0.0;
                    if (d.size() == 5 && d[4] != "") {
                        thisVar.a1_weight = stringToDouble(d[4]);
                    }
                    m[d[0]+":"+d[1]+":"+d[2]+":"+d[3]]=thisVar;
                    if (!count(c.begin(), c.end(), d[0])) {
                        c.push_back(d[0]);
                    }
                }
                // if not ignore and print user message
                else {
                    cout << l << " ignored as position is not a number" << endl;
                }
                d.clear();
            }
	}
	variantFile.close();
    }
    else {
	cout << " error" << endl;
        cout << "Unable to open " << f << " for reading" << endl;
        exit(EXIT_FAILURE);
    }
}


void ReadSubjectIncListFile(string& f, vector<string>& v) {
    /*
    Function to read IDs from bgen sample list to vector
    Parameters:
      f : ref to string variable holding sample inclusion list filename
      v : ref to vector to fill with sample ids - fid:iid format
    Returns:
      NA
    */
    string l, iid, fid;
    ifstream sampleFile(f.c_str());
    if (sampleFile.is_open()) {
        while (sampleFile.good()) {
            getline (sampleFile, l);
            if (!l.empty()) {
                istringstream iss;
                iss.str(l);
                iss >> fid;
                iss >> iid;
                v.push_back(fid+":"+iid);
            }
        }
        sampleFile.close();
    }
    else {
        cout << " error" << endl;
        cout << "Unable to open " << f << " for reading" << endl;
        exit(EXIT_FAILURE);
    }
}


void FlagSubjectsforExclusion(vector<string>& v, map<int, SAMPLE>& s) {

    /* This method takes the vector of subjects to include if provided
       and sets use=false for any subject not in the vector held in the sample map
       Parameters:
         v: ref to vector of strings holding sample IDs to include in outputs
         s: ref to map holding sample index -> SAMPPLE STRUCT
       Returns:
         NA
    */

    for (map<int, SAMPLE>::iterator it = s.begin(); it != s.end(); ++it) {
        if (find(v.begin(), v.end(), it->second.id) == v.end()) {
            it->second.use = false;
        }
    }

}
