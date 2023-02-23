/*
   BGEN Extraction Tool (bgex)
   FileReader.cpp
   Functions related to reading human readable files
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#include "FileWriter.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

void OutputProbs(map<int,SAMPLE>& s, map<string,VARIANT>& v, vector<string>& o, string& f) {
    /* This function outputs the genotype probabilities for each individual
       Parameters:
         s: ref to map of sample file index -> SAMPLE struct
         v: ref to map of user defined chr:pos:a1:a2 -> VARIANT struct
         o: ref to vector containing order of variant extraction based on user defined chr:pos:a1:a2
         f: ref to string holding file-prefix for output
       Returns:
         NA
    */
    ofstream outFile;
    outFile.open(f+".probs", ios_base::out);
    if (outFile.is_open()) {
        // generate the header
        outFile << "n_eid";
        for (size_t i = 0; i < o.size(); ++i) {
            if (v[o[i]].use) {
                outFile << "\t" << v[o[i]].chr+":"+v[o[i]].pos+":"+v[o[i]].bgen_a1+":"+v[o[i]].bgen_a2;
            }
        }
        outFile << endl;
        // now cycle through the samples and [robs
        for (map<int,SAMPLE>::iterator it=s.begin(); it!=s.end(); ++it) {
            // check this is an individual wanted within data extraction
            if (it->second.use) {
                // output n_eid
                outFile << it->second.id;
                // cycle through prob vector, 2 at a time, converting the bytes to probs
                int var_index = 0;
                for (size_t i = 0; i < it->second.prob_bytes.size(); i=i+2) {
                    // based on the ordering of the probs, check whether the prob pair need to be extracted and output
                    if (v[o[var_index]].use) {
                        outFile << "\t" << it->second.prob_bytes[i]  /(pow(2.0,8.0)-1)
                                << ","  << it->second.prob_bytes[i+1]/(pow(2.0,8.0)-1);
                    }
                    var_index++;
                }
                outFile << endl;
            }
        }
        outFile.close();
    }
    else {
        cout << "Can not open " << f << " for writing" << endl;
        exit(EXIT_FAILURE);
    }
}


void OutputDosages(map<int,SAMPLE>& s, map<string,VARIANT>& v, vector<string>& o, string& f) {
    /* This function outputs the genotype dosages for each individual (aligned to bgen A2)
       Parameters:
         s: ref to map of sample file index -> SAMPLE struct
         v: ref to map of user defined chr:pos:a1:a2 -> VARIANT struct 
         o: ref to vector containing order of variant extraction based on user defined chr:pos:a1:a2
         f: ref to string holding file-prefix for output
       Returns:
         NA
    */
    ofstream outFile;
    outFile.open(f+".dosages", ios_base::out);
    if (outFile.is_open()) {
        // generate the header
        outFile << "n_eid";
        for (size_t i = 0; i < o.size(); ++i) {
            if (v[o[i]].use) {
                outFile << "\t" << v[o[i]].chr+":"+v[o[i]].pos+":"+v[o[i]].bgen_a1+":"+v[o[i]].bgen_a2;
            }
	}
	outFile << endl;
        // now cycle through the samples and [robs
        for (map<int,SAMPLE>::iterator it=s.begin(); it!=s.end(); ++it) {
            // check this is an individual wanted within data extraction
            if (it->second.use) {
                // output n_eid
                outFile << it->second.id;
                // cycle through prob vector, 2 at a time, converting the bytes to probs
                int var_index = 0;
                for (size_t i = 0; i < it->second.prob_bytes.size(); i=i+2) {
                    // based on the ordering of the probs, check whether the prob pair need to be extracted and dosage output
                    if (v[o[var_index]].use) {
                        double pab = it->second.prob_bytes[i+1]/(pow(2.0,8.0)-1);
                        double pbb = (1-  (it->second.prob_bytes[i]/(pow(2.0,8.0)-1) + it->second.prob_bytes[i+1]/(pow(2.0,8.0)-1)) );
                        outFile << "\t" << pab+(2*pbb);
                    }
                    var_index++;
                }
                outFile << endl;
            }
	}
	outFile.close();
    }
    else {
	cout << "Can not open " << f << " for writing" << endl;
        exit(EXIT_FAILURE);
    }
}


void OutputPS(map<int,SAMPLE>& s, string& f) {
    /* This function outputs polygenic scores for each individual
       Parameters:
         s: ref to map of sample file index -> SAMPLE struct
         f: ref to string holding file-prefix for output
       Returns:
         NA
    */
    ofstream outFile;
    outFile.open(f+".pscores", ios_base::out);
    if (outFile.is_open()) {
        // generate the header
        outFile << "n_eid\tpscore\n";
        // now cycle through the samples and [robs
        for (map<int,SAMPLE>::iterator it=s.begin(); it!=s.end(); ++it) {
            // check this is an individual wanted within data extraction
            if (it->second.use) {
                // output n_eid
                outFile << it->second.id << "\t" << it->second.ps << endl;
            }
	}
	outFile.close();
    }
}


void OutputInfoScores(map<string,VARIANT>& v, string& f) {
    /* This function outputs info scores calculated for each variant
       Parameters:
         s: ref to map of user defeined chr:pos:a1:a2 -> VARIANT struct
         f: ref to string holding file-prefix for output
       Returns:
         NA
    */
    ofstream outFile;
    outFile.open(f+".infoscores", ios_base::out);
    if (outFile.is_open()) {
        // generate the header
        outFile << "variant\tinfo_score\n";
        // now cycle through the samples and [robs
        for (map<string,VARIANT>::iterator it=v.begin(); it!=v.end(); ++it) {
            // check this is a variant passinging min info threshold
            if (it->second.use) {
                // output id and info score
                outFile << it->first << "\t" << it->second.info_score << endl;
            }
        }
	outFile.close();
    }

}

