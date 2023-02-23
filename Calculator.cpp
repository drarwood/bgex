/*
   BGEN Extraction Tool (bgex)
   Calculator.cpp
   Functions related to calculations involving genotypes
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#include "FileReader.h"
#include "Calculator.h"
#include <math.h>


using namespace std;

void CalculatePS(map<int,SAMPLE>& s, map<string,VARIANT>& v, vector<string>& o) {

    /*
    This function calculates and assigns polygenic predictor to each individual
    Parameters:
      s : ref to map holdibg bgen sample index -> SAMPLE struct
      v : ref to map holiding user defined chr:pos:a1:a2 -> VARIANT struct
      o : ref to vector holding respective order variants were read and therefore defines repsecitive order of probs
    */

    // cycle through individuals
    for (map<int,SAMPLE>::iterator it=s.begin(); it!=s.end(); ++it) {

        // check this is an individual wanted within data extraction
        if (it->second.use) {

           // initialise ps
           it->second.ps = 0.0;

            // cycle through prob vector, 2 at a time, converting the bytes to probs
            int var_index = 0;
            for (size_t i = 0; i < it->second.prob_bytes.size(); i=i+2) {

                // based on the ordering of the probs and respetive variant, check whether the prob pair need to be extracted and output
                if (v[o[var_index]].use && v[o[var_index]].a1_weight != 0.0) {

                    // extract the probabilities
                    double pab = it->second.prob_bytes[i+1]/(pow(2.0,8.0)-1);
                    double pbb = (1-  (it->second.prob_bytes[i]/(pow(2.0,8.0)-1) + it->second.prob_bytes[i+1]/(pow(2.0,8.0)-1)) );
                    double dosage = pab+(2*pbb);

                    // check whether user defined a1 = bgen_a2 - we will
                    if (v[o[var_index]].a1_weight > 0 && v[o[var_index]].a1 == v[o[var_index]].bgen_a2)  {
                        it->second.ps += dosage*v[o[var_index]].a1_weight;
                    }
                    else if (v[o[var_index]].a1_weight > 0 && v[o[var_index]].a1 == v[o[var_index]].bgen_a1) {
                        it->second.ps += (2-dosage)*v[o[var_index]].a1_weight;
                    }
                    else if (v[o[var_index]].a1_weight < 0 && v[o[var_index]].a1 == v[o[var_index]].bgen_a1) {
                         it->second.ps += dosage*abs(v[o[var_index]].a1_weight);
                    }
                    else if (v[o[var_index]].a1_weight < 0 && v[o[var_index]].a1 == v[o[var_index]].bgen_a2) {
                         it->second.ps += (2-dosage)*abs(v[o[var_index]].a1_weight);
                    }
                    else {
                         cerr << "oops: something unexpected has been encountered - will not include in predictor" << endl;
                    }

                } // ignore variant if not to be used in calculation

                var_index++;
            }
        } // ignore individual if not to be used
    }
}

