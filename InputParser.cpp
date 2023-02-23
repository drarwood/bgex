/*
   BGEN Extraction Tool (bgex)
   InputParser.cpp
   Functions related to pasring command line
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#include "InputParser.h"
#include <string>
#include <iostream>

using namespace std;

void ShowOptions() {
    cout << endl;
    cout << "Usage:" << endl;
    cout << "  --bgens     -b  [bgen list file]"                       << endl;
    cout << "  --samples   -s  [bgen sample file]"                     << endl;
    cout << "  --variants  -v  [file of variants to extract]"          << endl;
    cout << "  --extract   -e  [file of samples to keep (optional)]"   << endl;
    cout << "  --min-info  -m  [min INFO for genotype extraction/use]" << endl;
    cout << "  --dosages   -d  [flag to extract dosages]"              << endl;
    cout << "  --probs     -p  [flag to extract probabilities]"        << endl;
    cout << "  --pscore    -g  [flag to extract polygenic score]"      << endl;
    cout << "  --info      -q  [flag to extract info score]"           << endl;
    cout << "  --out       -o  [file prefix for outputs]"              << endl;
    cout << endl;
}

void ParseCommand(int argc, char** argv, ARGS *a) {

    if (argc > 1) {

        // defaults
        a->minInfo          = 0.0;
        a->extractDosages   = false;
        a->extractProbs     = false;
        a->extractPScore    = false;
        a->extractInfoScore = false;
        a->ok               = true;

        for (int i = 1; i != argc; ++i) {
            if (string(argv[i]) == "--bgens" || string(argv[i]) == "-b") {
                if (i+1 <= argc-1 && string(argv[i+1]).substr(0,1) != "-") {
                    a->bgenListFileString = argv[i+1];
                }
                else {
                    cout << "Not entered --bgen (-b) argument" << endl << endl;
                    a->ok = false;
                }
            }
            else if (string(argv[i]) == "--samples" || string(argv[i]) == "-s") {
                if (i+1 <= argc-1 && string(argv[i+1]).substr(0,1) != "-") {
                    a->sampleFileString = argv[i+1];
                }
                else {
                    cout << "Not entered --samples (-s) argument" << endl << endl;
                    a->ok = false;
                }
            }
            else if (string(argv[i]) == "--variants" || string(argv[i]) == "-v") {
                if (i+1 <= argc-1 && string(argv[i+1]).substr(0,1) != "-") {
                    a->variantFileString = argv[i+1];
                }
                else {
                    cout << "Not entered --variants (-v) argument" << endl << endl;
                    a->ok = false;
                }
            }
            else if (string(argv[i]) == "--extract" || string(argv[i]) == "-e") {
                if (i+1 <= argc-1 && string(argv[i+1]).substr(0,1) != "-") {
                    a->subjectIncFileString = argv[i+1];
                }
                else {
                    cout << "Not entered --extract (-e) argument" << endl << endl;
                    a->ok = false;
                }
            }
            else if (string(argv[i]) == "--min-info" || string(argv[i]) == "-i") {
                if (i+1 <= argc-1 && string(argv[i+1]).substr(0,1) != "-") {
                    a->minInfo = stod(argv[i+1]);
                }
                else {
                    cout << "Not entered --min-info (-i) argument" << endl << endl;
                    a->ok = false;
                }
            }
            else if (string(argv[i]) == "--dosages" || string(argv[i]) == "-d") {
                    a->extractDosages = true;
            }
            else if (string(argv[i]) == "--probs" || string(argv[i]) == "-p") {
                    a->extractProbs = true;
            }
            else if (string(argv[i]) == "--pscore" || string(argv[i]) == "-g") {
                    a->extractPScore = true;
            }
            else if (string(argv[i]) == "--info" || string(argv[i]) == "-q") {
                    a->extractInfoScore = true;
            }
            else if (string(argv[i]) == "--out" || string(argv[i]) == "-o") {
                if (i+1 <= argc-1 && string(argv[i+1]).substr(0,1) != "-") {
                    a->outFileString = argv[i+1];
                }
                else {
                    cout << "Not entered --out (-o) argument" << endl << endl;
                    a->ok = false;
                }
            }
        }


        if (!a->extractDosages && !a->extractProbs && !a->extractPScore && !a->extractInfoScore) {
            cout << "You must select at least one of dosage, probability, polygenic score, info score extraction" << endl;
            a->ok = false;
        }

    }
    else {
        a->ok = false;
    }

}
