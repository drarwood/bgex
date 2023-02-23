/*
   BGEN Extraction Tool (bgex)
   InputParser.h
   Functions related to pasring command line
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#ifndef __INPUTPARSER__
#define __INPUTPARSER__

#include <string>

struct ARGS
{
    std::string bgenListFileString;
    std::string sampleFileString;
    std::string variantFileString;
    std::string subjectIncFileString;
    double      minInfo;
    bool        extractDosages;
    bool        extractProbs;
    bool        extractPScore;
    bool        extractInfoScore;
    std::string outFileString;
    bool        ok;
};

void ShowOptions();
void ParseCommand(int, char**, ARGS*);

#endif
