/*
   BGEN Extraction Tool (bgex)
   FileReader.h
   Functions related to reading human readable files
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#ifndef __FILEREADER__
#define __FILEREADER__

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "StringOperations.h"

struct VARIANT
{
    std::string chr;
    std::string pos;
    std::string a1;
    std::string a2;
    double a1_weight;
    std::string bgen_a1;
    std::string bgen_a2;
    double info_score;
    bool use = false; // later true if found and >= min info threshold
};

struct SAMPLE
{
    std::string id;
    std::vector<unsigned char> prob_bytes;
    std::vector<double> probs;
    double ps;
    bool use = true; // later false if not in an inclusion list
};

void ReadBGENListFile(std::string&, std::map<std::string,std::string>&);
void ReadBGENSampleFile(std::string&, std::map<int,SAMPLE>&);
void ReadVariantListFile(std::string&, std::map<std::string,VARIANT>&, std::vector<std::string>&);
void ReadSubjectIncListFile(std::string&, std::vector<std::string>&);
void FlagSubjectsforExclusion(std::vector<std::string>&, std::map<int,SAMPLE>&);

#endif

