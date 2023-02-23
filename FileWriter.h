/*
   BGEN Extraction Tool (bgex)
   FileWriter.h
   Functions related to writing human readable files
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#ifndef __FILEWRITER__
#define __FILEWRITER__

#include "FileReader.h"
#include "StringOperations.h"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

void OutputProbs(std::map<int,SAMPLE>&, std::map<std::string,VARIANT>&, std::vector<std::string>&, std::string&);
void OutputDosages(std::map<int,SAMPLE>&, std::map<std::string,VARIANT>&, std::vector<std::string>&, std::string&);
void OutputPS(std::map<int,SAMPLE>&, std::string&);
void OutputInfoScores(std::map<std::string,VARIANT>&, std::string&);

#endif

