/*
   BGEN Extraction Tool (bgex)
   Calculator.h
   Functions related to calculations involving genotypes
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#ifndef __CALCULATOR__
#define __CALCULATOR__

#include <map>
#include <string>
#include <vector>

void CalculatePS(std::map<int,SAMPLE>&, std::map<std::string,VARIANT>&, std::vector<std::string>&);

#endif
