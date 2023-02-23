/*
   BGEN Extraction Tool (bgex)
   StringOperations.h
   Functions related to converting to or from strings
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#ifndef ___STRINGOPERATIONS___
#define ___STRINGOPERATIONS___

#include <string>
#include <vector>

double      stringToDouble(std::string);
int         stringToInt(std::string);
std::string intToString(int);
std::string doubleToString(double);
void        split(std::string&, std::string, std::vector<std::string>&);

#endif


