/*
   BGEN Extraction Tool (bgex)
   BGIProcessor.h
   Functions related to finding byte addresses of variants in BGEN files
   Also maintains user defined chr:pos:a1:a2
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#ifndef __BGIPROCESSOR__
#define __BGIPROCESSOR__

#include <string>
#include <vector>
#include <map>
#include "FileReader.h"
#include <stdint.h>
#include <algorithm>

std::map<uint64_t, std::string>& GetListofVariantStartBytes(std::string&, std::map<std::string,VARIANT>&, std::string&);
void MapToUserVarIDs(std::map<std::string, VARIANT>&, std::map<uint64_t, std::string>&);


#endif
