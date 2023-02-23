/*
   BGEN Extraction Tool (bgex)
   BGENProcessor.h
   Functions related to reading and storing BGEN data (UK Biobank)
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#ifndef __BGENPROCESSOR__
#define __BGENPROCESSOR__

#include <string>

struct BGEN
{
    uint32_t    first_var_byte_add;
    uint32_t    header_length;
    uint32_t    var_count;
    uint32_t    sample_count;
    std::string magic_number;
    uint32_t    free_area_size;
    int         compressed_probs;
    int         layout;
    int         samples_in_bgen;
};



void ReadBGENHeader(std::string&, struct BGEN&);
void ExtractGenotypeData(std::string&,
                         std::map<uint64_t, std::string>&,
                         std::map<int,SAMPLE>&, std::map<std::string,VARIANT>&,
                         double&, std::vector<std::string>&);

#endif
