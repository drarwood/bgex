/*
   BGEN Extraction Tool
   Main.cpp
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#include <iostream>
#include "InputParser.h"
#include "FileReader.h"
#include "BGIProcessor.h"
#include "BGENProcessor.h"
#include "Calculator.h"
#include "FileWriter.h"
#include <string>
#include <vector>

using namespace std;

int main(int argc, char** argv) {

    ARGS theseArgs;
    ParseCommand(argc, argv, &theseArgs);
    if(!theseArgs.ok) {
        ShowOptions();
        exit(EXIT_FAILURE);
    }
    // Feedback input
    cout << "Usage: " << endl;
    cout << "  --bgens    -b " << theseArgs.bgenListFileString   << endl;
    cout << "  --samples  -s " << theseArgs.sampleFileString     << endl;
    cout << "  --variants -v " << theseArgs.variantFileString    << endl;
    cout << "  --extract  -e " << theseArgs.subjectIncFileString << endl;
    cout << "  --min-info -i " << theseArgs.minInfo              << endl;
    cout << "  --dosages  -d " << theseArgs.extractDosages       << endl;
    cout << "  --probs    -p " << theseArgs.extractProbs         << endl;
    cout << "  --pscore   -g " << theseArgs.extractPScore        << endl;
    cout << "  --info     -q " << theseArgs.extractInfoScore     << endl;
    cout << "  --out      -o " << theseArgs.outFileString        << endl;
    cout << endl;

    // Read in the list of bgens
    // chromosome id -> bgen file id
    map<string,string> bgen_files;
    ReadBGENListFile(theseArgs.bgenListFileString, bgen_files);

    // Read in the bgen sample file
    // subject index (0-based) -> SAMPLE struct
    map<int,SAMPLE> bgen_samples;
    ReadBGENSampleFile(theseArgs.sampleFileString, bgen_samples);

    // Read in the list of subjects to extract dosages for (if provided)
    // fid:id vector
    vector<string> subject_list;
    if (!theseArgs.subjectIncFileString.empty()) {
        ReadSubjectIncListFile(theseArgs.subjectIncFileString, subject_list);
        FlagSubjectsforExclusion(subject_list, bgen_samples);
    }

    // Read in the variant file to extract/use
    // user defined chr:pos:a1:a2 -> VARIANT struct
    // will use this id when going through bgen to set info
    map<string,VARIANT> variant_list;
    vector<string> chrs;
    ReadVariantListFile(theseArgs.variantFileString, variant_list, chrs);



    // Create vector hold order of variants read across all bgens - will need later when checking to output or not
    vector<string> variants_read_order;

    // Cycle variants per chr to extract and obtain dosages and/or probs
    for (vector<string>::iterator it=chrs.begin(); it!=chrs.end(); ++it) {

        // Get starting byte address for variants for this chromosome
        map<uint64_t,string>& chr_var_byte_starts = GetListofVariantStartBytes(bgen_files[*it], variant_list, *it);
        // need to get byte aaddress to map to user defined chr:pos:ref:alt so can then link to variant_list map
        MapToUserVarIDs(variant_list, chr_var_byte_starts);

        // Read header data of bgen file into struct
        BGEN bgen_info;
        ReadBGENHeader(bgen_files[*it], bgen_info);

        // check format as expected before going on:
        if (bgen_info.sample_count == bgen_samples.size() && bgen_info.layout == 2 && bgen_info.compressed_probs == 1) {
            // Get the genotype probability bytes
            ExtractGenotypeData(bgen_files[*it], chr_var_byte_starts, bgen_samples, variant_list, theseArgs.minInfo, variants_read_order);
        }
        else {
            cout << "Skipping bgen as unexpected format" << endl;
        }

    }

    // cycle through the samples and output probs / dosages / polygenic score flagged as for use
    if (theseArgs.extractProbs) {
        OutputProbs(bgen_samples, variant_list, variants_read_order, theseArgs.outFileString);
    }

    if (theseArgs.extractDosages) {
        OutputDosages(bgen_samples, variant_list, variants_read_order, theseArgs.outFileString);
    }

    if (theseArgs.extractPScore) {
        CalculatePS(bgen_samples, variant_list, variants_read_order);
        OutputPS(bgen_samples, theseArgs.outFileString);
    }

    if (theseArgs.extractInfoScore) {
        OutputInfoScores(variant_list, theseArgs.outFileString);
    }

    return 0;
}
