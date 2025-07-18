/*
   BGEN Extraction Tool (bgex)
   BGENProcessor.cpp
   Functions related to reading and storing BGEN data (UK Biobank)
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <math.h>
#include "zlib.h"
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include "zstd/lib/zstd.h"
#include "zstd/examples/common.h"

#include "FileReader.h"
#include "BGENProcessor.h"

using namespace std;

uint32_t get_uint32(char* b) {
    uint32_t x = int((unsigned char)(b[3]) << 24 |
                     (unsigned char)(b[2]) << 16 |
                     (unsigned char)(b[1]) << 8  |
                     (unsigned char)(b[0]));
    return x;
}

uint16_t get_uint16(char* b) {
    uint16_t x = int((unsigned char)(b[1]) << 8 |
                     (unsigned char)(b[0]));
    return x;
}


string chars2string(char* a, int size) {
    string s = "";
    for (int i = 0; i < size; ++i) {
        s += a[i];
    }
    return s;
}



void ReadBGENHeader(string& f, struct BGEN& b) {


    // Get byte address of first variant relative to the 5th byte in the file
    char buff4[4];

    // Open the bgen file:
    FILE* bgen = fopen(f.c_str(), "rb");
    if (bgen != NULL) {

        (void)!fread(&buff4, 1, sizeof(buff4), bgen);
        b.first_var_byte_add = get_uint32(buff4)+4;
        //cout << "First byte address: " << b.first_var_byte_add << endl;

        // Move to the header block
        fseek(bgen, 4, SEEK_SET);

        // Get total size of header block
        (void)!fread(&buff4, 1, sizeof(buff4), bgen);
        b.header_length = get_uint32(buff4);
        //cout << "Size of header block (bytes): " << b.header_length << endl;

        // Get number of variant blocks in bgen
        (void)!fread(&buff4, 1, sizeof(buff4), bgen);
        b.var_count = get_uint32(buff4);
        //cout << "Number of variant blocks stored in file: " << b.var_count << endl;

        // Get sample count
        (void)!fread(&buff4, 1, sizeof(buff4), bgen);
        b.sample_count = get_uint32(buff4);
        //cout << "Sample count in header: " << b.sample_count << endl;

        // Get magic number ('b','g','e','n')
        (void)!fread(&buff4, 1, sizeof(buff4), bgen);
        b.magic_number = chars2string(buff4, sizeof(buff4));
        //cout << "Magic number: " << b.magic_number << endl;

        // skip over free data block area if required: go to 32+freesizeblock
        b.free_area_size = b.header_length - 20;
        if (b.free_area_size > 0) {
            fseek(bgen, 32+b.free_area_size, SEEK_SET);
        }

	// read final 4 bytes of the header which will contain various flags
        (void)!fread(&buff4, 1, sizeof(buff4), bgen);
        // bits 0-1 flag for decompressed (0) or compressed (1) genotype probabilities
        b.compressed_probs = buff4[0] & 0x0003;
        //cout << "Compressed probability flag: " << b.compressed_probs << endl;

        // bits 2-5 provide layout code:
        // (0) allows only single-character alleles,
        // (1) allows multiple characters in alleles,
        // (2) allows multiple alleles, phased and unphased genotypes, explicit specification of ploidy and missing data, and configurable l$
        b.layout =  (buff4[0] >> 2) & 0x000F;
        //cout << "Layout flag: " << b.layout << endl;

        // bit 31 flags whether sample identifiers are present in the bgen file (0) = no, (1) = yes and sample block follows
        b.samples_in_bgen = (buff4[3] >> 7) & 0x0001;
        //cout << "Samples in bgen flag: " << b.samples_in_bgen << endl;

    }

    fclose(bgen);

}


void ExtractGenotypeData(string& f,
                         map<uint64_t,string>& b,
                         map<int, SAMPLE>& s,
                         map<string,VARIANT>& vl,
                         double& mi,
                         vector<string>& vo,
                         unsigned int cFlag) {

    /* This method reads the bgen files and extracts the dosages / probabilities.
       Dosages and probabilities are stored within the SAMPLE structs
       Parameters:
         f:     ref to string holding name of bgen file
         b:     ref to vector holding byte start addresses of variants
         s:     ref to map containing ordered sample index to SAMPLE struct
         vl:    ref to map of user defined chr:pos:a1:A2 -> VARIANT STRUCT for this chromosome
         mi:    ref to double for min info threshold for variant to be output / included
         vo:    ref to vector to push_back on to hold user variant ids for genotypes extracted ordered as read
         cFlag: integer: 1=zlib compressed, 2=zstd compressed.
    */

    // Open the bgen file:
    FILE* bgen = fopen(f.c_str(), "rb");
    if (bgen != NULL) {

        // Get byte address of first variant relative to the 5th byte in the file
        char buff4[4];
        char buff2[2];

        // AT THIS POINT WE MUST GO TO THE BYTE ADDRESS OF THE NEXT VARIANT WE WANT TO GET INFORMATION FOR
        // Cycle through:
        for (map<uint64_t,string>::iterator it=b.begin(); it!=b.end(); ++it) {

            // push back next user defined variant id onto vector for later use
            vo.push_back(it->second);
            //cout << "Next variant in bgen processed" << endl;
            //cout << it->second << endl;

            // go to starting byte
            fseek(bgen, it->first, SEEK_SET);

            // set up buffer for variant identifier and fetch
            (void)!fread(&buff2, 1, sizeof(buff2), bgen);
            uint16_t varid_length = get_uint16(buff2);
            char varid[varid_length];;
            (void)!fread(&varid, 1, sizeof(varid), bgen);
            string varid_s = chars2string(varid, varid_length);
            //cout << "variant id: " << varid_s << endl;

            // set up buffer for rsid
            (void)!fread(&buff2, 1, sizeof(buff2), bgen);
            uint16_t rsid_length = get_uint16(buff2);
            char rsid[rsid_length];
            (void)!fread(&rsid, 1, sizeof(rsid), bgen);
            string rsid_s = chars2string(rsid, sizeof(rsid));
            //cout << "rsid:       " << rsid_s << endl;
            // set up buffer for chromosome
            (void)!fread(&buff2, 1, sizeof(buff2), bgen);
            uint16_t chr_length = get_uint16(buff2);
            char chr[chr_length];
            (void)!fread(&chr, 1, sizeof(chr), bgen);
            string chr_s = chars2string(chr, sizeof(chr));
            if (chr_s.substr(0,1) == "0") {
                chr_s = chr_s.substr(1,1);
            }
            //cout << "chr:        " << chr_s << endl;
            // position of variant
            (void)!fread(&buff4, 1, sizeof(buff4), bgen);
            uint32_t bp = get_uint32(buff4);
            //cout << "bp:         " << bp << endl;
            // number of alleles
            (void)!fread(&buff2, 1, sizeof(buff2), bgen);
            //uint16_t allele_n = get_uint16(buff2);
            //cout << "alleles:    " << allele_n << endl;
            // get allele 1
            (void)!fread(&buff4, 1, sizeof(buff4), bgen);
            uint32_t a1_length = get_uint32(buff4);
            char a1[a1_length];
            (void)!fread(&a1, 1, sizeof(a1), bgen);
            string a1_s = chars2string(a1, sizeof(a1));
            //cout << "A1:         " << a1_s << endl;
            // get allele 2
            (void)!fread(&buff4, 1, sizeof(buff4), bgen);
            uint32_t a2_length = get_uint32(buff4);
            char a2[a2_length];
            (void)!fread(&a2, 1, sizeof(a2), bgen);
            string a2_s = chars2string(a2, sizeof(a2));
            //cout << "A2:         " << a2_s << endl;

            // define bgen varid based on chr:pos:a1a2 here
            string bgen_id = chr_s + ":" + to_string(bp) + ":" + a1_s + ":" + a2_s;
            //cout << "bgen id: " << bgen_id << endl;
            // Now we are at the compressed block
            // 1. Total length of the remainder of the block
            (void)!fread(&buff4, 1, sizeof(buff4), bgen);
            unsigned long int remaining_block_length = int((unsigned char)(buff4[3]) << 24 | (unsigned char)(buff4[2]) << 16 | (unsigned char)(buff4[1]) << 8 | (unsigned char)(buff4[0]));
            //cout << "remaining_block_length:    " << remaining_block_length << endl;
            // 2. Total length of the prob data block uncompressed
            (void)!fread(&buff4, 1, sizeof(buff4), bgen);
            unsigned long int uncompressed_block_length = int((unsigned char)(buff4[3]) << 24 | (unsigned char)(buff4[2]) << 16 | (unsigned char)(buff4[1]) << 8 | (unsigned char)(buff4[0]));
            //cout << "uncompressed_block_length: " << uncompressed_block_length << endl;
            // 3. Size of compressed block
            unsigned long int compressed_block_length = remaining_block_length-4;
            //cout << "compressed_block_length:   " << compressed_block_length << endl;

            // get current read poisiton
            //auto read_pos = ftell(bgen);
            //cout << "Current read position: " << read_pos << endl;

            // Get compressed block
            unsigned char comp_block[compressed_block_length];
            (void)!fread(&comp_block, 1, sizeof(comp_block), bgen);

            unsigned char uncompressed_block[uncompressed_block_length];
            // decompress prob block based on compression algorithm
            if (cFlag == 1) {
                // ZLIB compression:
                // generate buffer to hold decompressed block
                //unsigned char uncompressed_block[uncompressed_block_length];
                // decompress using zlib
                uncompress(uncompressed_block, &uncompressed_block_length,  comp_block, sizeof(comp_block));
            }
            else if (cFlag == 2) {
                // ZSTD compression:
                ZSTD_DStream* const dstream = ZSTD_createDStream();
                ZSTD_initDStream(dstream);
                size_t const buffOutSize = uncompressed_block_length;
                void*  const buffOut     = malloc_orDie(buffOutSize);
                ZSTD_inBuffer input = { comp_block,  compressed_block_length, 0 };
                ZSTD_outBuffer output = { buffOut, uncompressed_block_length, 0 };
                //std::size_t toRead = ZSTD_decompressStream(dstream, &output, &input);
                (void)!ZSTD_decompressStream(dstream, &output, &input);
                //unsigned char* uncompressed_block = (unsigned char*)buffOut;
                unsigned char* t = (unsigned char*)buffOut;
                memcpy(uncompressed_block, t, uncompressed_block_length);
            }


            // get first 4 bytes from buffer to confirm N
            for (int i = 0; i < 4; ++i) {
                buff4[i] = uncompressed_block[i];
            }
            uint32_t varN = get_uint32(buff4);
            //cout << "Number of individuals " << varN << endl;

            // get number of alleles for this variant
            //for (int i = 4; i < 6; ++i) {
            //    buff2[i-4] = uncompressed_block[i];
            //}
            //uint16_t alleleN = get_uint16(buff2);
            //cout << "Number of alleles " << alleleN << endl;

            // get min ploidy
            //unsigned int min_ploidy = uncompressed_block[6];
            //cout << "Min ploidy " << min_ploidy << endl;
            // get max ploidy
            //unsigned int max_ploidy = uncompressed_block[7];
            //cout << "Max ploidy " << max_ploidy << endl;

            // phased flag
            ///unsigned int phased_flag = uncompressed_block[varN+8];
            //cout << "Phased flag: " << phased_flag << endl;

            // Get the number of bits representing each prob
            unsigned int bits_per_prob = uncompressed_block[8+varN+1];
            //cout << "Number of bits per probability " << bits_per_prob << endl;


            // cycle through the rest of the buffer getting the probabilites
            // probabilities are stored as: P1aa, P1ab, P2aa, P2ab,...PNaa,PNab
            int index             = 0; // index of subject within BGEN
            double info_numerator = 0; // numerator for INFO calculation
            double dosage_sum     = 0; // sum of dosages among individuals for which data extract
            double sample_n       = 0; // number of individuals for which genotypes extract

            // cycle through individuals
            for (size_t i = 10+varN; i < uncompressed_block_length; i = i+2) {

                // check if need to process this individual based on possible inclusion list
                if (s[index].use) {

                    // get genotype probabilities
                    double paa = uncompressed_block[i]   / (pow(2.0, bits_per_prob) - 1);
                    double pab = uncompressed_block[i+1] / (pow(2.0, bits_per_prob) - 1);

                    // update numerator for INFO-score calculation
                    info_numerator += pab + (4*paa) - (pow(pab+(2*paa),2.0) );
                    // add to sum of genotype dosages
                    dosage_sum += pab + (2*(1-(paa+pab)));
                    // increment sample N - needed for final INFO score
                    sample_n++;

                    // push back to sample
                    s[index].prob_bytes.push_back(uncompressed_block[i]);
                    s[index].prob_bytes.push_back(uncompressed_block[i+1]);

                }
                // update subject index
                index++;
            } // done for all individuals for this variant

            // Calculate INFO score
            double freq=dosage_sum/(2*sample_n);
            double info = 1.0;
            if (freq > 0 && freq < 1) {
                info = 1 - (info_numerator / ( 2*sample_n*freq*(1-freq)) );
            }

            vl[it->second].bgen_a1 = a1_s;
            vl[it->second].bgen_a2 = a2_s;
            vl[it->second].info_score = info;
            if (info >= mi) {
                vl[it->second].use = true;
            }
        }
    }
}
