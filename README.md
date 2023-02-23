# bgex - BEGEN Extraction Tool
This tool has been designed to extract data from binary bgen files from the UK Biobank and performs the following functions:
 * extract genotype dosages
 * extract genotype probabilties
 * calculate polygenic risk scores
 * calculate imputation quality scores (INFO-scores)

The above functions work on a list of individuals and a list of variants.

## Options
```
    --bgens     -b  [bgen list file]
    --samples   -s  [bgen sample file]
    --variants  -v  [file of variants to extract]
    --extract   -e  [file of samples to keep (optional)]
    --min-info  -m  [min INFO for genotype extraction/use]
    --dosages   -d  [flag to extract dosages]
    --probs     -p  [flag to extract probabilities]
    --pscore    -g  [flag to extract polygenic score]
    --info      -q  [flag to extract info score]
    --out       -o  [file prefix for outputs]
```

## Compiling manually (make file due)
```
# compile sqlite
cd sqlite/
gcc -Wall -pthread sqlite3.c -c
# compile main program
g++ -O2 -Wall -pthread *.cpp -c  
g++ -O2 -pthread -o bgex *.o sqlite/sqlite3.o -ldl -I./zlib-1.2.13 -lz
```

## Input Files
### BGEN file list
This should be a file containing chromosome and absolute file path to respective bgen: 
```
1	/full/path/to/ukb_imp_chr1_v3.bgen
2	/full/path/to/ukb_imp_chr2_v3.bgen
3	/full/path/to/ukb_imp_chr3_v3.bgen
...
```

### Variant file
This file should contain at least the chromosome, bp-position, allele1, allele2. A 5th column may be provided that specifies the beta or log(OR) aligned to allele1. The 5th column will only be used if `--pscore` specified:
```
1	10235	T	TA	0.2
16	53818708	T	G	-0.22
16	53798523	A	G	0.83
16	53831146	T	C	-0.1
18	57838401	A	G	-0.123
2	25150296	G	A	0.828
2	632348	A	G	0.009
1	177889480	A	G	-0.013
11	27679916	T	C	0.88
```
Note that polygenic scores will be derived based on weighting phenotype raising alleles to the absolute values of the 5th column.


### Optional subjects inclusion file
This file should contain a family id and individual id that will be searched for in the bgen sample file
```
1234567	1234567
2323433	2323433
...
```
Note, INFO-scores are calculated based on individuals included in this file.

