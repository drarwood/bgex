# bgex - BEGEN Extraction Tool
This tool has been designed to extract data from bgen files from the UK Biobank and performs the following functions:
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
