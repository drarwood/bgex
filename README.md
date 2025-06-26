# bgex - BGEN Extraction Tool
### Written by Andrew R Wood (University of Exeter)
This C++ tool has been designed to extract data from binary bgen files from the UK Biobank and performs the following functions:
 * extract genotype dosages
 * extract genotype probabilties
 * calculate polygenic risk scores
 * calculate imputation quality scores (INFO-scores)

The above functions work on a list of individuals and a list of variants.
Note this version expects layout 2 with compression blocks derived through either zlib or zstd libraries.
The present version has been tested on UK Biobank bgen files storing imputed data derived from HRC+UK10 and TOPMed imputation panels on a cloud workstation on the DNAnexus platform. 
## Options
```
    --bgens     -b  [bgen and sample file list]
    --variants  -v  [file of variants to extract]
    --samples   -s  [.sample file]
    --extract   -e  [file of samples to keep (optional)]
    --min-info  -m  [min INFO for genotype extraction/use]
    --dosages   -d  [flag to extract dosages]
    --probs     -p  [flag to extract probabilities]
    --pscore    -g  [flag to extract polygenic score]
    --info      -q  [flag to extract info score]
    --out       -o  [file prefix for outputs]
```
At least one of `--dosages` `--probs` `--pscore` `--info`  must be provided.

## Compiling and install
```
# create binary
make
# install system-wide (if permissions)
make install
```

## Consider using DXFUSE to avoid downloading the .bgen and .bgi files if using UKB data
A major bottleneck in working with UKB genetic data is downloading it to a cloud workstation first. 
To avoid this, you should consider using streaming the data through `DXFUSE` instead of downloading prior to using `bgex` (although you do not have to do this to use `bgex`).
It's relatively quick and straightforward to set this up - here are some commands to faciliate this.

```
# Initialise  workspace (if not done already:
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

# Get latest dxfuse binary and install
wget https://github.com/dnanexus/dxfuse/releases/download/v1.5.0/dxfuse-linux
sudo mv dxfuse-linux /usr/bin/dxfuse
sudo chmod 777 /usr/bin/dxfuse

# create directory to mount UKB project to (e.g. "project/")
mkdir project/

# mount project to directory
dxfuse project "UKB_500k_WGS"

```
The example input files listing bgens within the `example_input/` directory reflect this mounting for the example UKB-RAP project called `UKB_500k_WGS`. 


## Input Files
### BGEN file list
This should be a tab-delimited file containing chromosome and absolute file path to respective bgen: 
```
1	/home/dnanexus/project/UKB_500k_WGS/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3.bgen
2	/home/dnanexus/project/UKB_500k_WGS/Bulk/Imputation/UKB imputation from genotype/ukb22828_c2_b0_v3.bgen
3	/home/dnanexus/project/UKB_500k_WGS/Bulk/Imputation/UKB imputation from genotype/ukb22828_c3_b0_v3.bgen
...
```
If the file path contains spaces, do not try to escape them or enclose the filename in quotes (see example files provided)
Note the respective `.bgi` is required and expected to be in the same directory as the `.bgen` file.

### Variant file
This tab-delimited file should contain at least the chromosome, bp-position, allele1, allele2. A 5th column may be provided that specifies the beta or log(OR) aligned to allele1. The 5th column will only be used if `--pscore` specified:
```
1	4414033	C	T	0.006406851
1	7821917	C	CT	0.006961472
1	16045250	A	G	0.00648666
1	17958038	C	T	0.016447903
1	32180167	T	A	0.00900976
...
```
Note that polygenic scores will be derived based on weighting phenotype raising alleles to the absolute values of the 5th column.


### Optional subjects inclusion file
This file should contain a family id and individual id that will be searched for in the bgen sample file
```
1234567	1234567
2323433	2323433
...
```
Note, info-scores are calculated based on individuals included in this file.


## Output Files
### .dosages (requires `--dosages` flag)
The format of the `.dosages` file is:
```
fid:iid           var1                var2               var3               ...  varN
1234567:1234567   p(a1a2)+2p(aa2a2)   p(a1a2)+2p(aa2a2)  p(a1a2)+2p(aa2a2)  ...  p(a1a2)+2p(aa2a2)
...
```
The variant IDs are of the form `chr:pos:allele_1:allele_2` where `allele_1` and `allele_2` are defined by the bgen format - not the user. The dosage increasing allele is `allele_2`.

### .probs (requires `--probs` flag)
The `.probs` file contains genotype probability pairs for `allele_1/allele_1` and `allele_1/allele_2`. The probability of being homozygous for `allele_2` can be derived by substracting the sum of the two probabilities from 1. The format of the `.probs' file is:
```
fid:iid           var1              var2	      var3              ...   varN
1234567:1234567   p(a1a1),p(a1a2)   p(a1a1),p(a1a2)   p(a1a1),p(a1a2)   ...   p(a1a1),p(a1a2)
...
```
The variant IDs are of the form `chr:pos:allele_1:allele_2` where `allele_1` and `allele_2` are defined by the bgen format - not the user.

### .pscores (requires `--pscores` flag)
The `.pscores` file contains the derived polygenic scores based on alleles and weights provided in the variant file (see above). Polygenic scores are calculated by summing the number of trait raising alleles multiplied by the respective abs(weight). 
The format of the `.pscores` file is:
```
fid:iid          pscore
1234567:1234567  0.8939
...
```

### .infoscores (requires `--info` flag)
The `.infoscores` file contains the calculated INFO-score metric commonly used as a measure of genotype imputation quality. The INFO-scores for a given variant are calcuated based on either all individuals represented in the bgen or based on those listed in the subject inclusion file (see above).
The format of the `.infoscore` file is:
```
variant       info_score
1:10235:T:TA  0.8990
...
```
Note the variant ID is based on the original user chr, bp position, a1, and a2 in the variant file input by the user.


## Example command:
Assuming the `DXFUSE` configuration above, the following command will produce all output files for all individuals in BGEN. Info scores will be based on all individuals.
```
./bgex \
  --bgens    example_input/hrc_uk10k_bgens.txt \
  --samples  "/home/dnanexus/project/UKB_500k_WGS/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3.sample"  \
  --variants example_input/var_list_b37.txt \
  --out      my_output_file_prefix \
  --probs \
  --dosages \
  --pscore \
  --info
```

The following command will produce all output files for all individuals in the optional subject list file (this is not provided in the `example_input/` directory). 
Info scores will be based on individuals extracted. Only variants with an INFO score >=0.4 will be output or used for polygenic score deriviation.
```
./bgex \
  --bgens    example_input/hrc_uk10k_bgens.txt \
  --samples  "/home/dnanexus/project/UKB_500k_WGS/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3.sample"  \
  --variants example_input/var_list_b37.txt \
  --out      my_output_file_prefix \
  --extract  my_subject_list.txt \
  --probs \
  --dosages \
  --pscore \
  --info \
  --min-info 0.4 
```

