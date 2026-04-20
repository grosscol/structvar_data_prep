# Structural Variants Data Processing

Process structural variant data into formats suitable for presentation in BRAVO

- Structural variant calls themselves (bcfs already given)
- Coverage of reads of carriers in each SV region (from crams)
- Data to power SAMPlot style visualization (from crams)

## Data Set Summary Statistics
sample ids:             138,134
max homs per variant:   137,795 (291 hets)
max hets per variant:   137,517 (1 hom)
variants with samples:  482,950

For the variant with the most homozygous calls, there were only 51 samples that were not carriers of the structural variant.

## Stage 1. Generate initial mapping file
This pipeline produces the mapping file that will be subsequently used by downstream pipelines.
Holds INFO fields region and crams for each variant.

1. randomly select 5 of each hets and homs per each variant
2. create bedfile spec for each variant +/- 100 bp
3. merge het/homs and bedfile into single tsv.
4. omit variants with no carriers called
5. make list of all unique subject ids
6. determine which ids do not have locatable crams
7. remove ids that don't have crams
8. substitute in paths to crams for ids.
9. Trim single representative cram per SV to bedfile region specified in the mapping file.

## Stage 2: Cram Processing
1. Continued processing variants with metadata one by one. (See variant processing)

## Stage 3

### Variant Processing
1. Use processed cram file for each variant.
2. Create SAMPLot style summaries of reads.
3. Compile samplot style data into single data set.

### Depth Processing
1. Use processed cram file for each variant.
1. Create pileup and bravo depth data (using subset of bravo pipeline)
1. Concatenate depth data in similar fashion to bravo depth.


## Notes
The position column is 0-based and the variant id column is 1 based
```txt
chr2    0       DEL_2:1-126300          N       <DEL>   -99     126400
chr2    10000   DUP_2:10001-32200       C       <DUP>   9901    32300
```
