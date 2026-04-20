# Initial Data Cleaning and Reshaping
This pipeline produces two products: 
- bcf that has been filtered and retagged
- tsv of the INFO fields from the bcf

## Retag and trim
Remove variants with no called genotype.  Recalculate the counts and frequencies from present genotypes.
It is likely the list of samples that were included in the study changed.
This meant some variants that previously had genotypes called, no longer had samples supporting that call.

