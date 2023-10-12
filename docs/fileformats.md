## Summary Statistics

easyfinemap requires only one input file, which is the GWAS summary statistics file. This file typically does not have a fixed format and is essentially a tabular file where each row corresponds to the information and results of the association analysis for a specific SNP. To facilitate processing and speed up analysis, we have designed a unified format based on tabix. The file should have exactly 10 columns, and CHR and BP should be indexed using tabix.

The 10 columns in the file are as follows:

1. CHR: Chromosome (integer)
2. BP: Base pair position (integer)
3. rsID: rsID of the SNP (string, allows null values)
4. EA: Effective allele (string)
5. NEA: Non-effective allele (string)
6. EAF: Effective allele frequency (float, allows null values)
7. MAF: Minor allele frequency (float, allows null values)
8. BETA: Effect size (float)
9. SE: Standard error (float, non-zero values)
10. P: p-value (positive float)

By following this format and indexing CHR and BP using tabix, you can ensure compatibility and efficient processing of the GWAS summary statistics file with easyfinemap.

Users can easily convert summary statistics from other formats into this format using Smunger.

## LD reference (Optional)
To perform LD-based fine-mapping, users need to provide individual genotype data to calculate LD. Ideally, these genotypes should be matched to the sample of summary statistics. However, it is common practice to use publicly available reference panels such as 1000 Genomes (1000G). Since easyfinemap uses PLINK v1.9 to calculate LD, the required genotype data format is PLINK's bfile format.

Users need to split the genotype data by chromosome and convert it to the PLINK bfile format. Then, they can use the `easyfinemap validate-ldref` command to format the LD reference.
