easyfinemap supports three types of fine-mapping methods: LD-free, LD-based, and function-based. LD-free methods do not require LD information, LD-based methods rely on LD information, and function-based methods require functional annotation information. Below are the specific software supported for each of these three methods in easyfinemap. The principles of these software are not explained here and can be referred to their respective documentation.

* LD-free fine-mapping
    * [aBF](https://www.nature.com/articles/ejhg20161)
* LD-based fine-mapping
    * [PAINOR](https://github.com/gkichaev/PAINTOR_V3.0)
    * [CAVIARBF](https://bitbucket.org/Wenan/caviarbf)
    * [FINEMAP](http://christianbenner.com/)
    * [SuSiE](https://github.com/stephenslab/susieR)
* Function-based fine-mapping
    * [PolyFun](https://github.com/omerwe/polyfun) + [FINEMAP](http://christianbenner.com/)
    * [PolyFun](https://github.com/omerwe/polyfun) + [SuSiE](https://github.com/stephenslab/susieR)

## LD-free fine-mapping
```
$ easyfinemap fine-mapping \
    -m abf \
    --credible-threshold 0.95 \
    PH251.txt.gz \
    PH251.distance.loci.txt \
    PH251.distance.leadsnp.txt \
    PH251.abf.txt
```
easyfinemap requires four mandatory parameters:

* GWAS summary statistics file
* loci file
* lead SNPs file
* output file

The ABF method was referenced from the article by Asimit, J. L. et al. in Eur J Hum Genet (2016).

$$
ABF = \sqrt{\frac{se}{se + w^2}}\exp\left(\frac{w^2}{se + w^2}\left(\frac{beta^2}{se^2}\right)/2\right)
$$

where w is variance prior (parameter `--var-prior`), usually set to 0.15 for quantitative traits
and 0.2 for binary traits.
the posterior probability of each variant being causal is calculated
using the formula:
PP(causal) = ABF / sum(all_abfs)

Option `--credible-threshold` is used to specify the credible threshold, which is used to determine the credible set. The default value is 0.95. The credible set is determined by the posterior probability of each variant being causal. The variants with the highest posterior probability are included in the credible set until the sum of the posterior probability of the variants in the credible set is greater than or equal to the credible threshold. Will output all variants if the credible threshold is set to None.
## LD-based fine-mapping
```
$ easyfinemap fine-mapping \
    -m finemap \
    --ldref ./EUR.valid.chr{chrom} \
    --use-ref-eaf \
    --credible-threshold 0.95 \
    -n 1000 \
    ./PH251.txt.gz \
    PH251.distance.loci.txt \
    PH251.distance.leadsnp.txt \
    PH251.abf
```
The example mentioned above utilized FINEMAP for fine-mapping, which requires LD information. You can generate LD reference files based on the tutorial provided in the "Validate LD reference" section. The easyfinemap tool automatically generates the necessary input files for each loci and retrieves the results.

Please note that if the summary statistics do not include EAF (Effect Allele Frequency) information, you need to use the `--use-ref-eaf` parameter to use the LD reference's EAF as the EAF in order to avoid any errors.

## Function-based fine-mapping
Users can also perform fine-mapping using the prior probabilities provided by PolyFun. This method requires users to provide functional annotation information, such as gene expression levels and chromatin states. The specific usage details can be found in the PolyFun documentation.

To perform PolyFun's fine-mapping using easyfinemap, it is straightforward. Simply specify the `--prior-file` parameter in the command line, for example:
```
$ easyfinemap fine-mapping \
    -m polyfun_finemap \
    --prior-file ./prior.txt \
    --credible-threshold 0.95 \
    PH251.txt.gz \
    PH251.distance.loci.txt \
    PH251.distance.leadsnp.txt \
    PH251.abf.txt
```
Users can utilize the prior probability file provided by PolyFun or generate their own. For specific usage instructions, please refer to the PolyFun documentation.

The prior probability file provided by PolyFun is in Parquet format. To facilitate GitHub uploading, it is divided into two files: [chr1-7](https://github.com/omerwe/polyfun/blob/master/snpvar_meta.chr1_7.parquet) and [chr8-22](https://github.com/omerwe/polyfun/blob/master/snpvar_meta.chr8_22.parquet). Users need to merge these two files into one and then specify the path to this file using the `--prior-file` parameter.

The format is as follows:
```
CHR     BP      SNP     A1      A2      snpvar_bin
1       10177   rs367896724     A       AC      0.0
1       10352   rs201106462     T       TA      0.0
1       10511   rs534229142     G       A       0.0
1       10616   1:10616_CCGCCGTTGCAAAGGCGCGCCG_C        CCGCCGTTGCAAAGGCGCGCCG  C       0.0
1       11008   rs575272151     C       G       0.0
1       11012   rs544419019     C       G       0.0
1       13110   rs540538026     G       A       0.0
1       13116   rs62635286      T       G       0.0
1       13118   rs62028691      A       G       0.0
```


## Fine-mapping with conditional analysis
If users have used conditional analysis to determine the lead SNP of a locus, they can utilize the `--conditional` parameter to perform fine-mapping using the results of conditional analysis. For example:
```
$ easyfinemap fine-mapping \
    -m finemap \
    --credible-threshold 0.95 \
    --conditional \
    --ldref ./EUR.valid.chr{chrom} \
    --use-ref-eaf \
    -n 1000 \
    PH251.txt.gz \
    PH251.conditional.loci.txt \
    PH251.conditional.leadsnp.txt \
    PH251.abf.txt
```
easyfinemap will employ GCTA COJO to perform conditional analysis. For each locus, conditional analysis generates new summary statistics by conditioning on other lead SNPs within a specified interval. This interval is determined by the '--cond-snps-wind-kb' parameter, which has a default value of 10000 (10 Mb).

## Fine-mapping with multiple methods
In summary, easyfinemap supports multiple fine-mapping methods. Users can specify a single method or multiple methods to perform fine-mapping simultaneously by using multiple `-m` parameters in the command line. For example:
```
$ easyfinemap fine-mapping \
    -m abf -m finemap \
    --ldref ./EUR.valid.chr{chrom} \
    --use-ref-eaf \
    --credible-threshold 0.95 \
    --credible-method finemap \
    -n 1000 \
    ./PH251.txt.gz \
    PH251.distance.loci.txt \
    PH251.distance.leadsnp.txt \
    PH251.abf
```
You can also specify `-m` all to perform fine-mapping using all supported methods. For example:
```
$ easyfinemap fine-mapping \
    -m all \
    --ldref ./EUR.valid.chr{chrom} \
    --use-ref-eaf \
    --credible-threshold 0.95 \
    --credible-method finemap \
    -n 1000 \
    ./PH251.txt.gz \
    PH251.distance.loci.txt \
    PH251.distance.leadsnp.txt \
    PH251.abf
```
Please note that when using multiple methods for fine-mapping, you need to specify the `--credible-method` parameter to determine which method should be used to determine the credible set.

## Output
The output file is a tab-delimited text file. The first line is the header, which contains the following columns:
```
$ head PH251.finemap.txt
SNPID   CHR     BP      rsID    EA      NEA     EAF     MAF     BETA    SE      P       PP_FINEMAP      LEAD_SNP
22-29451671-A-G 22      29451671        rs4823006       G       A                       -0.0241 0.0037  7.99e-11        0.917195        22-29451671-A-G
22-29449477-A-G 22      29449477        rs2294239       G       A                       -0.0243 0.004   8.26e-10        0.0721721       22-29451671-A-G
```
The first 11 columns are the same as the GWAS summary statistics file. The `PP_FINEMAP` column is the posterior probability of each variant being causal. The `LEAD_SNP` column is the lead SNP of the locus.

## Other parameters
`--max-causal` parameter is used to specify the maximum number of causal SNPs. The default value is 1.

`--cond-snps-wind-kb` parameter is used to specify the conditional SNPs window size, in kb. The default value is 10000.

`--threads` parameter is used to specify the number of threads. The default value is 1.

## Help information

```
$ easyfinemap fine-mapping -h
────────────────────────────────── EasyFinemap ───────────────────────────────────
                                  Version: 0.4.0
                               Author: Jianhua Wang
                          Email: jianhua.mert@gmail.com

 Usage: easyfinemap fine-mapping [OPTIONS] SUMSTATS_PATH LOCI_PATH
                                 LEAD_SNPS_PATH OUTFILE

 Fine mapping.

╭─ Arguments ────────────────────────────────────────────────────────────────────╮
│ *    sumstats_path       TEXT  The path to the GWAS summary statistics file.   │
│                                [default: None]                                 │
│                                [required]                                      │
│ *    loci_path           TEXT  The path to the loci file, generated by         │
│                                get-loci command.                               │
│                                [default: None]                                 │
│                                [required]                                      │
│ *    lead_snps_path      TEXT  The path to the lead SNPs file, generated by    │
│                                get-loci command.                               │
│                                [default: None]                                 │
│                                [required]                                      │
│ *    outfile             TEXT  The output file. [default: None] [required]     │
╰────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────────────╮
│ *  --methods             -m      [abf|finemap|paintor|  The methods to use.    │
│                                  caviarbf|susie|polyfu  [default: None]        │
│                                  n_susie|polyfun_finem  [required]             │
│                                  ap|all]                                       │
│    --var-prior                   FLOAT                  The prior variance for │
│                                                         the aBF method.        │
│                                                         [default: 0.2]         │
│    --conditional         -c                             Whether to use         │
│                                                         conditional mode.      │
│    --prior-file                  TEXT                   The path to the prior  │
│                                                         file.                  │
│                                                         [default: None]        │
│    --sample-size         -n      INTEGER                The sample size for    │
│                                                         conditional mode.      │
│                                                         [default: None]        │
│    --ldref                       TEXT                   The path to the LD     │
│                                                         reference file.        │
│                                                         [default: None]        │
│    --cond-snps-wind-kb           INTEGER                The conditional SNPs   │
│                                                         window size, in kb.    │
│                                                         [default: 10000]       │
│    --max-causal                  INTEGER                The maximum number of  │
│                                                         causal SNPs.           │
│                                                         [default: 1]           │
│    --credible-threshold          FLOAT                  The credible           │
│                                                         threshold.             │
│                                                         [default: None]        │
│    --credible-method             TEXT                   The fine-mapping       │
│                                                         method for credible    │
│                                                         set.                   │
│                                                         [default: None]        │
│    --use-ref-eaf                                        Whether to use the     │
│                                                         reference panel EAF.   │
│    --threads             -t      INTEGER                The number of threads. │
│                                                         [default: 1]           │
│    --help                -h                             Show this message and  │
│                                                         exit.                  │
╰────────────────────────────────────────────────────────────────────────────────╯
```