Because fine-mapping is typically conducted locally, determining the regions for fine-mapping a summary statistic requires identifying the regions based on significant loci in a GWAS.

EasyFinemap provides three methods to determine the regions for fine-mapping: distance, clumping, and conditional. The distance method merges nearby variants within a given distance threshold into a locus and selects a lead SNP from each locus. The clumping method merges variants within the same clumping window based on a specified window size and r2 threshold, and selects a lead SNP from each merged locus. The conditional method merges variants within the same cojo window based on a specified window size and cojo collinear threshold, and selects a lead SNP from each merged locus. The key difference among these methods lies in the criteria used for merging variants: distance and clumping methods are based on physical distances between variants, while the conditional method is based on their correlation. Both distance and clumping methods utilize GWAS summary statistics, whereas the conditional method incorporates both GWAS summary statistics and LD reference.

Here is an example of using EasyFinemap for locus identification:

Download the example data (ignore if already downloaded):
```
git clone https://github.com/Jianhua-Wang/easyfinemap.git
cd easyfinemap/exampledata
ls
```
PH251.txt.gz is the GWAS summary statistics file, and PH251.txt.gz.tbi is the tabix index file.

PH251.sig.txt.gz is the significant SNPs file which contains the SNPs with p-value less than 5e-08.

To expedite the process, you can directly use the "sig.txt.gz" file, and the results will be the same as those from "PH251.sig.txt.gz".

### By distance
```
$ easyfinemap get-loci -m  distance PH251.sig.txt.gz PH251.distance
────────────────────────────────── EasyFinemap ───────────────────────────────────
                                  Version: 0.3.9
                               Author: Jianhua Wang
                          Email: jianhua.mert@gmail.com
[22:33:22] INFO     root - Loading PH251.sig.txt.gz...
           INFO     Loci - Save 1 independent loci to PH251.distance.loci.txt
           INFO     Loci - Save 1 independent lead snps to PH251.distance.leadsnp.txt
```

### By clumping

```
$ easyfinemap get-loci -m clumping --ldref ./EUR.valid.chr{chrom} ./PH251.sig.txt.gz PH251.clumping
────────────────────────────────── EasyFinemap ───────────────────────────────────
                                  Version: 0.3.9
                               Author: Jianhua Wang
                          Email: jianhua.mert@gmail.com
[22:35:20] INFO     root - Loading PH251.sig.txt.gz...
           INFO     Loci - Save 1 independent loci to PH251.clumping.loci.txt
           INFO     Loci - Save 1 independent lead snps to
                    PH251.clumping.leadsnp.txt
```
Other parameters can be specified as follows:
```
--clump-kb: the clumping window size, in kb. Default: 500
--clump-r2: the clumping r2 threshold. Default: 0.1
```
### By conditional analysis
```
$ easyfinemap get-loci \
    -m conditional \
    --use-ref-eaf \
    -n 40000 \
    --ldref ./EUR.valid.chr{chrom} \
    PH251.sig.txt.gz \
    PH251.conditional
────────────────────────────────── EasyFinemap ───────────────────────────────────
                                  Version: 0.3.9
                               Author: Jianhua Wang
                          Email: jianhua.mert@gmail.com
[23:07:54] INFO     root - Loading PH251.sig.txt.gz...
Run cojo-slct ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 1/1 0:00:00
           INFO     Loci - Save 1 independent loci to PH251.conditional.loci.txt
           INFO     Loci - Save 1 independent lead snps to
                    PH251.conditional.leadsnp.txt
```
Parameters can be specified as follows:
```
--use-ref-eaf: whether to use the reference panel EAF. Default: False. 
If True, the reference panel EAF will be used for conditional analysis. 
Must be specified when the EAF is not available in the GWAS summary statistics file.

-n: the sample size for conditional analysis.
--cojo-window-kb: the cojo window size, in kb. Default: 10000
--cojo-collinear: the cojo collinear threshold. Default: 0.9
--diff-freq: the difference in allele frequency threshold. Default: 0.2

```
### Other parameters
```
--sig-threshold: the significance threshold. Default: 5e-08
--loci-extension: the extension from lead SNPs, in kb. Default: 500
--ldblock: the path to the LD block file. Default: None, bed file with 3 columns: chr, start, end
```
### Output files
.leadsnp.txt, contains the lead SNPs for each locus. The first line is the header, and the following lines are the lead SNPs for each locus. The columns are:
```
$ head PH251.clumping.leadsnp.txt
SNPID   CHR     BP      rsID    EA      NEA     EAF     MAF     BETA    SE      P
22-29451671-A-G 22      29451671        rs4823006       G       A                       -0.0241 0.0037  7.99e-11
```
.loci.txt, contains the loci for each lead SNP. The first line is the header, and the following lines are the loci for each lead SNP. The columns are:
```
$ head PH251.clumping.loci.txt
CHR     START   END     LEAD_SNP        LEAD_SNP_P      LEAD_SNP_BP
22      28951671        29951671        22-29451671-A-G 7.99e-11        29451671
```
### Help information
```
$ easyfinemap get-loci -h
────────────────────────────────── EasyFinemap ───────────────────────────────────
                                  Version: 0.3.9
                               Author: Jianhua Wang
                          Email: jianhua.mert@gmail.com

 Usage: easyfinemap get-loci [OPTIONS] SUMSTATS_PATH OUTPUT

 Get the loci from the GWAS summary statistics file.

╭─ Arguments ────────────────────────────────────────────────────────────────────╮
│ *    sumstats_path      PATH  The path to the GWAS summary statistics file.    │
│                               [default: None]                                  │
│                               [required]                                       │
│ *    output             TEXT  The output prefix. [default: None] [required]    │
╰────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────────────╮
│ --sig-threshold      -s      FLOAT                    The significance         │
│                                                       threshold.               │
│                                                       [default: 5e-08]         │
│ --loci-extension     -l      INTEGER                  The extension from lead  │
│                                                       SNPs, in kb              │
│                                                       [default: 500]           │
│ --ldblock                    TEXT                     The path to the LD block │
│                                                       file.                    │
│                                                       [default: None]          │
│ --merge-loci         -i                               Whether to merge the     │
│                                                       loci, not recommanded    │
│                                                       for conditional mode.    │
│ --method             -m      [distance|clumping|cond  The method to identify   │
│                              itional]                 the lead SNPs.           │
│                                                       [default:                │
│                                                       LociMethod.distance]     │
│ --distance           -d      INTEGER                  The distance threshold   │
│                                                       for distance method, in  │
│                                                       kb.                      │
│                                                       [default: 50]            │
│ --ldref                      TEXT                     The path to the LD       │
│                                                       reference file.          │
│                                                       [default: None]          │
│ --clump-kb           -c      INTEGER                  The clumping window      │
│                                                       size, in kb.             │
│                                                       [default: 500]           │
│ --clump-r2           -r      FLOAT                    The clumping r2          │
│                                                       threshold.               │
│                                                       [default: 0.1]           │
│ --sample-size        -n      INTEGER                  The sample size for      │
│                                                       conditional method.      │
│                                                       [default: None]          │
│ --cojo-window-kb             INTEGER                  The cojo window size, in │
│                                                       kb.                      │
│                                                       [default: 10000]         │
│ --cojo-collinear             FLOAT                    The cojo collinear       │
│                                                       threshold.               │
│                                                       [default: 0.9]           │
│ --diff-freq                  FLOAT                    The difference in allele │
│                                                       frequency threshold.     │
│                                                       [default: 0.2]           │
│ --use-ref-eaf                                         Whether to use the       │
│                                                       reference panel EAF.     │
│ --only-use-sig-snps                                   Whether to only use the  │
│                                                       significant SNPs.        │
│ --threads            -t      INTEGER                  The number of threads.   │
│                                                       [default: 1]             │
│ --help               -h                               Show this message and    │
│                                                       exit.                    │
╰────────────────────────────────────────────────────────────────────────────────╯
```