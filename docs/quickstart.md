To perform the simplest fine-mapping using easyfinemap, including determining independent loci using the distance-based method and calculating the posterior probability for each SNP using the abf method, you can follow these steps:

1. Download the example data (ignore if already downloaded):
```
git clone https://github.com/Jianhua-Wang/easyfinemap.git
cd easyfinemap/exampledata
```

2. Determine independent loci using the distance-based method:
```
$ easyfinemap get-loci -m  distance PH251.sig.txt.gz PH251.distance
```

3. Calculate the posterior probability for each SNP using the abf method:
```
$ easyfinemap fine-mapping \
    -m abf \
    --credible-threshold 0.95 \
    PH251.txt.gz \
    PH251.distance.loci.txt \
    PH251.distance.leadsnp.txt \
    PH251.abf.txt
```

For more advanced features and parameter explanations, please refer to the [documentation](https://jianhua-wang.github.io/easyfinemap/fileformats/)。