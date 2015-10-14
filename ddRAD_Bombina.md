#Bombina bioinformatics pipeline for 4 test samples

CHR1 & 2
ZIN 1& 2


I received the data back from FGCZ 14 Oct 2015. The samples were run in H20, which had a total of 30 samples

Overall, it looks like the Bombina samples were very under represented. I added equal amounts of DNA into the library. Ideally Bombina should have had 3 times as much DNA in the library as the Rana temporaria samples (since the genome is ~3x bigger), but due to limited DNA, this wasn't possible. 

Subsequently, 

```

```

```
vcftools --vcf subset.final.vcf --maf 0.25 --recode --recode-INFO-all --out subset.final.maf25
```
