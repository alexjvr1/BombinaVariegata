#Bombina bioinformatics pipeline for 4 test samples

CHR1 & 2
ZIN 1& 2


I received the data back from FGCZ 14 Oct 2015. The samples were run in H20, which had a total of 30 samples

Overall, it looks like the Bombina samples were very under represented. I added equal amounts of DNA into the library. Ideally Bombina should have had 3 times as much DNA in the library as the Rana temporaria samples (since the genome is ~3x bigger), but due to limited DNA, this wasn't possible. 

Subsequently, there are many loci sequenced, but 94% of the sequences are singletons. 

This shouldn't be a problem if samples are equally represented (i.e. all Bombina). 

We should also use a different restriction enzyme so that less loci are sequenced. The size selection could possibly be adjusted. 


Demultiplexing was done using process_radtags from Stacks. 

Adapter dimers removed using Trimmomatic

Both on FGCZ server
```

```

De novo assembly was conducted using pyRAD (on fgcz) at 94% clustering.

This still needs to be optimised for Bombina. 


Subsequent SNP filtering done using VCFtools

```
/usr/local/ngseq/stow/vcftools-v-0.1.14-2/bin/vcftools --vcf BomN4Clust94.vcf --maf 0.25 --recode --recode-INFO-all --out Bom.test

VCFtools - 0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BomN4Clust94.vcf
	--recode-INFO-all
	--maf 0.25
	--out Bom.test
	--recode

Eighth Header entry should be INFO: INFO    
After filtering, kept 4 out of 4 Individuals
Outputting VCF file...
After filtering, kept 535 out of a possible 555 Sites
Run Time = 0.00 seconds

```

Can't filter for MAC of 3, since there are only 4 individuals

```
/usr/local/ngseq/stow/vcftools-v-0.1.14-2/bin/vcftools --vcf BomN4Clust94.vcf --max-missing 0 --mac 1 --recode --recode-INFO-all --out BomN4test 

VCFtools - 0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BomN4Clust94.vcf
	--recode-INFO-all
	--mac 1
	--out BomN4test
	--recode

Eighth Header entry should be INFO: INFO    
After filtering, kept 4 out of 4 Individuals
Outputting VCF file...
After filtering, kept 555 out of a possible 555 Sites
Run Time = 0.00 seconds
```

Filter for a depth of 3, but this has already been filtered in pyRAD

```
vcftools --vcf BomN4test.recode.vcf --minDP 3 --recode --recode-INFO-all --out BomN4testMinDepth3

VCFtools - 0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BomN4test.recode.vcf
	--recode-INFO-all
	--minDP 3
	--out BomN4testMinDepth3
	--recode

After filtering, kept 4 out of 4 Individuals
Outputting VCF file...
After filtering, kept 555 out of a possible 555 Sites
Run Time = 0.00 seconds
```


Assess the missingness
```
vcftools --vcf BomN4test.recode.vcf --missing-indv  ##creates and imiss file which we can then plot

mawk '!/IN/' out.imiss | cut -f5 > totalmissing

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin( $1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
```
Results show that CHR samples are missing 80% of the data (out of 555 sites). 


If I use vcftools to calculate Fst between populations, we get: 
```
--vcf BomN4min10.recode.vcf --weir-fst-pop CHRpop.txt --weir-fst-pop ZINpop.txt --out fstZIN.CHR

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BomN4min10.recode.vcf
	--weir-fst-pop CHRpop.txt
	--weir-fst-pop ZINpop.txt
	--out fstZIN.CHR

After filtering, kept 4 out of 4 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.049769
Weir and Cockerham weighted Fst estimate: 0.28993
After filtering, kept 555 out of a possible 555 Sites
```

Well, that looks pretty good! Higher than I expected. 

Lets see what the Fst is when I mix the individuals up: pop1=CHR01 and ZIN01 vs pop2= CHR02 and ZIN02

```
--vcf BomN4min10.recode.vcf --weir-fst-pop CHR2.ZIN2.txt --weir-fst-pop CHR1.ZIN1.txt --out fstmix1vs2

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BomN4min10.recode.vcf
	--weir-fst-pop CHR2.ZIN2.txt
	--weir-fst-pop CHR1.ZIN1.txt
	--out fstmix1vs2

After filtering, kept 4 out of 4 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: -0.41164
Weir and Cockerham weighted Fst estimate: -0.29295
After filtering, kept 555 out of a possible 555 Sites
```

What do the negative Fst values mean?



If I filter for 0.75 missingness across loci (i.e. 3/4 of the individuals scored per locus), I end up with 69 sites: 
```
--vcf BomN4min10.recode.vcf --max-missing 0.75 --recode --recode-INFO-all --out maxmissloci.75

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BomN4min10.recode.vcf
	--recode-INFO-all
	--max-missing 0.75
	--out maxmissloci.25
	--recode

After filtering, kept 4 out of 4 Individuals
Outputting VCF file...
After filtering, kept 69 out of a possible 555 Sites
Run Time = 0.00 seconds
```


If I then run the same Fst calculations as before: 

```
vcftools --vcf maxmissloci.25.recode.vcf --weir-fst-pop CHRpop.txt --weir-fst-pop ZINpop.txt --out fstZIN.CHR

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf maxmissloci.25.recode.vcf
	--weir-fst-pop CHRpop.txt
	--weir-fst-pop ZINpop.txt
	--out fstZIN.CHR

After filtering, kept 4 out of 4 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.049769
Weir and Cockerham weighted Fst estimate: 0.28993
After filtering, kept 69 out of a possible 69 Sites

```

And for the mixed pops:
```
--vcf maxmissloci.25.recode.vcf --weir-fst-pop CHR2.ZIN2.txt --weir-fst-pop CHR1.ZIN1.txt --out fstmix1vs2

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf maxmissloci.25.recode.vcf
	--weir-fst-pop CHR2.ZIN2.txt
	--weir-fst-pop CHR1.ZIN1.txt
	--out fstmix1vs2

After filtering, kept 4 out of 4 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: -0.41164
Weir and Cockerham weighted Fst estimate: -0.29295
After filtering, kept 69 out of a possible 69 Sites
```
These values are identical. So I wonder if the Fst calculations are based only on loci scored across all individuals??



After reading a bit, I'm worried about the Fst calculation in vcftools. I'm going to transfer the data to my computer, and run some analyses in R.

```

```



```

```
