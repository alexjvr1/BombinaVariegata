#Bombina Final dataset analysis


pyRAD was run on BV2 and 3 combined. 

Some of these samples are from 2012 and might have poor quality DNA. 

1. Determine the quality of the dataset

2. What is the minimum number of reads if I want to subset the data. 

3. Which samples need to be included in run BV4?


###Filter BV23


50% genotyping rate. And MAC of 3. (across 48 indivs = 3/96 = 3.1%)

```
vcftools --vcf BV2and3.vcf --max-missing 0.5 --mac 3 --recode --recode-INFO-all --out BV2and3.s2


Parameters as interpreted:
	--vcf BV2and3.vcf
	--recode-INFO-all
	--mac 3
	--max-missing 0.5
	--out BV2and3.s2
	--recode

Eighth Header entry should be INFO: INFO    
After filtering, kept 48 out of 48 Individuals
Outputting VCF file...
After filtering, kept 29368 out of a possible 717200 Sites
Run Time = 8.00 seconds
```

Calculate missingness: 

```
vcftools --vcf BV2and3.s2.recode.vcf --missing-indv

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
