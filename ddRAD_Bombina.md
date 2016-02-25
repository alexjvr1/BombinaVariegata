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

Convert the .vcf into a Plink file

```
vcftools --vcf maxmissloci.75.recode.vcf --out Bomb69.PLINK --plink

plink --file Bomb69.PLINK --out Bomb69test --recodeA


Writing this text to log file [ Bomb69test.log ]
Analysis started: Thu Oct 15 13:36:29 2015

Options in effect:
	--file Bomb69.PLINK
	--out Bomb69test
	--recodeA

69 (of 69) markers to be included from [ Bomb69.PLINK.map ]
Warning, found 4 individuals with ambiguous sex codes
Writing list of these individuals to [ Bomb69test.nosex ]
4 individuals read from [ Bomb69.PLINK.ped ] 
0 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
0 cases, 0 controls and 4 missing
0 males, 0 females, and 4 of unspecified sex
Before frequency and genotyping pruning, there are 69 SNPs
4 founders and 0 non-founders found
Total genotyping rate in remaining individuals is 0.789855
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 69 SNPs
After filtering, 0 cases, 0 controls and 4 missing
After filtering, 0 males, 0 females, and 4 of unspecified sex
Writing recoded file to [ Bomb69test.raw ] 
```


copy everything over to the mac
```
 scp -r alexjvr@gdcsrv1.ethz.ch:/gdc_home4/alexjvr/Bombina/Filtered/* .
```


##2 Nov 2015

Aim: Prepare 1 library for HiSeq with the following conditions: 

1. 24 individually barcoded samples 
2. PstI + EcoRI enzymes
3. Selective Illumina adapters? - The selective Illumina primer that I'm using selectively sequences only 25% of the sequences (A/G/C/T). With the new enzyme combination, I'm not sure how many loci to expect, so I'm not sure whether the selective primer is necessary to further reduce the number of sequences. 


I have to design new P2 adapters for the PstI cutsite. I've left a T after the 4bp overhang, so that I can use the Illumina primers that I already have available: 

PstIFwd
5' AGATCGGAAGAGCACACGTCTGAACTCCAGTCA 3'

PstIRev
5' TGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTGCA 3'





```

```

AMOVA


##Samples to include in first test run


##ddRAD H1 (BVH1) 15 Nov 2015

###Aim

1. Test the new ddRAD protocol (new enzymes & adapters; less samples per lane)

2. Characterise population structure

3. Determine whether individuals from the derived populations can be assigned to the original populations


###Samples

Sam chose samples to represent 3 of the original populations (5-6 individuals each), and several individuals from the derived populations to test whether they can be assigned to the original populations. 

Samples from 2012 had much lower DNA yields than 2015 samples. DNA may have degraded during storage, or swabbing protocols have changed between years. 

Libraries are prepared starting with 12ug of DNA. i.e. 500ng per sample. This equates to at least 5.7ng/ul of DNA made up to 88ul with UHQ. 

The lowest amount of DNA is 400ng, so I will decrease the starting amount to 400ng x 24 = 9.6ug, i.e. I need 4.6ng/ul. 

I will need to concentrate the following samples: 

1. 

2. 12A_IBA_005 (E01/B05) - combined with duplicate extract

3. 12A_IBA_006 (E01/B06) - combined with duplicate extract

4. 12A_IBA_007 (Box3)

5. 12A_IBA_009 (Box3)

6. 12A_UNT_002 (E01/A08) - combined with duplicate extract

7. 12A_UNT_003 (E01/A09) - combined with duplicate extract 

8. 12A_WIL_002 (E01/C08) - combined with duplicate extract

9. 12A_WIL_002 (E01/C08) - combined with duplicate extract

10. 12A_WIL_003 (E01/C09) - combined with duplicate extract

11. 12A_WIL_006 (E01/C12) - combined with duplicate extract



##Qubit of speedvac samples

The concentrations of the speedvac samples are not much higher than before, and the volumes are now <50% of what they were before. A possible explanation is that the DNA is sticking to the sides of the tube, so I am not getting an accurate Qubit measurement. I will use all of the speedVac samples and make all the samples up to 88ul with UHQ. I will reassess after the cleanup & Qubit measurement of Restriction Digested samples. 


21 Nov 2015. 

I ran a restriction digest

Cleanup with Qiagen MinElute columns. 

Qubit HS measurement. 

The SpeedVac concentrated samples have ~1.2-2ng/ul DNA, but the other samples have ~0.6ng/ul DNA. I realise that I didn't adjust the volume needed to 500ng/sample. I need to run another Restriction Digest for 14 samples. 

4 samples need to be concentrated up with the speedvac: 


#BV2: Test with 24 samples

Received HiSeq data back from the FGCZ on 20 Jan 2016. 

Demultiplexed samples using process_radtags from Stacks. 

```
146273859 total sequences;
  47928413 ambiguous barcode drops;
  0 low quality read drops;
  0 ambiguous RAD-Tag drops;
98345446 retained reads.
```

i.e. 32.8% ambiguous barcodes dropped. There was ~20% PhiX added, so ~13% of Bombina data dropped. This is comparable to the R. temporaria data. 


Removal of adapter dimer using Trimmomatic

```
screen -S TrimSubset -L
for i in *.fq; do  java -jar /usr/local/ngseq/src/Trimmomatic-0.33/trimmomatic-0.33.jar SE $i $i.trim ILLUMINACLIP:/usr/local/ngseq/src/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;done
```

QC: nr of reads per sample




#pyRAD BV2: optimise clustering threshold

I will test 90%-98% clustering threshold. I'm expecting an optimum around 93-94%

1. Move demultiplexed and trimmed data to the gdc server. 

2. Set up pyrad 

3. run pyrad on both servers (i.e. 2 runs at a time). 


##pyRAD optimising run results: 

data in pyRADopt_20160217.xlsx

![alt_txt][paralog]
[paralog]:https://cloud.githubusercontent.com/assets/12142475/13322436/ad257048-dbd4-11e5-8b9a-a1150dc4d656.png

![alt_txt][polyfreq]
[polyfreq]:https://cloud.githubusercontent.com/assets/12142475/13322438/ad27b5a6-dbd4-11e5-9b7d-f1b3739659da.png

![alt_txt][nloci]
[nloci]:https://cloud.githubusercontent.com/assets/12142475/13322439/ad29c954-dbd4-11e5-979e-1c677c668d43.png

![alt_txt][Het]
[Het]:https://cloud.githubusercontent.com/assets/12142475/13322638/a680f9d2-dbd5-11e5-9bbc-495b19abcfa4.png


Based on the frequency of polymorphisms (polyfreq), and the number of loci (nloci), I decided on an optimum clustering threshold of 93%. Althought the results based on Het and paralogs are less clear, the other two graphs show that over-splitting seems to happen from ~94-95% clustering threshold. 

As an interesting aside, it seems like the frequency of polymorphisms is only ~30% that of Rana temporaria (0.001 vs 0.003). 



###Analyses

Compare clustering threshold to: 

1. Het (from Pi_ file)

2. Polyfreq (s5_)
 
3. number of loci

4. nr loci filtered as paralogues (s5_)


##Data from training run

FGCZ tested their new illumina machine, and used the BV2 library for the test run. We got data back for Fwd & Rev reads. I will use only the Fwd reads to see what the data looks like, and perhaps add the data to the original BV2 run. 


Demultiplex the data

```
/srv/kenlab/alexjvr_p1795/Bombina/BV2.2$ /usr/local/ngseq/stow/stacks-1.28/bin/process_radtags -i gzfastq -f /srv/gstore4users/p1795/HiSeq_20160121_trainingRun1_DataDelivery/20160121.B-BV2_R1.fastq.gz -o demultiplexed_BV2.2/ -y fastq -b barcodes_BV2 --disable_rad_check -r -D
```


