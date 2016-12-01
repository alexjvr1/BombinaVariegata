#Bombina variegata project  

Full dataset:

72 individuals, 3 HiSeq lanes sequenced

Filtering: 

Step1: Loci genotyped <50% of individuals

rename all samples

SFS & LD

Step2: Heterozygosity (HWE)

SFS & LD

Step3: LD

SFS & LD

Step4: thin 

LD & SFS

step5: MAF0.05

LD & SFS

Checks:  
 
 1. Missingness per population
 
 2. Overall SFS
 
 3. Loci fixed across all pops
 
 4. Loci variable across all pops
 
 5. Overall LD
 
 6. Per pop SFS
 
 7. Per pop LD
 
##Filter

###Step1: Missingness

/Users/alexjvr/2016RADAnalysis/Bombina/BV234/Analyses_20161128

Starting with BV234.vcf

Filter all loci genotyped in <80% of individuals
```
vcftools --vcf BV234.vcf --max-missing 0.8 --recode --recode-INFO-all --out BV234.s1.maxmiss0.8

VCFtools - v0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BV234.vcf
	--recode-INFO-all
	--max-missing 0.8
	--out BV234.s2.maxmiss0.8
	--recode

Eighth Header entry should be INFO: INFO    
After filtering, kept 72 out of 72 Individuals
Outputting VCF file...
After filtering, kept 14046 out of a possible 1065852 Sites
Run Time = 15.00 seconds

```

Rename samples & plot missingness per population
```
bcftools reheader BV234.s1.maxmiss0.8.recode.vcf -s newnames.BV234 --o BV234.s1.names.vcf

vcftools --vcf BV234.s1.names.vcf --missing-indv

##R

library(ggplot2)
BV234.72.s1 <- read.table("out.imiss", header=T)
pop <- read.table("72.popnames", header=T)

BV.s1.df <- as.data.frame(BV234.72.s1)
BV.s1.df.sort <- BV.s1.df[order(BV.s1.df$INDV),]

pop.sort <- pop[order(pop$indiv),]
BV.s1.df.sort$pop <- pop.sort$pop
BV.s1.df.sort$order <- pop.sort$order
head(BV.s1.df.sort)

BV.s1.df.sort <- BV.s1.df.sort[order(BV.s1.df.sort$order),]
BV.s1.df.sort$pop <- factor(BV.s1.df.sort$pop, levels=BV.s1.df.sort$pop)
qplot(pop, F_MISS, data=BV.s1.df.sort, geom=c("boxplot", "jitter"))
```

![alt_txt][missing.s1]
[missing.s1]:https://cloud.githubusercontent.com/assets/12142475/20677691/3dd31034-b594-11e6-83c3-365a89582a14.png


SFS and LD (r2) for the full dataset:
```
##linux

vcftools --vcf BV234.s1.names.vcf --plink --out BV.s1.plink
plink --file BV.s1.plink --recode --recodeA

plink --file BV.s1.plink --freq --out BV.s1
plink --file BV.s1.plink --r2 --out BV.s1

###R

BV.s1.freq <- read.table("BV.s1.frq", header=T)
hist(BV.s1.freq$MAF, main="BV s1 (max.miss0.8; 14046; 72) SFS")

BV.s1.r2 <- read.table("BV.s1.ld", header=T)
hist(BV.s1.r2$R2, main="BV s1 (max.miss0.8; 14046; 72) R2")

```


![alt_txt][SFS.s1]
[SFS.s1]:https://cloud.githubusercontent.com/assets/12142475/20667114/fc08a000-b567-11e6-878a-31246213c243.png

![alt_txt][LD.s1]
[LD.s1]:https://cloud.githubusercontent.com/assets/12142475/20667122/027f4042-b568-11e6-9969-0450396789d7.png


###Step2: identify loci out of HWE

Subset the plink file to all the populations. the --keep file should have two columns, duplicated in this case (FamID, IndID)

```
for i in $(ls plink.pops/); do plink --file BV.s1.plink --keep plink.pops/$i --recode --recodeA --out subset.data/$i.plink ; done
```

And keep only SNPs that are variable in the population (because otherwise the A1 & A2 columns have missing data, and it is impossible to read it into R)
```
for i in $(ls plink.pops/); do plink --file subset.data/$i.plink --maf 0.01 --recode --recodeA --out subset.data/$i.SNPonly; done
```

Calculate HWE within each pop
```
for i in $(ls plink.pops/); do plink --file subset.data/$i.SNPonly --hardy --out $i; done
```
Read into R and
```
BRU.hwe <- read.table("subset.data/1.BRU.hwe", header=T)

ZIN.hwe <- read.table("subset.data/2.ZIN.hwe", header=T)

ZIV.hwe <- read.table("subset.data/3.1.ZIV.hwe", header=T)
NAG.hwe <- read.table("subset.data/3.2.NAG.hwe", header=T)
WIL.hwe <- read.table("subset.data/3.3.WIL.hwe", header=T)
HOP.hwe <- read.table("subset.data/3.4.HOP.hwe", header=T)

CHR.hwe <- read.table("subset.data/4.1.CHR.hwe", header=T)
IBA.hwe <- read.table("subset.data/4.2.IBA.hwe", header=T)
HOC.hwe <- read.table("subset.data/4.3.HOC.hwe", header=T)
UNT.hwe <- read.table("subset.data/4.4.UNT.hwe", header=T)
```

Keep only the Test=All columns
```
BRU.hwe <- subset(BRU.hwe, TEST=="ALL")
ZIN.hwe <- subset(ZIN.hwe, TEST=="ALL")
ZIV.hwe <- subset(ZIV.hwe, TEST=="ALL")
NAG.hwe <- subset(NAG.hwe, TEST=="ALL")
WIL.hwe <- subset(WIL.hwe, TEST=="ALL")
HOP.hwe <- subset(HOP.hwe, TEST=="ALL")
CHR.hwe <- subset(CHR.hwe, TEST=="ALL")
IBA.hwe <- subset(IBA.hwe, TEST=="ALL")
HOC.hwe <- subset(HOC.hwe, TEST=="ALL")
UNT.hwe <- subset(UNT.hwe, TEST=="ALL")
```

Split into loci out of HWE + O.Het > 0.5, and those to keep
```
BRU.hwe.remove <- subset(BRU.hwe, P<0.050001)
BRU.hwe.remove <- rbind(subset(BRU.hwe, O.HET.>0.5))
BRU.hwe.keep <- subset(BRU.hwe, O.HET. <0.51)
BRU.hwe.keep <- subset(BRU.hwe.keep, P>0.05)

ZIN.hwe.remove <- subset(ZIN.hwe, P<0.050001)
ZIN.hwe.remove <- rbind(subset(ZIN.hwe, O.HET.>0.5))
ZIN.hwe.keep <- subset(ZIN.hwe, O.HET. <0.51)
ZIN.hwe.keep <- subset(ZIN.hwe.keep, P>0.05)

ZIV.hwe.remove <- subset(ZIV.hwe, P<0.050001)
ZIV.hwe.remove <- rbind(subset(ZIV.hwe, O.HET.>0.5))
ZIV.hwe.keep <- subset(ZIV.hwe, O.HET. <0.51)
ZIV.hwe.keep <- subset(ZIV.hwe.keep, P>0.05)

NAG.hwe.remove <- subset(NAG.hwe, P<0.050001)
NAG.hwe.remove <- rbind(subset(NAG.hwe, O.HET.>0.5))
NAG.hwe.keep <- subset(NAG.hwe, O.HET. <0.51)
NAG.hwe.keep <- subset(NAG.hwe.keep, P>0.05)

WIL.hwe.remove <- subset(WIL.hwe, P<0.050001)
WIL.hwe.remove <- rbind(subset(WIL.hwe, O.HET.>0.5))
WIL.hwe.keep <- subset(WIL.hwe, O.HET. <0.51)
WIL.hwe.keep <- subset(WIL.hwe.keep, P>0.05)

HOP.hwe.remove <- subset(HOP.hwe, P<0.050001)
HOP.hwe.remove <- rbind(subset(HOP.hwe, O.HET.>0.5))
HOP.hwe.keep <- subset(HOP.hwe, O.HET. <0.51)
HOP.hwe.keep <- subset(HOP.hwe.keep, P>0.05)

CHR.hwe.remove <- subset(CHR.hwe, P<0.050001)
CHR.hwe.remove <- rbind(subset(CHR.hwe, O.HET.>0.5))
CHR.hwe.keep <- subset(CHR.hwe, O.HET. <0.51)
CHR.hwe.keep <- subset(CHR.hwe.keep, P>0.05)

IBA.hwe.remove <- subset(IBA.hwe, P<0.050001)
IBA.hwe.remove <- rbind(subset(IBA.hwe, O.HET.>0.5))
IBA.hwe.keep <- subset(IBA.hwe, O.HET. <0.51)
IBA.hwe.keep <- subset(IBA.hwe.keep, P>0.05)

HOC.hwe.remove <- subset(HOC.hwe, P<0.050001)
HOC.hwe.remove <- rbind(subset(HOC.hwe, O.HET.>0.5))
HOC.hwe.keep <- subset(HOC.hwe, O.HET. <0.51)
HOC.hwe.keep <- subset(HOC.hwe.keep, P>0.05)

UNT.hwe.remove <- subset(UNT.hwe, P<0.050001)
UNT.hwe.remove <- rbind(subset(UNT.hwe, O.HET.>0.5))
UNT.hwe.keep <- subset(UNT.hwe, O.HET. <0.51)
UNT.hwe.keep <- subset(UNT.hwe.keep, P>0.05)

```


Bind all files together and plot the frequency of SNP names
```
BV.hwe.remove.file <- do.call(rbind, lapply(ls(pattern="remove$"), get))

test.table <- data.frame(table(BV.hwe.remove.file$SNP))
hist(test.table$Freq, xlab="Number of pops in which locus occurs", ylab="Frequency", main="Frequency of SNPs deviating from HWE and H.Obs>0.5", breaks=10, xlim=c(0,10), xaxt='n')
axis(side=1, at=seq(0.5,9.5, 1), labels=seq(1,10,1))
```

![alt_txt][HWE.s1]
[HWE.s1]:https://cloud.githubusercontent.com/assets/12142475/20668456/4c06ef10-b56f-11e6-92e3-c5f76367d3db.png


Sturgeon paper removes loci that deviate from HWE in >60% of the population. I will us 5pops as the cutoff
```
HWE.loci.remove.names <- subset(test.table, Freq>4)
write.table(HWE.loci.remove.names$Var1, quote=F, col.names=F, file="HWE.loci.remove.names", row.names=F)
```

And remove the 413 loci using plink
```
plink --file BV.s1.plink --exclude HWE.loci.remove.names --recode --recodeA --out BV.s2
```


###Step3: identify loci with high LD

On the same s1 plink file, calculate the per population LD & SFS

```
for i in $(ls plink.pops/); do plink --file subset.data/$i.plink --r2 --out subset.data/$i; done
for i in $(ls plink.pops/); do plink --file subset.data/$i.plink --freq --out subset.data/$i; done
```
R

```
BRU.freq <- read.table("subset.data/1.BRU.frq", header=T)
ZIN.freq <- read.table("subset.data/2.ZIN.frq", header=T)
ZIV.freq <- read.table("subset.data/3.1.ZIV.frq", header=T)
NAG.freq <- read.table("subset.data/3.2.NAG.frq", header=T)
WIL.freq <- read.table("subset.data/3.3.WIL.frq", header=T)
HOP.freq <- read.table("subset.data/3.4.HOP.frq", header=T)
CHR.freq <- read.table("subset.data/4.1.CHR.frq", header=T)
IBA.freq <- read.table("subset.data/4.2.IBA.frq", header=T)
HOC.freq <- read.table("subset.data/4.3.HOC.frq", header=T)
UNT.freq <- read.table("subset.data/4.4.UNT.frq", header=T)

par(mfrow=c(4,3))
my.bin.width <- 0.05
hist(BRU.freq$MAF, main="BRU s1 (max.miss0.8; 13886; 7) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(ZIN.freq$MAF, main="ZIN s1 (max.miss0.8; 13886; 8) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(ZIV.freq$MAF, main="ZIV s1 (max.miss0.8; 13886; 7) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(NAG.freq$MAF, main="NAG s1 (max.miss0.8; 13886; 6) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(WIL.freq$MAF, main="WIL s1 (max.miss0.8; 13886; 6) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(HOP.freq$MAF, main="HOP s1 (max.miss0.8; 13886; 6) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(CHR.freq$MAF, main="CHR s1 (max.miss0.8; 13886; 8) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(IBA.freq$MAF, main="IBA s1 (max.miss0.8; 13886; 8) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(HOC.freq$MAF, main="HOC s1 (max.miss0.8; 13886; 10) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(UNT.freq$MAF, main="UNT s1 (max.miss0.8; 13886; 6) SFS", breaks=seq(0,0.5, by=my.bin.width))


BRU.ld <- read.table("subset.data/1.BRU.ld", header=T)
ZIN.ld <- read.table("subset.data/2.ZIN.ld", header=T)
ZIV.ld <- read.table("subset.data/3.1.ZIV.ld", header=T)
NAG.ld <- read.table("subset.data/3.2.NAG.ld", header=T)
WIL.ld <- read.table("subset.data/3.3.WIL.ld", header=T)
HOP.ld <- read.table("subset.data/3.4.HOP.ld", header=T)
CHR.ld <- read.table("subset.data/4.1.CHR.ld", header=T)
IBA.ld <- read.table("subset.data/4.2.IBA.ld", header=T)
HOC.ld <- read.table("subset.data/4.3.HOC.ld", header=T)
UNT.ld <- read.table("subset.data/4.4.UNT.ld", header=T)

par(mfrow=c(4,3))
my.bin.width <- 0.05
hist(BRU.ld$R2, main="BRU s1 (max.miss0.8; 13886; 7) LD", breaks=seq(0,1, by=my.bin.width))
hist(ZIN.ld$R2, main="ZIN s1 (max.miss0.8; 13886; 8) LD", breaks=seq(0,1, by=my.bin.width))
hist(ZIV.ld$R2, main="ZIV s1 (max.miss0.8; 13886; 7) LD", breaks=seq(0,1, by=my.bin.width))
hist(NAG.ld$R2, main="NAG s1 (max.miss0.8; 13886; 6) LD", breaks=seq(0,1, by=my.bin.width))
hist(WIL.ld$R2, main="WIL s1 (max.miss0.8; 13886; 6) LD", breaks=seq(0,1, by=my.bin.width))
hist(HOP.ld$R2, main="HOP s1 (max.miss0.8; 13886; 6) LD", breaks=seq(0,1, by=my.bin.width))
hist(CHR.ld$R2, main="CHR s1 (max.miss0.8; 13886; 8) LD", breaks=seq(0,1, by=my.bin.width))
hist(IBA.ld$R2, main="IBA s1 (max.miss0.8; 13886; 8) LD", breaks=seq(0,1, by=my.bin.width))
hist(HOC.ld$R2, main="HOC s1 (max.miss0.8; 13886; 10) LD", breaks=seq(0,1, by=my.bin.width))
hist(UNT.ld$R2, main="UNT s1 (max.miss0.8; 13886; 6) LD", breaks=seq(0,1, by=my.bin.width))


##Check frequency of loci iwth R2>0.8

##Keep only loci R2>0.8 and plot frequency in R
BRU.ld.remove <- subset(BRU.ld, R2>0.799999)
ZIN.ld.remove <- subset(ZIN.ld, R2>0.799999)
ZIV.ld.remove <- subset(ZIV.ld, R2>0.799999)
NAG.ld.remove <- subset(NAG.ld, R2>0.799999)
WIL.ld.remove <- subset(WIL.ld, R2>0.799999)
HOP.ld.remove <- subset(HOP.ld, R2>0.799999)
CHR.ld.remove <- subset(CHR.ld, R2>0.799999)
IBA.ld.remove <- subset(IBA.ld, R2>0.799999)
HOC.ld.remove <- subset(HOC.ld, R2>0.799999)
UNT.ld.remove <- subset(UNT.ld, R2>0.799999)

BV.ld.remove.file <- do.call(rbind, lapply(ls(pattern="ld.remove$"), get))
test.table.ld <- data.frame(table(BV.ld.remove.file$SNP_A, BV.ld.remove.file$SNP_B))
test.table.ld.subset <- subset(test.table.ld, Freq>0)
test.table.ld.1.subset <- subset(test.table.ld, Freq>1)

LD.loci.remove.names <- subset(test.table.ld, Freq>4)
write.table(LD.loci.remove.names$Var1, quote=F, col.names=F, file="LD.loci.remove.names", row.names=F)

par(mfrow=c(1,2))
hist(test.table.ld.subset$Freq, main="Frequency of loci R2>0.8 across 10 pops")
hist(test.table.ld.1.subset$Freq, main="Frequency of loci R2>0.8 across 10 pops")
```

SFS

![alt_txt][SFS.s1]
[SFS.s1]:https://cloud.githubusercontent.com/assets/12142475/20670053/6f9fd0a6-b577-11e6-882f-3eeee155f17c.png


LD

![alt_txt][LD.s1]
[LD.s1]:https://cloud.githubusercontent.com/assets/12142475/20670052/6f9f9f64-b577-11e6-9734-6070d03753e2.png


Freq of loci with R2>0.8
![alt_txt][LD.freq.s1]
[LD.freq.s1]:https://cloud.githubusercontent.com/assets/12142475/20670267/7cfeb72a-b578-11e6-8f0c-6b17e014c727.png


##Step 4: thin

Remove multiple SNPs per locus. 


```
vcftools --vcf BV234.s1.names.vcf --thin 200 --recode --recode-INFO-all --out BV.s4

VCFtools - v0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BV234.s1.names.vcf
	--recode-INFO-all
	--thin 200
	--out BV.s4
	--recode

After filtering, kept 72 out of 72 Individuals
Outputting VCF file...
After filtering, kept 8740 out of a possible 14046 Sites
Run Time = 1.00 seconds

```

###Step5: MAF 0.05

low MAF can bias pop structure and Fst outlier analyses. I will filter for a minimum MAF of 0.05 across the entire dataset. This is a minimum minor allele count of 4 in 72 individuals. i.e. this is also a filter for any spurious SNP calls that may have been kept through the filtering process thus far. 
```
vcftools --vcf BV.s4.recode.vcf --maf 0.05 --recode --recode-INFO-all --out BV.s5.maf0.05

VCFtools - v0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf BV.s4.recode.vcf
	--recode-INFO-all
	--maf 0.05
	--out BV.s5.maf0.05
	--recode

After filtering, kept 72 out of 72 Individuals
Outputting VCF file...
After filtering, kept 1931 out of a possible 8740 Sites
Run Time = 0.00 seconds

```


##Final Checks
Convert to plink and check that all s2 & s3 loci have been removed
```
vcftools --vcf BV.s5.maf0.05.recode.vcf --plink --out BV.s5.plink

plink --file BV.s5.plink --exclude HWE.loci.remove.names --recode --recodeA --out BV.s5.hwe.ld
plink --file BV.s5.hwe --exclude LD.loci.remove.names --recode --recodeA --out BV.s5.hwe.ld
```

266 loci removed for deviation from HWE, 0 for LD

Final dataset: 

71 individuals  (see below, 1 individual is removed)

10 populations

86.9% genotyping rate

1665 loci


Convert .ped to .vcf using pgdspider (sequences)

Missingness per population
```
vcftools --vcf BV.s5.vcf --missing-indv --out BV.s5

#R
library(ggplot2)
BV.s5 <- read.table("BV.s5.imiss", header=T)
BV.s5.sort <- BV.s5[order(BV.s5$INDV),]
pop <- read.table("72.popnames.alphabetical", header=T)
BV.s5.sort$pop <- pop$pop
BV.s5.sort$pop.order <- pop$order

BV.s5.sort2 <- BV.s5.sort[order(BV.s5.sort$pop.order),]

BV.s5.sort2$pop <- factor(BV.s5.sort2$pop, levels=BV.s5.sort2$pop)   ##sort pop.nr. Numbers from south to North

qplot(pop, F_MISS, data=BV.s5.sort2, geom=c("boxplot", "jitter"))
```

And remove individuals with >45% missing data

Only IBA.05 needs to be removed. 
```
vcftools --vcf BV.s5.vcf --remove lowDP.indiv --recode --recode-INFO-all --out BV.71.1665
```

![alt.txt][final.missing]
[final.missing]:https://cloud.githubusercontent.com/assets/12142475/20675736/99feff00-b58d-11e6-87e4-b2ecfb6cf0bf.png


SFS & LD per population
```
for i in $(ls plink.pops/); do plink --file BV.s5.hwe.ld --keep plink.pops/$i --recode --recodeA --out final.subsets/$i.1665; done

for i in $(ls plink.pops/); do plink --file final.subsets/$i.1665 --freq --out final.subsets/$i; done
for i in $(ls plink.pops/); do plink --file final.subsets/$i.1665 --r2 --out final.subsets/$i; done

BRU.freq <- read.table("final.subsets/1.BRU.frq", header=T)
ZIN.freq <- read.table("final.subsets/2.ZIN.frq", header=T)
ZIV.freq <- read.table("final.subsets/3.1.ZIV.frq", header=T)
NAG.freq <- read.table("final.subsets/3.2.NAG.frq", header=T)
WIL.freq <- read.table("final.subsets/3.3.WIL.frq", header=T)
HOP.freq <- read.table("final.subsets/3.4.HOP.frq", header=T)
CHR.freq <- read.table("final.subsets/4.1.CHR.frq", header=T)
IBA.freq <- read.table("final.subsets/4.2.IBA.frq", header=T)
HOC.freq <- read.table("final.subsets/4.3.HOC.frq", header=T)
UNT.freq <- read.table("final.subsets/4.4.UNT.frq", header=T)

par(mfrow=c(4,3))
my.bin.width <- 0.05
hist(BRU.freq$MAF, main="BRU final (1668; 7) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(ZIN.freq$MAF, main="ZIN final (1668; 8) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(ZIV.freq$MAF, main="ZIV final (1668; 7) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(NAG.freq$MAF, main="NAG final (1668; 6) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(WIL.freq$MAF, main="WIL final (1668; 6) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(HOP.freq$MAF, main="HOP final (1668; 6) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(CHR.freq$MAF, main="CHR final (1668; 8) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(IBA.freq$MAF, main="IBA final (1668; 7) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(HOC.freq$MAF, main="HOC final (1668; 10) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(UNT.freq$MAF, main="UNT final (1668; 6) SFS", breaks=seq(0,0.5, by=my.bin.width))


BRU.ld <- read.table("final.subsets/1.BRU.ld", header=T)
ZIN.ld <- read.table("final.subsets/2.ZIN.ld", header=T)
ZIV.ld <- read.table("final.subsets/3.1.ZIV.ld", header=T)
NAG.ld <- read.table("final.subsets/3.2.NAG.ld", header=T)
WIL.ld <- read.table("final.subsets/3.3.WIL.ld", header=T)
HOP.ld <- read.table("final.subsets/3.4.HOP.ld", header=T)
CHR.ld <- read.table("final.subsets/4.1.CHR.ld", header=T)
IBA.ld <- read.table("final.subsets/4.2.IBA.ld", header=T)
HOC.ld <- read.table("final.subsets/4.3.HOC.ld", header=T)
UNT.ld <- read.table("final.subsets/4.4.UNT.ld", header=T)

par(mfrow=c(4,3))
my.bin.width <- 0.05
hist(BRU.ld$R2, main="BRU final (1668; 7) LD", breaks=seq(0,1, by=my.bin.width))
hist(ZIN.ld$R2, main="ZIN final (1668; 8) LD", breaks=seq(0,1, by=my.bin.width))
hist(ZIV.ld$R2, main="ZIV final (1668; 7) LD", breaks=seq(0,1, by=my.bin.width))
hist(NAG.ld$R2, main="NAG final (1668; 6) LD", breaks=seq(0,1, by=my.bin.width))
hist(WIL.ld$R2, main="WIL final (1668; 6) LD", breaks=seq(0,1, by=my.bin.width))
hist(HOP.ld$R2, main="HOP final (1668; 6) LD", breaks=seq(0,1, by=my.bin.width))
hist(CHR.ld$R2, main="CHR final (1668; 8) LD", breaks=seq(0,1, by=my.bin.width))
hist(IBA.ld$R2, main="IBA final (1668; 7) LD", breaks=seq(0,1, by=my.bin.width))
hist(HOC.ld$R2, main="HOC final (1668; 10) LD", breaks=seq(0,1, by=my.bin.width))
hist(UNT.ld$R2, main="UNT final (1668; 6) LD", breaks=seq(0,1, by=my.bin.width))
```
![alt_txt][SFS.final]
[SFS.final]:https://cloud.githubusercontent.com/assets/12142475/20675734/99feba54-b58d-11e6-8e07-40dda1b717f6.png

![alt_txt][LD.final]
[LD.final]:https://cloud.githubusercontent.com/assets/12142475/20675735/99feae24-b58d-11e6-9312-5f1383a1627b.png




#Summary statistics

Nr fixed loci per population

```
#And with only variable loci: 

BRU.frq.var <- subset(BRU.freq, MAF>0.001)
ZIN.frq.var <- subset(ZIN.freq, MAF>0.001)
ZIV.frq.var <- subset(ZIV.freq, MAF>0.001)
NAG.frq.var <- subset(NAG.freq, MAF>0.001)
WIL.frq.var <- subset(WIL.freq, MAF>0.001)
HOP.frq.var <- subset(HOP.freq, MAF>0.001)
CHR.frq.var <- subset(CHR.freq, MAF>0.001) 
IBA.frq.var <- subset(IBA.freq, MAF>0.001) 
HOC.frq.var <- subset(HOC.freq, MAF>0.001) 
UNT.frq.var <- subset(UNT.freq, MAF>0.001) 


##Loci variable in x populations
BV.var.loci.freq <- do.call(rbind, lapply(ls(pattern=".frq.var$"), get))
summary(BV.var.loci.freq)
BV.var.loci.freq.keep <- data.frame(table(BV.var.loci.freq$SNP)) 
#SE.var.loci.freq.keep <- subset(SE.var.loci.freq.keep, Freq>0) 
hist(BV.var.loci.freq.keep$Freq, xlab="Nr pops", ylab="Frequency", main="BV: Frequency of variable loci across pops", breaks=seq(0.5,9.5, by=1.0))
```

Fixed loci

```
###Looking only at fixed loci

BRU.frq.fix <- subset(BRU.freq, MAF<0.001)  ##399 loci
ZIN.frq.fix <- subset(ZIN.freq, MAF<0.001)  ##409
ZIV.frq.fix <- subset(ZIV.freq, MAF<0.001)  ##380
NAG.frq.fix <- subset(NAG.freq, MAF<0.001)  ##385
WIL.frq.fix <- subset(WIL.freq, MAF<0.001)  ##341
HOP.frq.fix <- subset(HOP.freq, MAF<0.001)  ##429
CHR.frq.fix <- subset(CHR.freq, MAF<0.001) ##346
IBA.frq.fix <- subset(IBA.freq, MAF<0.001)  ##476
HOC.frq.fix <- subset(HOC.freq, MAF<0.001)  ##289
UNT.frq.fix <- subset(UNT.freq, MAF<0.001)  ##457



BV.region.frq.fix.table <- do.call(rbind, lapply(ls(pattern="fix$"), get))

BV.fix.region.table.keep <- data.frame(table(BV.region.frq.fix.table$SNP))
BV.fix.region.table.keep <- subset(BV.fix.region.table.keep, Freq>0)
my.bin.width=1
hist(BV.fix.region.table.keep$Freq, xlab="Number of pops", ylab="Frequency", main="BV: Frequency of fixed loci found in increasing number of populations", breaks=seq(0.5,9.5, by=1))
```

![alt_txt][fixed.loci]
[fixed.loci]:https://cloud.githubusercontent.com/assets/12142475/20677295/e8f0f4ec-b592-11e6-996b-190caf37656f.png

![alt_txt][variable.loci]



##AvgHet

Based on the .phy output file on gdcsrv1. 
Looking at the .loci file, the data looks to have assembled fairly well. 

This estimate should be fairly accurate, since it is on a per individual basis, and excludes all Ns and missing data (dashes). Although it might be slightly higher than expected since singletons are still included. 

```
#Missing data
tr -d -c 'N\n'< BV234.phy |awk '{print length; }'
##gaps
sed 's/[^-]//g' BV234.phy |awk '{print length}'

#Transitions
tr -d -c 'R\n'< BV234.phy |awk '{print length; }'
tr -d -c 'Y\n'< BV234.phy |awk '{print length; }'

#Transversions
tr -d -c 'S\n'< BV234.phy |awk '{print length; }'
tr -d -c 'W\n'< BV234.phy |awk '{print length; }'
tr -d -c 'K\n'< BV234.phy |awk '{print length; }'
tr -d -c 'M\n'< BV234.phy |awk '{print length; }'

##Print out all the names in the phylip file
grep -Eo '^[^ ]+' BV234.phy 


```




##Fst
```
BV.71 <- read.structure("BV.71.1665.str")
BV.71

pop.BV.71 <- read.table("pops4pgdspider", header=F)  ##make sure the populations are numbered "01.DE.B", etc.
BV.10pops.factor <- as.factor(pop.BV.71$V2)
BV.71@pop <- BV.10pops.factor

hier.BV <- genind2hierfstat(BV.71)

BV.fst <- pairwise.fst(BV.71, pop=NULL, res.type=c("dist", "matrix"))

m <- as.matrix(BV.fst)
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2)<- c("c1","c2", "distance")

library(gplots)

shadesOfGrey <- colorRampPalette(c("grey100", "grey0"))  ##define the colourpalette. 

Dend <- read.table("heatmap.popcolours", header=T)  ##list of colour names for each population based on R colour palatte. In alphabetical order (as in genind file)
Dend.Colours <- as.character(Dend$colours.pop)

par(oma=c(1,1,2,1))
heatmap.2(as.matrix(BV.fst), trace="none", RowSideColors=Dend.Colours, ColSideColors=Dend.Colours, col=shadesOfGrey, labRow=F, labCol=F, key.ylab=NA, key.xlab=NA, key.title="Fst Colour Key", keysize=0.9, main="Pairwise Fst and dendrograms of BV: 10pops, 3regions, 1665loci")  ##RowSideColors is for the dendrogram on the row, ColSideColors for the upper dendrogram. Colour order should be the same as the input. The pop order is alphabetical in the output. 
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

popnames.all <- as.character(c("BRU", "ZIN", "Middle", "East"))
legend("bottom", popnames.all, xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty="o", pch=15, col=c("darkorange3", "darkorange2", "darkorchid1", "deepskyblue1"), title="Regions")

```
![alt_txt][Fst.all]
[Fst.all]:https://cloud.githubusercontent.com/assets/12142475/20679243/5c7c4bc6-b59a-11e6-93fa-3e54711291e9.png


Frequency of per locus Fst
```
BV.perloc <- stats.BV71$perloc
qplot(BV71.perloc$Fst, geom="histogram", binwidth=0.005)
```

![alt_txt][Fst.freq]
[Fst.freq]:https://cloud.githubusercontent.com/assets/12142475/20748674/656a241e-b6ef-11e6-8054-9654d0211e83.png


Outlier locus:
```
BV.highFst <- subset(BV71.perloc, Fst>0.4)
BV.highFst
            Ho     Hs     Ht    Dst  Htp   Dstp    Fst   Fstp Fis   Dest
X240046.110  0 0.0582 0.1138 0.0556 0.12 0.0618 0.4887 0.5151   1 0.0656
```


##IBD

Fst/(1-Fst) vs log.dist(km) - according to Rousset et al. 1997, this is correlated with the effective population density (De) and effective dispersal distance (variance)

```
##Fst 15 pops -> Fst/(1-Fst)

library(reshape)
library(fields)


m <- as.matrix(SEall.fst)
m
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2) <- c("c1", "c2", "distance")
m2
m2$IBD <- m2$distance/(1-m2$distance)


BV.pop.coords <- read.table("BV.pop.coords", header=T)
BVpop_lon.lat <- cbind(BV.pop.coords$Long, BV.pop.coords$Lat)
distance.matrix.BVpop <- rdist.earth(BVpop_lon.lat, miles=F)  ##great circle dist based on the coordinates
m.dist <- as.matrix(distance.matrix.BVpop)
summary(m.dist)

m2.dist <- melt(m.dist)[melt(upper.tri(m.dist))$value,]
names(m2.dist) <- c("c1", "c2", "distance")
m2.dist$m <- (m2.dist$distance*1000)
summary(m2.dist)
m2.dist$log.m <- log(m2.dist$m)


library(MASS)
#dens <- kde2d(m2$IBD,m2.dist$log.m, n=10)
#myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(m2$IBD~m2.dist$log.m, pch=20,cex=.5, xlab="log Geographic distance (m)", ylab="Fst/(1-Fst)")
#image(dens, col=transp(myPal(10),.7), add=TRUE)
abline(fit <- lm(m2$IBD~m2.dist$log.m))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - BVall")
```

And within Middle and East: 

```
##Middle

m2.middle <- subset(m2, c1>2)
m2.middle <- subset(m2.middle, c1<7)
m2.middle <- subset(m2.middle, c2<7)
m2.middle <- subset(m2.middle, c2>2)

m2.dist.middle <- subset(m2.dist, c2<7)
m2.dist.middle <- subset(m2.dist.middle, c2<7)
m2.dist.middle <- subset(m2.dist.middle, c2>2)
m2.dist.middle <- subset(m2.dist.middle, c1>2)

library(MASS)
plot(m2.middle$IBD~m2.dist.middle$log.m, pch=20,cex=.5, xlab="log Geographic distance (m)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.middle$IBD~m2.dist.middle$log.m))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - BVmiddle")


##EAST
m2.east <- subset(m2, c1>6)
m2.east <- subset(m2.east, c2>6)

m2.dist.east <- subset(m2.dist, c2>6)
m2.dist.east <- subset(m2.dist.east, c1>6)

plot(m2.east$IBD~m2.dist.east$log.m, pch=20,cex=.5, xlab="log Geographic distance (m)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.east$IBD~m2.dist.east$log.m))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - BVeast")


##all 3 plots
par(mfrow=c(2,2))
plot(m2.middle$IBD~m2.dist.middle$log.m, pch=20,cex=.5, xlab="log Geographic distance (m)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.middle$IBD~m2.dist.middle$log.m))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - BVmiddle")

plot(m2.east$IBD~m2.dist.east$log.m, pch=20,cex=.5, xlab="log Geographic distance (m)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.east$IBD~m2.dist.east$log.m))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - BVeast")


plot(m2$IBD~m2.dist$log.m, pch=20,cex=.5, xlab="log Geographic distance (m)", ylab="Fst/(1-Fst)")
#image(dens, col=transp(myPal(10),.7), add=TRUE)
abline(fit <- lm(m2$IBD~m2.dist$log.m))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - BVall")
```

![alt_txt][IBD.BV]
[IBD.BV]:https://cloud.githubusercontent.com/assets/12142475/20680390/2949b004-b59f-11e6-9d88-8b810c5494af.png




##AMOVA



##DAPC
tutorial-dapc: A tutorial for Discriminant Analysis of Principal Components (DAPC) using adegenet 2.0.0

total variance = (variance between groups) + (variance within groups)

or more simply, denoting X the data matrix:

VAR(X) = B(X) +W(X)

Usual approaches such as Principal Component Analysis (PCA) or Principal Coordinates Analysis (PCoA / MDS) focus on V AR(X). That is, they only describe the global diversity, possibly overlooking differences between groups. On the contrary, DAPC optimizes B(X) while minimizing W(X): it seeks synthetic variables, the discriminant functions, which show differences between groups as best as possible while minimizing variation within clusters.

##1. estimate the number of clusters

Using k-means. Which finds the number of clusters with minimises W(X) and maximises B(X). Compare using BIC

Run algorithm on PCA transformed data. I.e. reduce the dataset so that it can run faster. 

```
grp.BVall <- find.clusters(BV.71, max.n.clust=40)

> choose nr of PCs: 200  ##I try to keep all the PCs

> choose k: 5 ##see figure below


names.15 <- c("DE.B", "DE.K", "DE.W", "Sk.SF", "SK.SL", "Upp.Gra", "Upp.K", "Upp.O", "Um.Gr", "Um.Taf", "LT1", "LT3", "Kir.G", "Kir.L", "FIN")
names.15 <- as.character(names.15)
table.value(table(pop(BV.71), grp.BVall$grp), col.lab=paste("inf", 1:6), row.lab=1:10)

dapc1.BVall <- dapc(BV.71, grp.BVall$grp)
scatter(dapc1.BVall)
```




##PCA

PCAdapt in R:
```
##convert .vcf to plink 
##linux

vcftools --vcf BV.71.1665.FINAL.vcf --plink --out BV.71.plink

plink --file BV.71.plink --recode --recodeA

##R
library(pcadapt)

BV.71 <- read.pcadapt("BV.71.plink.ped", type="ped")
Summary:

        - input file      BV.71.plink.ped
        - output file     BV.71.plink.pcadapt

	- number of individuals detected:	71
	- number of loci detected:		1665

File has been sucessfully converted.


##Check the nr of PCs

x <- pcadapt(BV.71, K=20)

Reading file BV.71.plink.pcadapt...
Number of SNPs: 1665
Number of individuals: 71
Number of SNPs with minor allele frequency lower than 0.05 ignored: 1
14320 out of 118215 missing data ignored.

plot(x, option="screeplot")  ##PC for pop structure = on the steep curve
```
![alt_txt][screeplot]
[screeplot]:https://cloud.githubusercontent.com/assets/12142475/20680600/078e327c-b5a0-11e6-90f8-06da3ee281c9.png



Plot the PCA using population information
```
popnames <- read.table("popnames", header=T)
head(popnames)
                  indiv     pop region
1 NAG.09_12A_NAG.09_12A 3.2.NAG Middle
2 WIL.15_12A_WIL.15_12A 3.3.WIL Middle
3 WIL.19_12A_WIL.19_12A 3.3.WIL Middle
4 BRU.08_12B_BRU.08_12B   1.BRU    BRU
5 NAG.01_12B_NAG.01_12B 3.2.NAG Middle
6 CHR.01_15A_CHR.01_15A 4.1.CHR   EAST

poplist <- popnames[,2]

plot(x,option="scores",pop=poplist)

poplist <- as.character(popnames[,3]) ##select regions
plot(x,option="scores",pop=poplist)
```




##hierarchical fastStructure

K1-10, 5 reps each

Run fastStructure on gdcserver

First convert to PLINK. This way it is easier to keep track of the sample order. 

```
rawdataAnalysis/BV.71.1665.recode.vcf fastSTRUCTURE/
cd fastSTRUCTURE/
vcftools --vcf BV.71.1665.recode.vcf --plink --out BV.71.1665.plink
plink --file BV.71.1665.plink --recode --recodeA --out BV.71.1665.plink
plink --file BV.71.1665.plink --out BV.71.1665 --make-bed
scp -r BV.71.1665.FULLdata/ alexjvr@gdcsrv1.ethz.ch:/gdc_home4/alexjvr/Bombina/Analysis.BV234/fastStructure/
```


On GDCserver: 


```
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 2 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K2.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 2 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K2.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 2 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K2.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 2 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K2.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 2 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K2.5

/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 3 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K3.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 3 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K3.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 3 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K3.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 3 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K3.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 3 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K3.5

/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 4 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K4.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 4 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K4.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 4 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K4.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 4 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K4.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 4 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K4.5

/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 5 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K5.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 5 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K5.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 5 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K5.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 5 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K5.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 5 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K5.5

/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 6 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K6.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 6 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K6.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 6 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K6.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 6 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K6.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 6 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K6.5

/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 7 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K7.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 7 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K7.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 7 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K7.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 7 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K7.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 7 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K7.5

/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 8 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K8.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 8 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K8.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 8 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K8.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 8 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K8.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 8 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K8.5

/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 9 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K9.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 9 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K9.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 9 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K9.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 9 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K9.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 9 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K9.5

/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 10 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K10.1
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 10 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K10.2
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 10 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K10.3
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 10 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K10.4
/usr/bin/python2.6 /usr/local/fastStructure-20150714/structure.py -K 10 --format=bed --input=BV.71.1665 --output=K2-10.plinkinput/BV.K10.5
```

```
/usr/bin/python2.6 /usr/local/fastStructure-20150714/chooseK.py --input=K2-10.plinkinput/BV.K*
Model complexity that maximizes marginal likelihood = 2
Model components used to explain structure in data = 6
```

Summarise structure plots with CLUMPP

First, create a paramfile. This can be a copy of the example file. Specify Datatype 0 for individuals (vs 1 for pops), and number of individuals. Check the example from my run or from the example input for the chicken dataset. Most of the other CLUMPP parameters can be specified via command line, which will override the paramfile.

indfiles need to be created: paste all the structure outputs for a specific K below each other. No headers.
```
cat CHall.Data1_K2.1.2.meanQ CHall.Data1_K2.2.2.meanQ CHall.Data1_K2.3.2.meanQ CHall.Data1_K2.4.2.meanQ CHall.Data1_K2.5.2.meanQ > K2.Q.meanQ
cat CHall.Data1_K3.1.3.meanQ CHall.Data1_K3.2.3.meanQ CHall.Data1_K3.3.3.meanQ CHall.Data1_K3.4.3.meanQ CHall.Data1_K3.5.3.meanQ > K3.Q.meanQ
cat CHall.Data1_K4.1.4.meanQ CHall.Data1_K4.2.4.meanQ CHall.Data1_K4.3.4.meanQ CHall.Data1_K4.4.4.meanQ CHall.Data1_K4.5.4.meanQ > K4.Q.meanQ
cat CHall.Data1_K5.1.5.meanQ CHall.Data1_K5.2.5.meanQ CHall.Data1_K5.3.5.meanQ CHall.Data1_K5.4.5.meanQ CHall.Data1_K5.5.5.meanQ > K5.Q.meanQ
cat CHall.Data1_K6.1.6.meanQ CHall.Data1_K6.2.6.meanQ CHall.Data1_K6.3.6.meanQ CHall.Data1_K6.4.6.meanQ CHall.Data1_K6.5.6.meanQ > K6.Q.meanQ
cat CHall.Data1_K7.1.7.meanQ CHall.Data1_K7.2.7.meanQ CHall.Data1_K7.3.7.meanQ CHall.Data1_K7.4.7.meanQ CHall.Data1_K7.5.7.meanQ > K7.Q.meanQ
cat CHall.Data1_K8.1.8.meanQ CHall.Data1_K8.2.8.meanQ CHall.Data1_K8.3.8.meanQ CHall.Data1_K8.4.8.meanQ CHall.Data1_K8.5.8.meanQ > K8.Q.meanQ
cat CHall.Data1_K9.1.9.meanQ CHall.Data1_K9.2.9.meanQ CHall.Data1_K9.3.9.meanQ CHall.Data1_K9.4.9.meanQ CHall.Data1_K9.5.9.meanQ > K9.Q.meanQ
cat CHall.Data1_K10.1.10.meanQ CHall.Data1_K10.2.10.meanQ CHall.Data1_K10.3.10.meanQ CHall.Data1_K10.4.10.meanQ CHall.Data1_K10.5.10.meanQ > K10.Q.meanQ
```




I've specified Greedy option 2

```
CLUMPP paramfile -k 2 -i K2.Q.indfile.txt -o K2.meanQ.out
```


###Hierarchy 1: 

Structure on only Middle and East seperately: 
Mid
```
vcftools --vcf BV.71.1665.FINAL.vcf --keep BV.Middle/Middle.names --out BV.Middle/BV.25.Mid.vcf

```

East
```
vcftools --vcf BV.71.1665.FINAL.vcf --keep BV.East/East.names --out BV.East/BV.31.East.vcf
```

##TESS3

K1-10
