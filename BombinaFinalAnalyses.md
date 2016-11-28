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
BV234.72.s1$pop <- pop$pop
BV234.72.s1$pop.order <- pop$order

BV234.72.s1.sort <- BV234.72.s1[order(BV234.72.s1$pop.order),]

BV234.72.s1.sort$pop <- factor(BV234.72.s1.sort$pop, levels=BV234.72.s1.sort$pop)   ##sort pop.nr. Numbers from south to North

qplot(pop, F_MISS, data=BV234.72.s1.sort, geom=c("boxplot", "jitter"))
```

![alt_txt][missing.s1]
[missing.s1]:https://cloud.githubusercontent.com/assets/12142475/20667013/715f51e2-b567-11e6-9383-d57efe9f75c3.png


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


###Step2: HWE

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

And remove the 425 loci using plink
```
plink --file SEsubset.Nov.plink --exclude HWE.loci.remove.names --recode --recodeA --out SEsubset.Nov.s3
```




