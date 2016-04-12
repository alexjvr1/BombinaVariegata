#Analysis of BV2


Filter step 1

--vcf BV2N4Clust93.vcf --maf 0.25 --recode --recode-INFO-all --out 1.BV

Eighth Header entry should be INFO: INFO    
After filtering, kept 24 out of 24 Individuals
Outputting VCF file...
After filtering, kept 242490 out of a possible 415468 Sites
Run Time = 13.00 seconds


vcftools --vcf 1.BV.recode.vcf --minDP 3 --recode --recode-INFO-all --out 2.BV

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf 1.BV.recode.vcf
	--recode-INFO-all
	--minDP 3
	--out 2.BV
	--recode

After filtering, kept 24 out of 24 Individuals
Outputting VCF file...
After filtering, kept 242490 out of a possible 242490 Sites
Run Time = 9.00 seconds


vcftools --vcf 2.BV.recode.vcf --missing-indv


vcftools --vcf 2.BV.recode.vcf --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out 3.BV

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf 2.BV.recode.vcf
	--recode-INFO-all
	--maf 0.05
	--max-missing 0.8
	--out 3.BV
	--recode

After filtering, kept 24 out of 24 Individuals
Outputting VCF file...
After filtering, kept 1476 out of a possible 242490 Sites
Run Time = 3.00 seconds
[alexjvr@gdcsrv1 BV.93Filter]$ vcftools --vcf 3.BV.recode.vcf --missing-indv

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf 3.BV.recode.vcf
	--missing-indv

After filtering, kept 24 out of 24 Individuals
Outputting Individual Missingness
After filtering, kept 1476 out of a possible 1476 Sites
Run Time = 0.00 seconds
[alexjvr@gdcsrv1 BV.93Filter]$ mawk '!/IN/' out.imiss | cut -f5 > totalmissing                                                   
[alexjvr@gdcsrv1 BV.93Filter]$ gnuplot << \EOF                                                                                   
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin( $1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

                                         Histogram of % missing data per individual
  Number of Occurrences
      5 ++-----*****------+-----------------+------------------+-----------------+-----------------+----------------++
        +      *   *      +                 +                  'totalmissing' using (bin( $1,binwidth)):(1.0) ****** +
        |      *   *                                                                                                 |
    4.5 ++     *   *                                                                                                ++
        |      *   *                                                                                                 |
        |      *   *                                                                                                 |
      4 ++     *   *                                                                                                ++
        |      *   *                                                                                                 |
        |      *   *                                                                                                 |
    3.5 ++     *   *                                                                                                ++
        |      *   *                                                                                                 |
        |      *   *                                                                                                 |
      3 ++     *   *             *******                                                                            ++
        |      *   *             *     *                                                                             |
        |      *   *             *     *                                                                             |
    2.5 ++     *   *             *     *                                                                            ++
        |      *   *             *     *                                                                             |
        |      *   *             *     *                                                                             |
      2 ++  ****   *****         *     *                                                                            ++
        |   *  *   *   *         *     *                                                                             |
        |   *  *   *   *         *     *                                                                             |
    1.5 ++  *  *   *   *         *     *                                                                            ++
        |   *  *   *   *         *     *                                                                             |
        +   *  *   *   *  +      *     *    +                  +                 +                 +                 +
      1 ++--***************************************************************************************************-----++
        0                0.05              0.1                0.15              0.2               0.25              0.3
                                                      % of missing data




vcftools --vcf 3.BV.recode.vcf --out BV2.plink --plink 

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf 3.BV.recode.vcf
	--out BV2.plink
	--plink

After filtering, kept 24 out of 24 Individuals
Writing PLINK PED and MAP files ... 
	PLINK: Only outputting biallelic loci.
Done.
After filtering, kept 1476 out of a possible 1476 Sites
Run Time = 0.00 seconds


#correct the plink file

plink --file BV3.75.plink --recodeA


##29 Feb

Trying with 75% missing data
vcftools --vcf 2.BV.recode.vcf --max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out 3.BV

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf 2.BV.recode.vcf
	--recode-INFO-all
	--maf 0.05
	--max-missing 0.75
	--out 3.BV
	--recode

After filtering, kept 24 out of 24 Individuals
Outputting VCF file...
After filtering, kept 2683 out of a possible 242490 Sites
Run Time = 4.00 seconds


And convert to plink...

plink --file BV3.75.plink --recodeA

#copy everything to the BV2 folder /BV2/75

scp -r alexjvr@gdcsrv1.ethz.ch:/gdc_home4/alexjvr/Bombina/BV2/pyrad93/outfiles.1/BV.93Filter/BV3.75* .

#and the plink.raw file


##IN R

library("ade4")
library("adegenet")
library("pegas")




BV.75 <- read.PLINK("plink.raw") ##reads a plink file in as genlight object

BV.75 #check that it is a genlight object

indNames(BV.75)

indiv.names <- read.table("BV2.75.indiv.csv", header=T, quote="\"") #read in the indiv names file
indiv.names.factors <- as.factor(indiv.names$indiv) #and convert to a factor
summary(indiv.names.factors)

BV.75

pop(BV.75) <- (indiv.names.factors) #assign population names from a text file
pop(BV.75)  ##and check that they are correct


##Pop structure
pca1 <- glPca(BV.75) ##This displays a barplot of the eigenvalues and asks the user for a number of retained principal components

20

head(pca1)

s.class(pca1$scores, pop(BV.75), col=colors()[c(131,131,131,131,131,131,131,131,160,160,160,160,160,150,150,150,150,134,134,134,134,134,134,134)])
add.scatter.eig(pca1$eig,2,1,2)
abline(h=0,v=0,col="grey")

myCol <- colorplot(pca1$scores,pca1$scores,transp=T,cex=4)
abline(h=0,v=0,col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)


 

PCA of BV2 (25% missingness, 2675 SNPs)

Individuals with missing data: 

 


 


##This is a good tutorial for dapc & pca
http://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

BV2.75.genind

pop.names <- read.table("BV2.75.pop.csv", header=T, quote="\"") #read in the indiv names file
pop.names.factors <- as.factor(pop.names$pop) #and convert to a factor
summary(pop.names.factors)


pop(BV2.75.genind) <- (pop.names.factors) #assign population names from a text file
pop(BV2.75.genind)  ##and check that they are correct

dapc.BV2.75 <- dapc(BV2.75.genind, var.contrib = TRUE, scale = FALSE, n.pca = NULL, n.da = 4 - 1)

scatter(dapc.BV2.75, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

With n.pca = 8:
 


##It's important to choose npca and ndapc carefully (i.e. using BIC). nPCA greatly affects the outcome of the test!

##Now we can check whether any specific loci contribute the most to the differentiation: 

set.seed(4)
contrib <- loadingplot(dapc.BV2.75$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)

 


##And now check the DAPC



BV2.75.genind2 <- na.replace(BV2.75.genind, "mean")

set.seed(999)
xval1 <- xvalDapc(BV2.75.genind2$tab, pop.names.factors, n.pca.max = 300, n.da = NULL, training.set = 0.9, result = c("groupMean", "overall"),center=TRUE, scale = F, n.pca = NULL, n.rep = 30, xval.plot = T)

 
This is not very convincing at all! So is there not enough differentiation? 8 PC's looks like the best option, but there's still a large CI.. 


