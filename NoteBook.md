#*Bombina variegata* microsatellite genotyping

The aim of this project is to test and genotype 12 microsatellites in *Bombina variegata*. If necessary, additional markers will be developed and genotyped to supplement the dataset. 

1. Extract DNA from xx mouth swabs
2. Test developed *B. variegata* primers
3. Genotype all samples with optimised marker panel
4. Score msats
5. Order and develop additional primers if necessary
6. Genotype samples with additional primers



##Background info

- 12 *B. variegata* primers are already available
(papers of Hauswaldt et
al., 2007; Stuckas & Tiedemann, 2006).

Primer name |sequence | GenBank| Repeat| Ta| Allele size in *B.variegata*|Ref
:--:|:--:|:--:|:--:|:--:|:--:|:--:
BobomF2.for|5’AGC AGA GAT GAG AGG ACA GTG|AJ877247 |(GA)20|60|480, 456| Hauswaldt *et al*. 2007
BobomF2.rev|5’TCA GGG GTA GCA GAT TTT CAG||||
BobomF22.for |5′AGG CAA AGG ATT CTG AGA ATG|AJ877250| (GA)30|56| 132| Hauswaldt *et al.* 2007
BobomF22.rev| 5′CCT TCA AAG TCG AAA AAT ATT||||
BobomB13.for| 5′ATA TTT CTT GCT ATG TTG ATG|AJ877251| (GA)22|46 |126, 134|Hauswaldt *et al.* 2007
BobomB13.rev| 5′AAT TGT TTA ACT TAT TTT ATA||||
BobomB14.for| 5′ACT AAC CTG CCA CAT AAC TTG|AJ877254| (TC)4T(TC)7T|48| 184, 186
BobomB14.rev| 5′CTG GGT TTT TTA ATT GGA AGG| (TC)9GC(TC)7||||

- Some of these primers were problematic when applied to populations from Alsace region (2 problematic, 2 fixed)
- Similar for populations around Lac Léman (1 problematic, 3 fixed)

> So I can suggest to add the last 2 ones (B14 and 12F) to our primers list
and thus test with 12 microsatellites markers; but it's possible that one or
two will not work and 1-2 could be not so informative.

- NCBI sequences for potential msat primers: 

> In addition, some microsatellite sequences of BOBO are available on
GenBank. It would be possible to develop and test them, but it would need
more time (design the primers, test the amplification, see if we can improve
the PCR condition, and see if they are variable). All these steps are the
most time consuming aspect of developing microsat.

>>> there are more developped for BOVA or/and BOBO, but I should ask
Jean-Pierre Vacher (the first that conducted the genetic analyses on
Bombina) why he selected only those ones.


##ddRAD

Since the aim of the study is to infer movement between populations, it seems like using many markers would give the highest probability of finding these signatures in a population with low expected genetic divergence. 

We will try a ddRAD approach on the *B.variegata*

**B.variegata* has a very large genome: 9-10pg ~10Gb.

This is >2x as large as Rana temporaria. Since I am already getting rather low coverage for my samples, using a different enzyme combination is probably wise. 

- For Rana temporaria, I am currently using a 6-cutter (EcoRI) and a 4-cutter (MseI). 
- GeCKo protocol uses a 6-cutter (PstI) and a 4-cutter (MspI)
- Alan Brelsford has switched to an 8-cutter (SbfI) with a 4-cutter (MseI)


How many loci are needed? 

- I'm generating 100-500 000 loci per individual for *R.temp*
- We don't need as many for the population genetics study
- It depends on the expected genetic diversity in the population (large pop sizes & meta-population structure -> unlikely to be highly inbred)
- 

2 samples each were chosen from 2 locations. 

Aim: 

1. Find enzyme combination that will produce ~50000-100000 loci. 
2. Assess genetic diversity across all individuals, and per pop
3. Determine genetic distance between populations (per locus Fst, PCA, Structure)
4. 


#1. Which enzyme combination is best?

The worry is that we sequence too many loci & lose them due to low coverage across loci. This is not such a problem with many individuals per population (probabilitstic models for calling bases), but not good enough for the test. 

We are bound to EcoRI, because of our adapters (from the GDC). 

1. EcoRI + MseI (6 + 4-cutter: same as my protocol)
2. EcoRI + MspI (6 + 4-cutter: probably the same as in 1)
3. EcoRI + PstI (6 + 6-cutter: likely the best option)
