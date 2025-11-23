# Foundations of Genomics and Genetic Medicine 

## Table of contents
1. [Module Developers and Assistants](#devel)
2. [Learning Outcomes](#lo)
3. [Background and Summary/Objectives](#sum-obj)
4. [Overview of key resources needed](#overview)
5. [Activities and Exercises](#activites)
6. [Other resources](#other-res)
7. [References and Acknowledgements](#ackn)


## Module Developers and Assistants <a name="devel"></a>
- Cesar Fortes-Lima
- Shaun Aron

## Learning Outcomes <a name="lo"></a>

- 2.1 Introduction to the course. 
- 2.2 Foundations of Genetics and Genomics.
- 2.3 Human Genetic Variation.
- 2.4 Patterns of Inheritance.
- 2.5 Principles of Evolutionary and Population Genomics.
- 2.6 Databases.
- 2.7 Population Genomics in Africa.

## Background and Summary/Objectives <a name="sum-obj"></a>
- ...
- ...


## Overview of key resources needed <a name="overview"></a>
- ...
- ...


## Activities and Exercises <a name="activites"></a>
- To obtain an example of a genomic dataset, participants can download and analyse the mitochondrial DNA (mtDNA) and the Y chromosome data presented in the 1000 Genomes Project Phase 3 as follows:

#### Download the data for the mtDNA.
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
```

#### Check the size of the file (202K).
```
ls -lh ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
```

#### Extract individuals with African ancestry using a list of ID samples. This script also extracts only biallelic SNPs and removes duplicated SNPs in the raw dataset.
```
bcftools view -m2 -M2 -v snps -S list_of_AFRsamples_mtDNA.txt ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o 1KGP_African_individuals_mtDNA_BD.vcf.gz
```

#### Check the mtDNA haplogroups of the individuals using Haplogrep.
Upload the generated vcf file and classify

https://haplogrep.i-med.ac.at 


This list has 660 African and African-descendant individuals from 1KGP (GWD, MSL, YRI, ESN, LWK, ACB, and ASW). It takes less than 2 minutes to identify the mtDNA haplogroup using Haplogrep.


#### Download the data for the Y chromosome.
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz
```
#### Check the size of the file (5.4M).
```
ls -lh ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz
```

#### Extract individuals with African ancestry using a list of ID samples.
```
bcftools view -m2 -M2 -v snps -S list_of_AFRsamples_chrY.txt ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o 1KGP_African_individuals_chrY_BD.vcf.gz
```

#### Identify the Y-haplogroup using Y-LineageTracker:
```
LineageTracker classify -b 37 --filter auto --mut-info --vcf 1KGP_African_individuals_chrY_BD.vcf.gz --out 1KGP_chrY
```

More details about this program:
Chen et al. Y-LineageTracker: a high-throughput analysis framework for Y-chromosomal next-generation sequencing data. BMC Bioinformatics 22, 114 (2021). https://doi.org/10.1186/s12859-021-04057-z 


Other genomic data is also available here:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502 




## Other resources <a name="other-res"></a>
- ...
- ...


## References and Acknowledgements <a name="ackn"></a>
- ...
- ...


---

_Under development_
