# Genome wide Association Studies - Practical Session


## Quality control in genome-wide association studies

This practical is based on “Data quality control in genetic case-control association studies” 
[Anderson et al. 2010, Nature Protocols 5: 1564-73](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025522/pdf/ukmss-33586.pdf).

- Data for this practical is provided in the course data directory [here](../data) as `raw-GWA-data.tgz`.
- You may click [here](https://github.com/AGSCourse/AGSC_2025/raw/refs/heads/main/data/raw-GWA-data.tgz) to download it directly.

### Instructions
This practical will introduce you to performing quality control for Genome Wide Association Study (GWAS) 
datasets using a program called PLINK. You might finish quicker if you already have experience working on 
PLINK and sufficient programming experience. Please read the instruction before starting.


1.	This practical session assumes basic programming (eg bash scripting) and computing experience.

2.	You will need to work on a Unix/Linux environment (shell prompt)

3.	You will need to install PLINK if not already installed on your computer. See below on how to install PLINK.

4.	You should also have R installed

5.	Ideally, you should run the analysis on High-performance computing (HPC) environment.
   It would take longer if you run the analysis on your stand-alone computer/laptop


## Installing and running PLINK (in your working directory)

PLINK Website: https://www.cog-genomics.org/plink/2.0/
```
wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20210203.zip
Unzip plink2_linux_x86_64_20210203.zip
```

## PLINK general usage
Website: https://www.cog-genomics.org/plink/2.0/general_usage

The old PLINK manual is quite useful. It is here: https://zzz.bwh.harvard.edu/plink/dist/plink-doc-1.07.pdf 

**Remember that there are many versions of PLINK. The current version is PLINK2 which can be found on this link**: https://www.cog-genomics.org/plink/2.0/

## Quality control in genome-wide association studies
**Instruction**: 
In this practical, we will go through the steps in performing quality control (QC) of genotype data from a simulated 
genome-wide association study of 1000 cases and 1000 controls, typed for 317,503 autosomal and X chromosome SNPs.

We will begin by performing sample QC, including calculation of call rates, heterozygosity, and sex discordance.
We will then perform SNP QC, including calculation of call rates and deviation from Hardy-Weinberg equilibrium.
This practical is based on “Data quality control in genetic case-control association studies” (Anderson et al. 2010, Nature Protocols 5: 1564-73).

1.	Start by downloading and unpacking the data from the Amathuba platform into your working folder on your personal computer or an HPC. 

Unpack the data 
```
tar –xvzf raw-GWA-data.tgz
```

**Genotype Data**:
- `raw-GWA-data.ped`
- `raw-GWA-data.map`

You may also have noticed that there are some scripts (in `R` and `PERL`) we would use later in the downloaded data.

The first step is to convert your `raw-GWA-data.ped` and `raw-GWA-data.map` to binary file
```
./plink --file raw-GWA-data -–make-bed --out raw-GWA-data
```

**Note**: code above assumes you have downloaded and unpacked plink in the same directory as your uncompressed data

> **Questions**
> Q1. How many individuals do you have in the data?
> 
> Q2. How many SNPs?
> 
> Q3. How many missing phenotypes?
> 
> Q4. How many cases?
> 
> Q5. How many controls?


## Part A: Sample QC

**Step 1**: Identification of individuals with discordant sex information

At the shell prompt, type:
```
./plink --bfile raw-GWA-data --check-sex --out raw-GWA-data
```
File generated `raw-GWA-data.sexcheck`

This command will calculate the mean homozygosity rate across X-chromosome markers for each individual in the study. 
We can produce a list of individuals with discordant sex data by typing:

```
grep PROBLEM raw-GWA-data.sexcheck > raw-GWA-data.sexprobs awk '{if ($5=="PROBLEM")print $0}' raw-GWA-data.sexcheck | head
```

> **Question**:
>> Q6. How many discordant sex do we have? 

Open the file “raw-GWA-data.sexprobs” to obtain the family IDs (column 1) and individual ID (column 2) for these individuals. 
Column 3 denotes ascertained sex and column 4 denotes sex according to genotype datw, the genotype data are inconclusive 
regarding the sex of an individual and these are marked in column 4 by “0”.

In general, any discordance in sex should be reported to study co-ordinators to double check records for errors. 
In situations in which discrepancy cannot be resolved, add the family ID (FID) and individual ID (IID) of the samples 
to a file named `fail-sexcheck-qc.txt` (one individual per line, tab delimited). This file can be used to exclude these samples 
from downstream analyses.

**Step 2**: Identification of individuals with elevated missing data rates or outlying heterozygosity rate

At the shell prompt, type:
```
./plink --bfile raw-GWA-data --missing --out raw-GWA-data
```

This command will create the files “raw-GWA-data.imiss” and `raw-GWA-data.lmiss`. 
The fourth column in the file `raw-GWA-data.imiss` (N_MISS) denotes the number of missing SNPs and the sixth column 
(F_MISS) denotes the proportion of missing SNPs per individual.

At the shell prompt type:
```
./plink --bfile raw-GWA-data --het --out raw-GWA-data
 ```

This command will create the file `raw-GWA-data.het`, in which the third column denotes the observed number of 
homozygous genotypes `O(Hom)` and the fifth column denotes the number of non-missing genotypes `N(NM)` per individual.

You can calculate the observed heterozygosity rate per individual using the formula:

`Het = (N(NM) − O(Hom))/N(NM)`

Create a graph in which the observed heterozygosity rate per individual is plotted on the x axis and the proportion of 
missing SNPs per individuals is plotted on the y axis. This can be carried out using standard software such as Excel or R.

**Important**

First install the plotting tool `geneplotter` in `R`.

At the shell, type:
```
R
```

At the R console, type:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("geneplotter")
```

You will be prompted with:

`Update all/some/none? [a/s/n]:`

Type:
`n`

Then exit R:
```
q()
```

Prompt:

`Save workspace image? [y/n/c]:`

Type:
`y`


A script for calculating the heterozygosity rate and producing the graph using R is supplied: “imiss-vs-het.Rscript”. 
At the shell prompt, type:
```
R CMD BATCH imiss-vs-het.Rscript
```

This will create the graph `raw-GWA-data.imiss-vs-het.pdf`.

> **Question**
> 
> Q7. Examine the plot to decide reasonable thresholds at which to exclude individuals based on elevated
> missing genotype rates or extreme heterozygosity.

The red lines indicate suggested thresholds: missing data rate of 3% and heterozygosity rate
`±3 SD` from the mean. This file can easily be adapted to highlight different thresholds. Add the FID and IID 
of the samples failing this QC to a file named `fail-imisshet-qc.txt` (one individual per line, tab delimited). 
This file can be used to exclude these samples from downstream analyses.

**Step 3**: Identification of duplicated or related individuals

To minimize computational complexity, create an “independent” set of SNPs to generate the identity by descent (IBS) 
matrix. This can be done by “pruning” the data so that no pair of SNPs (within a given genomic interval) has an r2 
value greater than a given threshold, typically chosen to be 0.2. Typically, we would also exclude regions of strong 
LD, such as the MHC, which are listed in the file `high-LD-regions.txt`.

Please note that the data for this practical have been simulated without LD between SNPs. We have created a **pruned** 
list of SNPs to use called “raw-GWA-data.prune.in” which you should use for this practical. However, in general, 
you would use the following command to generate the pruned SNP list:
```
./plink --bfile raw-GWA-data --exclude high-LD-regions.txt --indep-pairwise 50 5 0.2 --out raw-GWA-data
```

To generate IBS between each pair of individuals, type the following command at the shell prompt:
```
./plink --bfile raw-GWA-data --extract raw-GWA-data.prune.in --genome --out raw-GWA-data
```

Please note that this step is computationally demanding and will take some time to run.

We have created a script file to identify all pairs of individuals with IBS > 0.185 which outputs the ID of the 
individual from the pair with lowest call rate (from the previously created file `raw-GWA-data.imiss`). 
This script file can easily be adapted to other data sets and thresholds for **related** individuals.

To run the script file, type the following command at the shell prompt:
```
perl run-IBD-QC.pl raw-GWA-data
```

The command creates a file called “fail-IBD-QC.txt” which can be used to exclude these samples from downstream analyses.

**Step 4**: Removal of all individuals failing sample QC

At the shell prompt, type the following command to concatenate all files listing individuals who have failed previous QC steps:
```
cat fail-* | sort -k1 | uniq > fail-qc-inds.txt
```

The file `fail-qc-inds.txt` should now contain a list of unique individuals failing the previous QC steps. 
To remove them from the data set, type the following command at the shell prompt:
```
./plink --bfile raw-GWA-data --remove fail-qc-inds.txt --make-bed --out clean-inds-GWA-data
```

The binary ped file `clean-inds-GWA-data` can be used for subsequent SNP QC analyses.


## Part B: SNP Quality Control

**Step 1**: Identification of all SNPs with an excessive missing data rate

To calculate the missing genotype rate for each SNP, type the following command at the shell prompt:
```
./plink --bfile clean-inds-GWA-data --missing --out clean-inds-GWA-data
```

Output file

The third column in the file `clean-inds-GWA-data.lmiss` (`N_MISS`) denotes the number of missing genotypes and 
the fifth column (`F_MISS`) denotes the proportion of missing genotypes per SNP.

A script for plotting a histogram of the missing genotype rate using R is supplied: `lmiss-hist.Rscript`. At the shell prompt, type:
```
R CMD BATCH lmiss-hist.Rscript
```

This will create the graph `clean-inds-GWA-data.lmiss.pdf`.

> **Question:**
>
> Q8. Examine the plot to decide a reasonable threshold at which to exclude SNPs based on elevated missing data.


The dotted line indicates a suggested threshold of a missing data rate of 5%. This file can easily be adapted to highlight different thresholds.

To remove SNPs with call rate less than 5%, simply add the `--geno 0.05` option to the PLINK command line. 
We will do this below when creating our final cleaned data set.

**Step 2**: Test SNPs for different genotype call rates between cases and controls

To test for differential call rates between cases and controls for each SNP, type the following command at the shell prompt:
```
./plink --bfile clean-inds-GWA-data --test-missing --out clean-inds-GWA-data
```

Output file

The output of this test can be found in the file “clean-inds-GWA-data.missing”. We have created a script to highlight all SNPs with
significant differences in case and control call rates (p<10-5) from this output file. This script can be easily 
adapted to other data sets and thresholds for differential call rates.

To run the script file, type the following command at the shell prompt:
```
perl run-diffmiss-qc.pl clean-inds-GWA-data
```

Output file

The command creates a file called `fail-diffmiss-qc.txt`, which can be used to exclude these SNPs from downstream association analyses.

**Step 3**: Removal of all SNPs failing QC

To remove low-quality SNPs, type the following command at the shell prompt:
```
./plink -- bfile clean-inds-GWA-data --exclude fail-diffmiss-qc.txt --geno 0.05 --hwe 0.00001 --make-bed --out clean-GWA-data
```

In addition to removing SNPs identified with differential call rates between cases and controls, this command removes 
SNPs with call rate less than 5% with `--geno` option and deviation from HWE (p<10-5) with the `--hwe` option. One additional 
option is `--maf` which excludes all SNPs with minor allele frequency less than a specified threshold.

This command will produce cleaned binary ped files for downstream association analyses:
- `clean-GWA-data.bed`
- `clean-GWA-data.bim`
- `clean-GWA-data.fam`


In the following **Practical Two** section, Use the quality controlled data (`clean-GWA-data.*`) generated above

## Practical Two

Basic analysis of genome-wide association studies

In this practical, we will perform basic analysis of genotype data from a simulated genome-wide association 
study of 2000 cases and 2000 controls, typed for 306,102 autosomal SNPs. You can assume that the data have 
already assessed for quality control, with low quality samples and SNPs removed. We will perform logistic 
regression analysis under a range of disease models and allow for adjustment for covariates.

This practical is based on “Basic statistical analysis in genetic case-control studies” (**Clarke et al. 2011, Nature Protocols 6: 121-33**).

In order to run this practical, you will require the following resources:

1.	Computing workstation with Unix or Linux operating system;
2.	PLINK software (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml);
3.	Haploview software (http://www.broadinstitute.org/haploview/haploview);
4.	Genome-wide association binary ped files (which can be downloaded from http://www.well.ox.ac.uk/dtc/).

Before starting the practical, you will need to unpack the genome-wide association binary ped files and additional 
analysis script files. You can do this by typing the following command at the shell prompt in the directory 
where you downloaded the file `raw-GWA-data.tgz`:
```
tar -xvzf gwa-data.gz
```

Test for association with disease status under an additive model

To test for association of SNPs with disease under an additive model (multiplicative on the odds scale), type the following command at the shell prompt:
```
./plink --bfile gwa --logistic --ci 0.95 --out additive.analysis
```

For each SNP, the output file “additive.analysis.assoc.logistic” contains the following information: ID, chromosome and position, 
minor allele, odds ratio for the minor allele and the corresponding 95% confidence interval, and the p-value for association.

To produce a Manhattan plot to summarise the output of this analysis:
- begin by starting the Haploview software. 
- Click the “PLINK Format” tab on the left hand side.
- Use the “Browse” buttons to select `additive.analysis.assoc.logistic` as the **Results File** and `gwa.bim` as the **Map File**.
- Then click **OK** to read in the data, which may take some time. Once the data is read in, it will appear in a new window.
- Click the `Plot` button and in the pop-up window select `P` for the **Y-Axis**, and `-log10` for the **Scale**.
- In the `Significant` box, you can optionally add the threshold for genome-wide significance, which on a `-log10` scale is `7.30`.
- Finally, click **OK** to produce the plot.

> **Question:**
>
> Do you find any evidence of association at genome-wide significance?

## Test for association with disease under a genotypic model

To test for association of SNPs with disease under a general genotypic model (two degree of freedom test), type the following command at the shell prompt:
```
./plink --bfile gwa --logistic genotypic --ci 0.95 --out genotypic.analysis
```

For each SNP, the output file “genotypic.analysis.assoc.logistic” contains three rows of results, one for the 
additive term from the general model, one for the dominance term, and one for the genotypic test. For the additive 
and dominance terms, the odds ratio and corresponding 95% confidence interval for the minor allele are presented. 
The p-value for association is also presented for the additive and dominance terms, as well as the genotypic test.

You can produce a Manhattan plot in the same way as before. However, after reading in the data, you will need to filter the 
data to focus on the association test results you want to present. For example, to present the general genotypic association 
test results: 
- select `TEST` from the `Filter` tab and then `=` and `GENO_2DF`, before clicking the `Filter` button.

> **Question:**
>
> How do the results compare between a genotypic and additive test of association? Is there any evidence of a deviation 
> from additivity (by considering the dominance test)?

## Test for association with disease allowing for covariates

To test for association under an additive model, allowing for non-genetic risk factors, it is necessary to provide 
PLINK with a covariate file, here “gwa.covar”, which provides one row per individual in the study, with covariates 
arranged in columns. To adjust for “AGE” from the covariate file, and also for sex (which is present in the fam file, 
so is not needed in the covariate file), type the following command at the shell prompt: 
```
./plink --bfile gwa --logistic sex hide-covar --ci 0.95 --covar gwa.covar --covar-name AGE --out additive.AGE.SEX.analysis
```

The option `hide-covar` suppresses printing of the parameter estimates for the covariate terms from the logistic regression 
model. The output file `additive.AGE.SEX.analysis` contains the same information for each SNP as before, but this 
time with the odds ratio for the minor allele, and p-value for association, adjusted for age and sex.

You can produce a Manhattan plot in the same way as before. Have the results changed after adjustment
