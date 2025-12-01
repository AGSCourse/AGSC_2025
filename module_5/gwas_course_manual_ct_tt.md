# GWAS Course Manual
Developed by Chris Kintu and Tsaone Tamuhla

# Table of Contents
- [What is GWAS?](#gwasdef)
  - [Terminologies](#terms)
  

## What is GWAS? <a name="gwasdef"></a>
A Genome-Wide Association Study (GWAS) is a study design that scans the whole genome
to find common genetic variants that influence traits or diseases.

**Goals**:
- Identify genetic risk factors for common, complex diseases (e.g. type 2 diabetes,
schizophrenia, heart disease) and for quantitative traits (e.g. cholesterol levels).
- Understand biology of disease and find potential therapeutic targets.
- Improve risk prediction and enable personalized or precision medicine, including
pharmacogenetics (e.g. warfarin dosing).

**Key early successes**:
- Complement Factor H variants in age-related macular degeneration (AMD).
- Variants in genes affecting warfarin dose and other drug response traits.

## Terminologies <a name="terms"></a>

| Term | Definition |
| :-: | :-: |
| Confounding | A type of bias in statistical analysis causing spurious or distorted findings caused by a correlation between an extraneous variable (the confounding variable) and both the dependent/exposure variable (e.g. the genotype at a given locus) and the independent/outcome variable (e.g. case- control status). |
| Failure rate | The proportion of missing genotypes. Genotypes are classified as missing if the genotype-calling algorithm cannot infer the genotype with sufficient confidence. Can be calculated across each individual and/or SNP. |
| False-negative | Occurs when a true disease-associated variant is not associated with disease in a given study. |
| False-positive | Occurs when a variant not in truth associated with disease status is significantly associated with disease in a given study. |
| Genotype calling algorithm | A statistical algorithm that, per marker and per individual, converts intensity data from two allelic probes into a single genotype for analysis. |
| Genotype Call Rate | The proportion of genotypes per marker with non-missing data. |
| HapMap | An international project to create a haplotype map of the human genome. The publicly available data consists of ~3.2 million SNPs genotyped across four different samples sets of 60-90 individuals of African, Asian or European Ancestry (stage II). HapMap stage III consists of ~1.5 million SNP genotypes from a greater number of individuals and populations. |
| Hardy-Weinberg equilibrium | Given a minor allele frequency of q the probabilities of the three possible genotypes (aa, Aa, AA) at a biallelic locus which is in Hardy-Weinberg equilibrium are ((1−q)2, 2q(1−q), q2). In a large, randomly mating, homogenous population these probabilities should be stable from generation to generation. |
| Heterozygosity rate | The proportion of heterozygous genotypes for a given individual. |
| Linkage Disequilibrium | Non-random association of alleles at two or more loci. |
| Pair-wise identity by state | The proportion of loci where a given pair of individuals share the same alleles. Given by (IBS2 + 0.5Å~IBS1) / (N SNP pairs) where IBS2 and IBS1 are the number of loci where the two individuals have 2 alleles and 1 allele in common, respectively and N SNP pairs is the number of common, non-missing, SNPs. |
| Population substructure | The presence of distinct groups of individuals with subtle differences in allele frequency such that genetic data can be used to cluster these individuals into separate groups. |
| Principal Component Analysis | A mathematical procedure for calculating a number of orthogonal latent variables that summarize a data matrix containing many potentially correlated variables. |
