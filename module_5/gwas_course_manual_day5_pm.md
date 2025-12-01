# GWAS Course Manual
Developed by Chris Kintu and Tsaone Tamuhla

# Table of Contents <a name="toc"></a>
1. [What is GWAS?](#gwasdef)
   -   [Terminologies](#terms)
2. [Basic genetic concepts behind GWAS](#concepts)
   -   [Single Nucleotide Polymorphisms (SNPs)](#snps)
3. [Why GWAS? From rare diseases to common diseases](#whygwas)
   -   [Linkage and its limitations](#linkage)
   -   [Common disease / common variant (CD/CV) hypothesis](#cdcv)
4. [Capturing common variation](#ccv)
   -   [The HapMap Project](#hapmap)
   -   [Linkage Disequilibrium (LD)](#ld)
   -   [Tag SNPs and indirect association](#tags)
5. [Genotyping technologies](#gtechs)
6. [Study design and phenotyping](#sdesign)
   -   [Phenotype types](#pheno)
   -   [Standardized phenotype criteria](#pheno-criteria)
   -   [Phenotypes from electronic medical records (EMRs)](#pheno-electro)
7. [Association analysis](#assoc)
   -   [Single-locus tests](#single-locus-test)
   -   [Genetic models (how we code genotype)](#models)
   -   [Covariate adjustment & population stratification](#covar)
   -   [Multiple testing](#multi-test)
   -   [Multi-locus (interaction) analysis](#multi-locus)
8. [Replication, meta-analysis and imputation](#repli-meta)
   -   [Meta-analysis](#meta)
   -   [Genotype Imputation](#impute)
9. [Looking beyond GWAS](#beyond-gwas)
10. [GWAS Quality Control](#qc)
    -   [Overview](#overview)
    -   [Learning Outcomes](#learn-obj)
    -   [Why QC Matters](#why-qc-matters)
    -   [Overview of the QC Workflow](#qc-work)
    -   [Per-Individual QC](#sample-qc)
    -   [Per-Marker QC](#snp-qc)
11. [References](#refs)
   
  

## 1. What is GWAS? <a name="gwasdef"></a>
A **Genome-Wide Association Study (GWAS)** is a study design that scans the whole genome
to find common genetic variants that influence traits or diseases.

**Goals**:
- Identify **genetic risk factors** for common, complex diseases (e.g. type 2 diabetes,
schizophrenia, heart disease) and for quantitative traits (e.g. cholesterol levels).
- Understand **biology of disease** and find potential therapeutic targets.
- Improve **risk prediction** and enable **personalized or precision medicine**, including
pharmacogenetics (e.g. warfarin dosing).

**Key early successes**:
- Complement Factor H variants in **age-related macular degeneration (AMD)**.
- Variants in genes affecting **warfarin dose** and other drug response traits.

## Terminologies <a name="terms"></a>

| Term | Definition |
| :----: | :----: |
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

[Back to top](#toc)

## 2. Basic genetic concepts behind GWAS <a name="concepts"></a>
### 2.1 Single Nucleotide Polymorphisms (SNPs) <a name="snps"></a>
- A **SNP** is a single base-pair change in the DNA sequence (e.g. A→G).
- Most GWAS use SNPs as **markers** that tag regions of the genome.
- Most SNPs are neutral, but some can:
   - Change amino acids in proteins.
   - Affect mRNA stability or splicing.
   - Alter transcription factor binding and gene regulation.
 
**Alleles and frequency**
- Most SNPs have two alleles in the population (e.g. A and G).
- We describe them using **minor allele frequency (MAF)**:
   - If G is present in 40% of chromosomes, MAF = 0.40.
- In GWAS we typically focus on common variants (MAF ≥ ~1–5%).

**SNPs vs mutations**
- In practice:
   - **SNP**: a common single-base variant.
   - **Mutation**: a very rare variant that often has a large functional effect (classic in Mendelian diseases like cystic fibrosis).
 

## 3. Why GWAS? From rare diseases to common diseases <a name="whygwas"></a>
### 3.1	Linkage and its limitations <a name="linkage"></a>
- For **rare, highly penetrant diseases** (e.g. cystic fibrosis, Huntington disease), family-based **linkage analysis** works well.
   - Track inheritance of markers and disease within families.
   - Identify genomic regions where markers and disease co-segregate.
- For **common complex diseases** (heart disease, cancer, diabetes), linkage has largely failed:
   -	Genetic effects are smaller.
   -	Many loci contribute → Need a different approach: population-based association.

## 3.2 Common disease / common variant (CD/CV) hypothesis <a name="cdcv"></a>
- Proposes that **common diseases are partly due to common variants**.
- Implications:
   1.	**Effect sizes are small**: A common allele with very large effect would make the disease much more prevalent.
   2.	**Many loci are involved**: Heritability of common disease is spread across many common variants, each explaining a small fraction of variance.
- Consequence: we need **large sample sizes** and **dense marker panels**, and we focus on **population-based association** rather than pedigree linkage.

## 4. Capturing common variation <a name="ccv"></a>
### 4.1 The HapMap Project <a name="hapmap"></a>
To do GWAS efficiently we must know:
- Where common variants are.
- How they differ across populations.
- How variants are correlated with each other.

The **International HapMap Project:**
- Catalogued millions of SNPs in multiple populations (initially European, Yoruba, Han Chinese, Japanese; later more).
- Measured **linkage disequilibrium (LD)**: correlation structure between SNPs along chromosomes.

### 4.2 Linkage Disequilibrium (LD) <a name="ld"></a>
- LD describes how alleles at two loci are **correlated** in a population.
- If two SNPs are in high LD (high r²), knowing one tells you about the other.
- Over many generations, **recombination** breaks up chromosomes; LD blocks become shorter.

**Population differences**:
- African-ancestry populations: shorter LD blocks (older population, more recombination).
- European/Asian-ancestry populations: longer LD blocks (founder events, fewer generations).

**Key measures**:
- **D′**: ranges 0–1, tells you whether recombination between two loci has occurred.
- **r²**: correlation measure, used in GWAS; high r² means two SNPs carry almost the same information.

### 4.3 Tag SNPs and indirect association <a name="tags"></a>
- Because of LD, we do not need to genotype every SNP.
- We choose **tag SNPs** that capture the variation of nearby SNPs (high r²).
- In GWAS, a significant SNP can be:
   - **Directly associated** (it is itself causal), or
   - **Indirectly associated** (it tags a nearby untyped causal variant).
-	Therefore, a GWAS hit marks a region; follow-up is required to find the true causal variant.

## 5. Genotyping technologies <a name="gtechs"></a>
GWAS became feasible thanks to **high-throughput SNP arrays**:
- **Illumina** and **Affymetrix** are the main platforms.
- They use different chemistries (bead-based vs probe-based hybridization) but similar principles: genotype hundreds of thousands to millions of SNPs in parallel.

**Important considerations**:
- Coverage differs by population: Because LD is shorter in Africans, more SNPs are needed to cover the genome.
- Technology is evolving: **next-generation sequencing** is gradually replacing arrays, giving complete sequence rather than a subset of SNPs.

## 6. Study design and phenotyping <a name="sdesign"></a>
Good genotyping is useless without good **phenotypes**.

### 6.1	Phenotype types <a name="pheno"></a>
1. Quantitative traits
   - Continuous measures (e.g. LDL levels, BMI, blood pressure).
   - Generally, **more powerful** for detecting genetic effects.
   - Interpretation is straightforward: change in trait per allele.
2. Categorical traits (typically case–control)
   - Disease status: affected vs unaffected.
   - Essential when no good quantitative proxy exists (e.g. multiple sclerosis).
   - Still very successful if phenotype definitions are robust.

### 6.2 Standardized phenotype criteria <a name="pheno-criteria"></a>
- Use clear, reproducible clinical definitions (e.g. McDonald criteria for multiple sclerosis).
- For multi-center studies, standardized rules reduce **misclassification** and **site effects**.
- Sometimes inter-rater reliability is examined to ensure clinicians apply criteria consistently.

### 6.3 Phenotypes from electronic medical records (EMRs) <a name="pheno-electro"></a>
- Biobanks linked to EMRs provide huge sample sizes.
- Challenges:
   - EMRs are built for healthcare and billing, not research.
   - Need algorithms that combine diagnostic codes, procedure codes and sometimes **natural language processing** of free text.
   - Algorithms are iteratively refined and validated by chart review (evaluate positive predictive value).

## 7. Association analysis <a name="assoc"></a>
Once we have genotypes and phenotypes, we test each SNP for association.

### 7.1 Single-locus tests <a name="single-locus-test"></a>
- Quantitative traits:
   - Usually analysed using **linear regression** or ANOVA.
   - Outcome: trait value; predictor: genotype plus covariates.
- Binary traits (case–control):
   - Either:
      - Chi-square / Fisher’s exact test on 2×3 or 2×2 genotype/allele count tables, or
      - **Logistic regression**, which models the log-odds of disease as a function of genotype and covariates.
   - Logistic regression is preferred because it:
      - Adjusts for covariates.
      - Provides **odds ratios** as effect sizes.
      - Has good diagnostics.

### 7.2 Genetic models (how we code genotype) <a name="models"></a>
For a biallelic SNP with alleles A and a, genotype can be encoded as:
- **Additive model**: 0, 1, 2 copies of A (assumes linear effect).
- **Dominant model (A)**: AA/Aa vs aa.
- **Recessive model (A)**: AA vs Aa/aa.
- **Multiplicative, genotypic, allelic** models, etc.

In practice:
- Most GWAS report additive models only:
   - Reasonable power for additive and many dominant effects.
   - Underpowered for purely recessive effects, but computationally simpler.

### 7.3 Covariate adjustment & population stratification <a name="covar"></a>
We typically adjust for:
- Demographics: age, sex.
- Study design: site, batch.
- Clinical factors: BMI, medication use, etc.
- Population structure:
   - Allele frequencies and disease prevalence differ across ancestries.
   - If not corrected, this **population stratification** can create spurious associations.
   - Approaches:
      - Use ancestry estimation tools (e.g. **STRUCTURE**, **EIGENSTRAT**).
      - Include principal components (PCs) of genotype data as covariates in regression.

### 7.4 Multiple testing <a name="multi-test"></a>
A GWAS tests hundreds of thousands to millions of SNPs → huge multiple-testing burden.

Methods:
- **Bonferroni correction**
   - Adjust significance threshold: α* = 0.05 / number_of_tests.
   - Very conservative because tests are not fully independent.
   - Leads to typical GWAS threshold around 5×10⁻⁸–10⁻⁷.
- **False Discovery Rate (FDR)**
   - Controls expected proportion of false positives among declared hits.
   - Useful when many associations are expected (e.g. eQTL, expression studies).
- **Permutation testing**
   - Empirical approach: repeatedly shuffle phenotypes to break real associations.
   - Build null distribution of test statistics, derive empirical p-values.
   - Computationally expensive but robust.
- **Genome-wide significance**
   - Based on the effective number of independent tests considering LD.
   - For European ancestry, commonly ~5×10⁻⁸ is used.
   - Appropriate only for truly genome-wide scans, not for focused candidate-gene or replication studies.

### 7.5 Multi-locus (interaction) analysis <a name="multi-locus"></a>
- GWAS enables investigation of:
   - **Gene–gene interactions (epistasis)**.
   - **Gene–environment interactions**.
- Challenges:
   - Combinatorial explosion: with 1 million SNPs, pairwise combinations are ~5×10¹¹.
   - Limited power, multiple testing, computational burden.

Strategies:
- Filter SNPs first (e.g. by main-effect p-value, biological pathway membership).
- Use tools that integrate biological knowledge (e.g. pathway databases) to narrow down search space (e.g. Biofilter, INTERSNP).
- Recognize that some interaction models may have **weak marginal effects**, so filtering only on main effects can miss “pure” epistasis.

## 8. Replication, meta-analysis and imputation <a name="repli-meta"></a>
### 8.1 Replication <a name="repli"></a>
A credible GWAS finding should be:
- Tested in an **independent sample**.
- Using **similar phenotype definitions** and **similar population**.
- With adequate power (often larger than the discovery sample, to overcome “winner’s curse”, where initial effect sizes are over-estimated).

Key points:
- The **unit of replication** is really the **region**, not necessarily the exact SNP, because of LD.
- It is acceptable to replicate with a different SNP in high LD, as long as effect direction and magnitude are consistent.

### 8.2 Meta-analysis <a name="meta"></a>
- Enables combining results from multiple GWAS without sharing raw data.
- Requires:
   - Harmonized phenotypes and covariate adjustments.
   - Common reference genome build and allele coding.
   - Careful QC to ensure samples are not duplicated across studies.
- Heterogeneity across studies is assessed using statistics like **Q** and **I²**.
- Meta-analysis has powered very large discoveries (e.g. 95 lipid loci from >100,000 individuals).

### 8.3 Genotype imputation <a name="impute"></a>
Problem: different studies use different SNP arrays → different sets of SNPs.

Solution: **imputation**:
- Use reference panels (HapMap, 1000 Genomes) with dense genotype/sequence data.
- Use LD and haplotype information to infer untyped genotypes in your study.
- Outputs are often **genotype probabilities** or **dosages** rather than hard calls.
- Analysis must account for uncertainty (often via Bayesian or likelihood approaches).

Requirements:
- Reference panel must match the study population’s ancestry reasonably well.
- Allele coding must be aligned between study and reference.

## 9. Looking beyond GWAS <a name="beyond-gwas"></a>
GWAS has:
- Identified thousands of loci for hundreds of traits.
- Shifted genetics from single-gene thinking to **genome-wide**, polygenic views.
- Opened doors to:
   - **Functional follow-up** (e.g. lab experiments, expression studies).
   - **Polygenic risk scores**.
   - **Pharmacogenomics** and precision medicine.

Emerging directions:
- **Whole-genome sequencing**: replacing array-based genotyping, capturing both common and rare variation.
- Integration of **multi-omics** (genomics, transcriptomics, proteomics, metabolomics) and environmental data.
- Linking genetics with high-dimensional phenotypes such as **imaging**, **wearables**, or detailed clinical trajectories.

## 10. GWAS Quality Control <a name="qc"></a>
**MATERIALS & EQUIPMENT**

**Data**
- Genome-wide SNP genotype data is provided on Amathuba or via this link [here](#)
- Associated documentation or analysis scripts are provided on Amathuba or via this link [here](#)

**Hardware**
- Computer workstation running a Unix/Linux operating system (use UCT HPC server if needed)

**Software**
- **PLINK** (for genome-wide association analysis):
Download: http://pngu.mgh.harvard.edu/_purcell/plink/download.shtml
- **SMARTPCA.pl** (for principal components analysis / population structure correction):
Download: http://genepath.med.harvard.edu/~reich/Software.htm
- **Statistical software** for downstream analysis and visualization, such as:
   - **R**: http://cran.r-project.org/

### 10.1. Overview <a name="overview"></a>
Genome-wide association studies (GWAS) and candidate-gene studies rely on high-quality genotype data. Poorly genotyped samples or SNPs can create **systematic bias**, leading to false-positive or false-negative associations.

This module introduces a practical quality control (QC) workflow for **case–control genetic association studies**, focusing on:
- Detecting and removing low-quality **individuals**
- Detecting and removing low-quality **markers (SNPs)**
- Identifying **population structure** and **relatedness** that can bias results
- Using standard tools such as **PLINK** and **SMARTPCA**

In practice, this QC pipeline is completed before any formal association testing and can be generalized, with a few changes, to quantitative traits.

### 10.2. Learning Outcomes <a name="learn-obj"></a>
By the end of this QC session, students should be able to:
1.	Explain why QC is essential in case–control genetic association studies.
2.	Describe the main per-individual and per-marker QC steps in GWAS.
3.	Interpret key QC metrics such as call rate, heterozygosity, Hardy–Weinberg equilibrium, minor allele frequency (MAF), and identity-by-descent (IBD).
4.	Recognise how population stratification and relatedness can bias association results.
5.	Outline how tools such as PLINK are used to implement QC.

### 10.3. Why QC Matters <a name="why-qc-markers"></a>
Even if cases and controls are well matched and genotyping is performed in a good laboratory, several problems can still occur:
- **Genotype calling errors** (e.g. poor clustering of intensity data)
- **Low DNA quality or concentration**, leading to high missingness
- **Sample mix-ups or mis-recorded sex**
- **Unrecognised relatedness** among participants
- **Population stratification**, where allele frequencies differ by ancestry rather than disease status

These issues can:
•	inflate false positives (spurious associations)
•	mask true signals (false negatives)

A standardised QC pipeline helps identify substandard **samples** and **markers** and remove them before association testing.

### 10.4. Overview of the QC Workflow <a name="qc-work"></a>
For GWAS, QC is usually performed in two stages:

#### 10.4.1.	Per-individual QC (sample-level) <a name="sample-qc-over"></a>
- Check concordance of **reported vs genetic sex**
- Identify individuals with **high missingness** or **outlying heterozygosity**
- Detect **duplicates or close relatives**
- Detect **ancestry outliers** (population stratification)

#### 10.4.2.	Per-marker QC (SNP-level) <a name="snp-qc-over"></a>
- Remove SNPs with **high missingness**
- Remove SNPs with **strong deviation from Hardy–Weinberg equilibrium** (in controls)
- Remove SNPs with **different call rates in cases vs controls**
- Remove SNPs with **very low MAF**

The general principle is:
- **First clean individuals, then clean SNPs**: Removing problematic samples first avoids discarding otherwise good SNPs that only look bad because of a few low-quality individuals.

#### 10.5. Per-Individual QC <a name="sample-qc"></a>
##### 10.5.1 Sex checks
**Goal**: identify discrepancies between recorded sex and genetic sex.
- Use **X-chromosome genotypes** to compute a homozygosity/heterozygosity measure.
- Males (one X chromosome) should appear almost completely homozygous; females should show substantial heterozygosity.
- Samples with genetic sex that does not match recorded sex may indicate:
   - sample mis-labelling
   - data entry error
   - or, rarely, biological sex differences

**Action**: Investigate. If discrepancies cannot be resolved, exclude the sample.

##### 10.5.2 Missingness and heterozygosity
**Goal**: remove samples with poor genotype quality.
- **Missing genotype rate per individual**: high values often reflect low DNA quality.
- **Heterozygosity rate**:
   - unusually **high** heterozygosity may indicate contamination (mixed DNA)
   - unusually **low** heterozygosity may suggest inbreeding or technical problems

Typical practice:
- Inspect distributions across all individuals.
- Remove individuals with:
   - missingness above a chosen threshold (e.g. ≥3–7%), and/or
   - heterozygosity more than a few standard deviations from the mean.

##### 10.5.3 Duplicates and related individuals
**Goal**: ensure a population-based case–control sample with (approximately) unrelated individuals.
- Calculate pairwise relatedness using genome-wide SNPs, typically via **identity-by-descent (IBD)** measures.
- Use a pruned set of SNPs (remove regions of high LD such as the HLA).
- Expected IBD:
   - ≈1.0 for duplicates/monozygotic twins
   - ≈0.5 first-degree
   - ≈0.25 second-degree
   - ≈0.125 third-degree

**Action**:
- Treat pairs with IBD close to 1 as duplicates and remove one sample.
- Often, any pair with IBD above ~0.1875 (between second- and third-degree) is considered too closely related, and one individual is removed.

##### 10.5.4 Ancestry and population stratification
**Goal**: identify individuals whose ancestry differs substantially from the main study population.
- Apply **principal components analysis (PCA)** or a similar method (e.g. SMARTPCA) using:
   - a pruned set of genome-wide markers
   - reference samples of known ancestry (e.g. HapMap, 1000 Genomes)
- Plot the first few principal components to visualize clusters by ancestry.

**Action**:
- Remove individuals who cluster far from the main study group or fall with clearly different ancestry groups.
- Fine-scale structure can be handled later in association testing by adjusting for principal components as covariates.

#### 10.6. Per-Marker QC <a name="snp-qc"><a/>
After removing low-quality individuals, QC focuses on the SNPs themselves.

##### 10.6.1 SNP missingness
**Goal**: remove markers that frequently fail to genotype.
- Compute call rate (1 – missingness) for each SNP.
- Plot the distribution of missingness across all SNPs.
Typical thresholds:
Exclude SNPs with call rate below 95–97%,
or stricter thresholds (e.g. 99%) for low-frequency variants.

##### 10.6.2 Hardy–Weinberg equilibrium (HWE)
**Goal**: detect markers with gross genotyping errors.
- Test each SNP for deviation from HWE in controls only.
- Strong deviation may indicate poor genotype clustering or technical artefacts.

Typical approach:
- Choose a p-value threshold (e.g. p < 1×10⁻6 or more stringent).
- Remove SNPs that fail, while manually inspecting cluster plots for any variants of special interest.

Note: genuine disease-associated variants can deviate from HWE in **cases**, so only controls are used for QC.

##### 10.6.3 Differential missingness between cases and controls
**Goal**: avoid artefacts where genotypes are missing more often in one group.
- Compare call rates between cases and controls for each SNP.
- A significant difference can create spurious association signals.

**Action**:
- Remove SNPs with evidence of differential missingness (using a pre-specified p-value threshold).

##### 10.6.4 Minor allele frequency (MAF)
**Goal**: avoid unstable results driven by very rare variants in standard GWAS frameworks.
- Remove SNPs with very low MAF (typically <1%), especially in modest sample sizes.
- Rare variants are harder to call accurately and offer low power in single-variant tests.

# 11. References <a name="refs"></a>
1.	[Marees A. T. et al., Nature Protocols, 2020](https://pubmed.ncbi.nlm.nih.gov/29484742/)
2.	[Visscher P. M. et al., Nature Reviews Genetics, 2021](https://www.nature.com/articles/s43586-021-00056-9).
3.	PLINK 1.9 and Regenie documentation pages.
4.	GWAS QC **Adapted from the tutorial by**: [Anderson et al., Nature Protocols 2010](https://www.nature.com/articles/nprot.2010.116)

