# Module 6.5: Meta-analysis in Genome-Wide Association Studies (GWAS)

## Table of Contents <a name="toc"></a>
- [Learning objectives](#lobj)
- [Expected learning outcomes](#expo)
- [Introduction](#intro)
- [Why do we need meta-analysis in GWAS?](#why-meta)
- [Basic idea of GWAS meta-analysis](#basic-id)
- [Setting up a GWAS meta-analysis: consortia and bias](#setup)
- [Information needed from each study](#info)
- [Heterogeneity: when studies disagree](#hetero)
- [Fixed-effect vs random-effects meta-analysis](#fere-meta)
- [Bayesian and sequential perspectives](#bayes)
- [Example: Type 2 diabetes (T2D)](#example)
- [Future directions and challenges](#future-dir)
- [Key points from the module](#key-points)


## Learning objectives <a name="lobj"></a>

In this session, participants will:

- Understand why meta-analysis is needed in the context of GWAS and
  polygenic traits.

- Become familiar with how GWAS meta-analyses are set up (consortia,
  protocols, and data sharing).

- Learn the key pieces of information required from each contributing
  study (study design, QC, analysis strategy, genome build, imputation,
  etc.).

- Understand how summary statistics are combined across studies, and the
  difference between fixed-effect and random-effects models.

- Recognize the concepts of heterogeneity and inconsistency across
  studies, and how these are quantified and interpreted.

- Appreciate the role of imputation and harmonization when combining
  data from different genotyping platforms.

- See, through real examples (e.g. type 2 diabetes), how meta-analysis
  enables discovery of loci that would be underpowered in single
  studies.

## Expected learning outcomes <a name="expo"></a>

By the end of this section, participants should be able to:

- Explain, in their own words, why meta-analysis is essential for
  detecting small genetic effects in GWAS.

- Describe the main steps in conducting a GWAS meta-analysis, from
  individual study QC to combined analysis of summary statistics.

- List the core study-level information that must be checked and
  harmonised before combining datasets (phenotype definitions,
  covariates, population structure, QC thresholds, genome build, strand,
  imputation quality).

- Distinguish between fixed-effect and random-effects models, and state
  when each might be appropriate in a GWAS context.

- Interpret basic heterogeneity metrics (Q and I²) and explain possible
  reasons for differences in effect sizes between studies.

- Outline how imputed genotypes allow different GWAS platforms to be
  combined on a common set of variants.

- Critically assess published GWAS meta-analyses with respect to design,
  QC, handling of heterogeneity, and replication strategy.

## Introduction <a name="intro"></a>

Genome-wide association studies (GWAS) have made it possible to find
many common genetic variants associated with complex diseases and
traits. However, the effect sizes of these variants are usually small.
Even large single GWAS are often underpowered to detect these effects
with confidence.

**Meta-analysis** is a statistical approach that combines results from
multiple GWAS to:

- Increase power to detect true associations

- Assess how consistent genetic effects are across different studies and
  populations

- Discover additional loci that would be missed by any single study

## Why do we need meta-analysis in GWAS? <a name="why-meta"></a>

For common complex traits:

- Most risk variants have **modest odds ratios** (often around
  1.1-1.3).

- Even well-powered single studies may detect only a fraction of the
  true signals.

- Different studies often use different genotyping platforms and sample
  designs.

By combining data across studies:

- The **effective sample size** increases.

- We can examine **more variants**, including imputed SNPs not directly
  typed in every study.

- We can test whether genetic effects are **consistent or
  heterogeneous** across datasets.

In practice, GWAS meta-analysis has become the standard route for moving
from early "hits" to a robust, well-supported set of loci for a given
trait.

## Basic idea of GWAS meta-analysis <a name="basic-id"></a>

Most GWAS meta-analyses combine **summary statistics**, not raw
genotypes:

- Each study runs its own GWAS using a pre-defined analysis plan.

- For each SNP, the study provides:

  - An effect estimate (e.g. odds ratio or beta)

  - A measure of uncertainty (standard error or confidence interval)

  - A p-value

- The meta-analysis then combines these estimates across all studies.

Modern practice also uses **imputation**:

- Different platforms type different sets of SNPs.

- Using reference panels (e.g. HapMap, 1000 Genomes), each study imputes
  genotypes at untyped variants.

- This allows meta-analysis on a common set of millions of SNPs, even if
  each study typed a different subset.

Meta-analysis can be **cumulative**: new studies can be added over time
as they become available, updating the evidence for each SNP.

## Setting up a GWAS meta-analysis: consortia and bias <a name="setup"></a>

Most large GWAS meta-analyses are carried out within **consortia**:

- Multiple research groups agree in advance to combine their data.

- They harmonize phenotypes, covariates, and analysis protocols as far
  as possible.

- They often have access to full GWAS summary results for all members.

Key concern: **selection bias**

- If only "interesting" or "positive" results are shared or published,
  the meta-analysis will be biased.

- To minimize this, consortia aim to:

  - Use **pre-specified protocols**.

  - Share **complete genome-wide results** (not only top hits).

  - Apply common quality control and analysis rules.

Public repositories (such as dbGaP) also support data-sharing for GWAS,
although they come with additional requirements around privacy and data
use.


[Back to top](#toc)


## Information needed from each study <a name="info"></a>

Before combining results, the meta-analysis team needs to understand and
standardize multiple aspects of each contributing study.

1.  **Study design and population**

    - Case-control or cohort design

    - Inclusion/exclusion criteria

    - Recruitment methods, ancestry, and potential sources of bias

2.  **Quality control**

> Typical checks include:

- Hardy-Weinberg equilibrium

- Call/missingness rates

- Minor allele frequency thresholds

- For imputed SNPs: imputation quality metrics

> Thresholds (e.g. HWE p-value cut-offs, call rate ≥ 95%) are usually
> defined in advance and applied consistently across studies.

3.  **Analysis methods and covariates**

    - Phenotype definitions (e.g. how "case" and "control" are defined)

    - Covariate adjustments (age, sex, ancestry principal components,
      etc.)

    - Genetic model (e.g. additive on the log-odds scale)

> Ideally, all studies use **harmonized definitions and covariates**.
> Where this is not possible, differences must be documented and their
> impact considered.

4.  **Independence and relatedness**

    - Within each study, related individuals must be handled properly
      (e.g. mixed models or family-based methods).

    - Population stratification must be controlled (e.g. via principal
      components or genomic control).

    - Overlap between studies (shared samples) must be identified; if
      present, this covariance needs to be accounted for.

5.  **Genome build, strand and variant representation**

    - All studies should report positions using the **same genome
      build** (e.g. hg19).

    - Strand issues (A/T and C/G SNPs) must be resolved so that the
      **same allele** is used as the effect allele across studies.

    - Directly genotyped vs imputed calls should be clearly labeled, and
      only imputed SNPs above a quality threshold should be included.

A failure at any of these steps can easily create spurious associations
or obscure true signals.

## Heterogeneity: when studies disagree <a name="hetero"></a>

Meta-analysis does not assume that all studies show exactly the same
effect. Instead, we can quantify and interpret **heterogeneity**.

Common measures:

- **Cochran's Q**

  - Tests whether there is more variability between study results than
    expected by chance.

  - Low power with few studies; almost always significant when many
    studies are combined.

- **I² statistic**

  - Estimates the **percentage** of total variation due to between-study
    differences rather than sampling error.

  - Rough interpretation:

    - \~0-25%: low heterogeneity

    - \~25-75%: moderate

    - 75%: high

- **Between-study variance (τ²)**

  - An absolute measure of how much true effect sizes vary across
    studies.

Possible reasons for heterogeneity include:

- Real biological differences between populations (e.g. differing LD
  patterns, environmental contexts)

- Differences in phenotype definitions or covariate adjustment

- Genotyping or imputation errors in some studies

- Pure chance fluctuations

Heterogeneity can be scientifically informative (e.g. pointing to
gene-environment interaction), but it also limits how widely a summary
effect can be generalized.

## Fixed-effect vs random-effects meta-analysis <a name="fere-meta"></a>

Two main models are used when combining effect sizes:

1.  **Fixed-effect model**

    - Assumes that there is **one true effect** shared by all studies.

    - Differences in study estimates arise only from sampling error.

    - More powerful when this assumption is reasonable.

    - Often used when the primary goal is to test "Is there any
      association?" at genome-wide scale.

2.  **Random-effects model**

    - Assumes that effects **vary across studies**, following some
      distribution.

    - Incorporates the between-study variance (τ²).

    - Produces wider confidence intervals when heterogeneity is present.

    - More conservative for hypothesis testing, but better reflects true
      uncertainty when effects really differ by study.

In the absence of heterogeneity, the two models give the same result.
When heterogeneity is substantial, random-effects models are more
cautious.

In GWAS meta-analysis:

- Many consortia report fixed-effect meta-analysis as the primary scan,
  often supplemented by heterogeneity testing and, sometimes,
  random-effects estimates for follow-up of specific loci.

## Bayesian and sequential perspectives <a name="bayes"></a>

Because GWAS meta-analysis is often **cumulative**, it fits naturally
into Bayesian and sequential frameworks:

- Bayesian meta-analysis treats earlier data as a **prior**, and each
  new study updates this to a **posterior**.

- Different prior assumptions can be explored (e.g. on effect size or
  heterogeneity), and robustness can be assessed.

- Sequential or "trial-style" approaches adjust significance thresholds
  depending on how many times the accumulating data have been examined.

In practice, most GWAS consortia still rely on frequentist fixed-effect
scans with predefined genome-wide significance thresholds, sometimes
complemented with simple Bayesian credibility metrics.

## Example: Type 2 diabetes (T2D) <a name="example"></a>

A classic example of GWAS meta-analysis is T2D:

- Three independent T2D GWAS (e.g. WTCCC, DGI, FUSION) were combined.

- Each study:

  - Imputed SNPs using a common reference panel.

  - Applied stringent QC for both typed and imputed SNPs.

  - Corrected for population stratification.

**Stage 1** - Meta-analysis of GWAS datasets:

- \~10,000 individuals

- \~2.2 million SNPs (typed or imputed)

- Produced a ranked list of association signals.

**Stage 2** - Follow-up of promising SNPs:

- Dozens of top SNPs were genotyped in \~20,000 additional samples.

- More associated loci emerged than expected by chance, supporting true
  signal enrichment.

**Stage 3** - Large-scale replication:

- A smaller set of the most promising SNPs was tested in \>50,000
  additional individuals.

- Several loci reached **genome-wide significance** only after this
  combined effort.

> Key observations:

- All T2D risk variants discovered had **small effect sizes** (allelic
  ORs around 1.1-1.15).

- No single study alone had enough power to reach genome-wide
  significance for most loci.

- Meta-analysis, followed by staged replication, was essential for
  building the current catalogue of T2D loci.

Even after meta-analysis, further work (fine-mapping, sequencing,
functional studies) is needed to identify the true causal variants and
understand their biology.


[Back to top](#toc)


## Future directions and challenges <a name="future-dir"></a>

Meta-analysis has become a central tool for GWAS, but several challenges
remain:

- **Missing heritability**\
  Even large GWAS meta-analyses explain only a modest fraction of the
  genetic variance for most traits.

- **Extension to diverse populations**\
  Most early work focused on European-ancestry cohorts. There is a
  growing need for:

  - GWAS and meta-analyses in non-European populations

  - Careful handling of differing LD patterns and allele frequencies

- **More complex data**\
  Future meta-analyses may need to integrate:

  - Rare variants (from sequencing data)

  - Structural variants and copy number variation

  - Multi-omics data (e.g. gene expression, methylation)

  - Gene-environment interactions

- **Clinical translation**\
  For risk prediction and "polygenic scores", we need:

  - Precise estimates of effect sizes

  - Consistent performance across populations and settings

  - Large, well-characterized validation cohorts

- **Data integration at scale**\
  As more GWAS, consortia, and biobanks appear, we will move towards:

  - "Meta-meta-analyses" that combine results across multiple consortia

  - Continuous updating of the evidence base for genetic associations.

## Key points from the module <a name="key-points"></a>

- GWAS meta-analysis is essential because **effect sizes are small**,
  and single studies are underpowered.

- It combines **summary statistics** from multiple studies, often using
  imputation to harmonize variant sets.

- Success depends on **rigorous QC, harmonized definitions, and careful
  control of bias and population structure**.

- Heterogeneity is expected and must be quantified and interpreted
  rather than ignored.

- Fixed-effect models are widely used for discovery, with heterogeneity
  and random-effects analyses used for interpretation.

- Real-world examples (e.g. type 2 diabetes) show that meta-analysis is
  often the only way to reach genome-wide significance for many loci.

- Future work will extend these methods to more diverse populations,
  more complex data types, and clinical applications.

[Back to top](#toc)

