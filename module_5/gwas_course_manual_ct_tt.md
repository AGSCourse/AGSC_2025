# GWAS Course Manual
Developed by Chris Kintu and Tsaone Tamuhla

## What is GWAS?
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

## Terminologies

| Term | Definition |
| :-: | :-: |
| Confounding | A type of bias in statistical analysis causing spurious or distorted findings caused by a correlation between an extraneous variable (the confounding variable) and both the dependent/exposure variable (e.g. the genotype at a given locus) and the independent/outcome variable (e.g. case- control status). |
| Failure rate | The proportion of missing genotypes. Genotypes are classified as missing if the genotype-calling algorithm cannot infer the genotype with sufficient confidence. Can be calculated across each individual and/or SNP. |
| False-negative | Occurs when a true disease-associated variant is not associated with disease in a given study. |
| False-positive | Occurs when a variant not in truth associated with disease status is significantly associated with disease in a given study. |
| Genotype calling algorithm | A statistical algorithm that, per marker and per individual, converts intensity data from two allelic probes into a single genotype for analysis. |

