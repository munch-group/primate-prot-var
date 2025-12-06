# gnomAD Variant Fields Reference

This document describes all fields retrieved by `scripts/primate_aa_variants.py` when using the `--method gnomad` option. The script queries the gnomAD v4 GraphQL API for missense variants with population allele frequency above the specified threshold (default: 10%).

## Table of Contents

1. [Basic Variant Information](#basic-variant-information)
2. [Amino Acid Change](#amino-acid-change)
3. [Global Allele Frequencies](#global-allele-frequencies)
4. [Filtering Allele Frequencies](#filtering-allele-frequencies)
5. [Population-Specific Frequencies](#population-specific-frequencies)
6. [In Silico Pathogenicity Predictors](#in-silico-pathogenicity-predictors)
7. [Transcript-Level Predictions](#transcript-level-predictions)
8. [Loss-of-Function Annotations](#loss-of-function-annotations)
9. [Quality Flags and Filters](#quality-flags-and-filters)

---

## Basic Variant Information

| Field | Type | Description |
|-------|------|-------------|
| `species` | string | Species identifier. Always "human" for gnomAD data. |
| `gene_id` | string | Ensembl gene ID (e.g., "ENSG00000012048" for BRCA1). |
| `gene_name` | string | HGNC gene symbol (e.g., "BRCA1"). |
| `transcript_id` | string | Ensembl transcript ID for the canonical transcript (e.g., "ENST00000357654"). |
| `variant_id` | string | gnomAD variant identifier in format `chrom-pos-ref-alt` (e.g., "17-43071077-T-C"). |
| `rsid` | string | dbSNP reference SNP ID (e.g., "rs1799966"). Null if not in dbSNP. |
| `chromosome` | string | Chromosome number or letter (e.g., "17", "X"). |
| `genomic_position` | integer | Genomic position on GRCh38 reference assembly. |
| `ref_allele` | string | Reference allele sequence. |
| `alt_allele` | string | Alternate (variant) allele sequence. |

---

## Amino Acid Change

| Field | Type | Description |
|-------|------|-------------|
| `aa_position` | integer | Position of the amino acid change in the protein sequence (1-based). |
| `ref_aa` | string | Reference amino acid at this position. Can be 1-letter (e.g., "S") or 3-letter code (e.g., "Ser"). |
| `alt_aa` | string | Alternate amino acid resulting from the variant. Can be 1-letter or 3-letter code. |
| `hgvsp` | string | HGVS protein notation (e.g., "p.Ser1613Gly"). Describes the amino acid change using standard nomenclature. |
| `hgvsc` | string | HGVS coding DNA notation (e.g., "c.4837A>G"). Describes the nucleotide change relative to the coding sequence. |
| `consequence` | string | Variant consequence type. For this script, always "missense_variant" since only missense variants are retrieved. |

---

## Global Allele Frequencies

These fields describe the overall allele frequency across all samples in gnomAD.

| Field | Type | Description |
|-------|------|-------------|
| `allele_frequency` | float | Maximum allele frequency between exome and genome datasets. This is the primary frequency used for filtering. Range: 0-1. |
| `exome_af` | float | Allele frequency in the exome dataset. Null if variant not observed in exomes. |
| `genome_af` | float | Allele frequency in the genome dataset. Null if variant not observed in genomes. |
| `exome_ac` | integer | Allele count in exomes (number of times the alternate allele was observed). |
| `exome_an` | integer | Allele number in exomes (total number of alleles genotyped at this site). |
| `exome_ac_hom` | integer | Number of homozygous individuals for the alternate allele in exomes. |
| `exome_ac_hemi` | integer | Number of hemizygous individuals (for X/Y chromosome variants in males) in exomes. |
| `genome_ac` | integer | Allele count in genomes. |
| `genome_an` | integer | Allele number in genomes. |
| `genome_ac_hom` | integer | Number of homozygous individuals in genomes. |
| `genome_ac_hemi` | integer | Number of hemizygous individuals in genomes. |

**Note**: Allele frequency is calculated as `AF = AC / AN`.

---

## Filtering Allele Frequencies

Filtering allele frequencies (FAF) provide conservative frequency estimates that account for sampling uncertainty. These are useful for clinical variant filtering.

| Field | Type | Description |
|-------|------|-------------|
| `faf95_max` | float | Maximum 95% confidence interval lower bound for allele frequency across populations. Conservative estimate accounting for sampling error. |
| `faf95_max_population` | string | Genetic ancestry group with the highest faf95 value (e.g., "nfe", "afr"). |
| `faf99_max` | float | Maximum 99% confidence interval lower bound for allele frequency across populations. Even more conservative than faf95. |
| `faf99_max_population` | string | Genetic ancestry group with the highest faf99 value. |
| `faf95_popmax` | float | Population maximum filtering allele frequency at 95% CI. Legacy field from earlier gnomAD versions. |
| `faf95_popmax_population` | string | Population with the highest faf95_popmax value. |

**Usage**: FAF values are commonly used in clinical genetics to filter out common variants. A variant with `faf95_max > 0.0001` is typically considered too common to cause a rare Mendelian disease.

---

## Population-Specific Frequencies

Allele frequencies calculated separately for each genetic ancestry group. These are computed from allele counts (AC/AN) for each population.

| Field | Type | Description |
|-------|------|-------------|
| `af_afr` | float | Allele frequency in African/African American population. |
| `af_amr` | float | Allele frequency in Admixed American/Latino population. |
| `af_asj` | float | Allele frequency in Ashkenazi Jewish population. |
| `af_eas` | float | Allele frequency in East Asian population. |
| `af_fin` | float | Allele frequency in Finnish population (separated from NFE due to founder effects). |
| `af_nfe` | float | Allele frequency in Non-Finnish European population. |
| `af_sas` | float | Allele frequency in South Asian population. |
| `af_mid` | float | Allele frequency in Middle Eastern population. |
| `af_ami` | float | Allele frequency in Amish population. |

**Note**: Null values indicate the variant was not observed in that population or the population was not well-represented at that site.

---

## In Silico Pathogenicity Predictors

Computational predictions of variant pathogenicity from various algorithms. These scores help prioritize variants for functional follow-up.

| Field | Type | Range | Description |
|-------|------|-------|-------------|
| `cadd_phred` | float | 0-60 | **CADD** (Combined Annotation Dependent Depletion) phred-scaled score. Integrates multiple annotations into a single deleteriousness score. Higher scores indicate more deleterious variants. Scores ≥20 represent the top 1% most deleterious substitutions; ≥30 represents top 0.1%. |
| `revel_score` | float | 0-1 | **REVEL** (Rare Exome Variant Ensemble Learner) score. Ensemble method combining 13 individual prediction tools specifically trained on rare missense variants. Higher scores indicate higher probability of pathogenicity. Suggested threshold: >0.5 for likely pathogenic. |
| `spliceai_score` | float | 0-1 | **SpliceAI** delta score (maximum across splice sites). Deep learning prediction of splicing impact. Higher scores indicate higher probability of affecting splicing. Thresholds: >0.2 (likely affects splicing), >0.5 (high confidence), >0.8 (very high confidence). |
| `alphamissense_score` | float | 0-1 | **AlphaMissense** pathogenicity score from DeepMind. Uses protein structure and evolutionary information. Thresholds: <0.34 (likely benign), 0.34-0.564 (ambiguous), >0.564 (likely pathogenic). |
| `pangolin_score` | float | -1 to 1 | **Pangolin** splice effect prediction. Positive scores indicate potential splice-altering effects. |
| `phylop_score` | float | -20 to 9.28 | **PhyloP** conservation score based on 241 placental mammals (Zoonomia). Positive scores indicate conservation (slower evolution than expected); negative scores indicate acceleration. Higher positive values = more conserved. |
| `sift_max` | float | 0-1 | **SIFT** (Sorting Intolerant From Tolerant) maximum score. Predicts whether amino acid substitution affects protein function based on sequence homology. Lower scores = more damaging. Threshold: <0.05 is "deleterious". |
| `polyphen_max` | float | 0-1 | **PolyPhen-2** (Polymorphism Phenotyping v2) maximum score. Predicts impact of amino acid substitution on protein structure/function. Higher scores = more damaging. Thresholds: >0.85 (probably damaging), 0.15-0.85 (possibly damaging), <0.15 (benign). |

### Interpretation Guidelines

| Predictor | Benign | Uncertain | Pathogenic |
|-----------|--------|-----------|------------|
| CADD | <10 | 10-20 | >20 |
| REVEL | <0.25 | 0.25-0.5 | >0.5 |
| AlphaMissense | <0.34 | 0.34-0.564 | >0.564 |
| SpliceAI | <0.1 | 0.1-0.5 | >0.5 |
| SIFT | >0.05 | - | <0.05 |
| PolyPhen | <0.15 | 0.15-0.85 | >0.85 |

---

## Transcript-Level Predictions

Categorical predictions from SIFT and PolyPhen for the specific transcript.

| Field | Type | Values | Description |
|-------|------|--------|-------------|
| `polyphen_prediction` | string | "benign", "possibly_damaging", "probably_damaging" | PolyPhen-2 categorical prediction for this transcript. Based on HumDiv model trained on human disease-causing mutations. |
| `sift_prediction` | string | "tolerated", "deleterious" | SIFT categorical prediction. "deleterious" corresponds to SIFT score <0.05. |

**Note**: These may be null even when `sift_max` or `polyphen_max` have values, as the categorical predictions are transcript-specific.

---

## Loss-of-Function Annotations

LOFTEE (Loss-Of-Function Transcript Effect Estimator) annotations. While this script retrieves missense variants, some may have LoF annotations if they affect splicing or have other impacts.

| Field | Type | Values | Description |
|-------|------|--------|-------------|
| `lof` | string | "HC", "LC", null | Loss-of-function confidence. "HC" = High Confidence (passes all LOFTEE filters). "LC" = Low Confidence (fails one or more filters). Null for non-LoF variants. |
| `lof_filter` | string | various | LOFTEE filter(s) that the variant failed. Examples: "END_TRUNC" (truncation in last 5% of transcript), "PHYLOCSF_WEAK" (weak conservation), "NON_CAN_SPLICE" (non-canonical splice site). |
| `lof_flags` | string | various | Additional LOFTEE flags. Examples: "SINGLE_EXON" (gene has single exon), "NAGNAG_SITE" (NAGNAG splice site), "PHYLOCSF_UNLIKELY_ORF" (unlikely to be protein-coding). |

---

## Quality Flags and Filters

Quality control information about the variant call.

| Field | Type | Description |
|-------|------|-------------|
| `flags` | string | Comma-separated list of variant-level flags. Examples: "lcr" (low complexity region), "segdup" (segmental duplication), "lcr" (low complexity region). Null if no flags. |
| `exome_filters` | string | Comma-separated list of filters applied to exome data. "PASS" or null indicates the variant passed all filters. Other values indicate QC failures (e.g., "AC0", "InbreedingCoeff"). |
| `genome_filters` | string | Comma-separated list of filters applied to genome data. Same interpretation as exome_filters. |
| `source` | string | Data source identifier. Always "gnomAD_v4" for this script. |

### Common Filter Values

| Filter | Description |
|--------|-------------|
| `PASS` | Variant passed all quality filters |
| `AC0` | Allele count is zero after filtering genotypes |
| `InbreedingCoeff` | Excess heterozygosity (inbreeding coefficient < -0.3) |
| `RF` | Failed random forest filtering |
| `AS_VQSR` | Failed allele-specific VQSR filtering |

---

## Example Output

```csv
species,gene_id,gene_name,transcript_id,variant_id,rsid,chromosome,genomic_position,ref_allele,alt_allele,aa_position,ref_aa,alt_aa,hgvsp,hgvsc,consequence,allele_frequency,exome_af,genome_af,exome_ac,exome_an,exome_ac_hom,exome_ac_hemi,genome_ac,genome_an,genome_ac_hom,genome_ac_hemi,faf95_max,faf95_max_population,faf99_max,faf99_max_population,faf95_popmax,faf95_popmax_population,af_afr,af_amr,af_asj,af_eas,af_fin,af_nfe,af_sas,af_mid,af_ami,cadd_phred,revel_score,spliceai_score,alphamissense_score,pangolin_score,phylop_score,sift_max,polyphen_max,polyphen_prediction,sift_prediction,lof,lof_filter,lof_flags,flags,exome_filters,genome_filters,source
human,ENSG00000012048,BRCA1,ENST00000357654,17-43071077-T-C,rs1799966,17,43071077,T,C,1613,Ser,Gly,p.Ser1613Gly,c.4837A>G,missense_variant,0.339841,0.339841,0.317829,496800,1461858,86242,0,48324,152044,8007,0,,,,,,,0.233124,0.320723,0.367922,0.353368,0.396417,0.327411,0.499258,0.37344,,8.64,0.252,0.0,,-0.01,0.012,0.11,0.038,,,,,,,,,gnomAD_v4
```

---

## References

- [gnomAD Browser](https://gnomad.broadinstitute.org/)
- [gnomAD v4 Release Notes](https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/)
- [CADD](https://cadd.gs.washington.edu/)
- [REVEL](https://sites.google.com/site/revelgenomics/)
- [SpliceAI](https://github.com/Illumina/SpliceAI)
- [AlphaMissense](https://alphamissense.hegelab.org/)
- [LOFTEE](https://github.com/konradjk/loftee)
- [SIFT](https://sift.bii.a-star.edu.sg/)
- [PolyPhen-2](http://genetics.bwh.harvard.edu/pph2/)
