# Exploratory Data Analysis Report
## Genomics of Drug Sensitivity in Cancer (GDSC) Dataset

| Field | Details |
|---|---|
| **Report Type** | Exploratory Data Analysis (EDA) |
| **Dataset** | Genomics of Drug Sensitivity in Cancer (GDSC) |
| **Source Files** | `data_exp.ipynb`, `notebook.ipynb` |
| **Language / Tools** | Python 3, pandas, matplotlib, seaborn, NumPy |
| **Report Date** | 6 May 2026 |
| **Status** | Draft — Internship Stage One |

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Background and Objectives](#2-background-and-objectives)
3. [Dataset Overview](#3-dataset-overview)
4. [Data Preprocessing and Quality Assessment](#4-data-preprocessing-and-quality-assessment)
5. [Exploratory Data Analysis](#5-exploratory-data-analysis)
   - 5.1 [Univariate Analysis — Continuous Variables](#51-univariate-analysis--continuous-variables)
   - 5.2 [Univariate Analysis — Categorical Variables](#52-univariate-analysis--categorical-variables)
   - 5.3 [Multivariate Analysis](#53-multivariate-analysis)
6. [Drug Sensitivity Patterns](#6-drug-sensitivity-patterns)
7. [Cancer Cell Line Analysis](#7-cancer-cell-line-analysis)
8. [Genomic, Transcriptomic, and Epigenomic Influence on Drug Response](#8-genomic-transcriptomic-and-epigenomic-influence-on-drug-response)
9. [Key Findings and Observations](#9-key-findings-and-observations)
10. [Limitations and Future Work](#10-limitations-and-future-work)
11. [Appendix: Variable Glossary](#11-appendix-variable-glossary)

## 1. Executive Summary

This report documents the exploratory data analysis (EDA) conducted on the **Genomics of Drug Sensitivity in Cancer (GDSC)** dataset as part of an AI for Genomics internship, Stage One. The dataset contains 162,103 drug–cell-line interaction records covering 246 unique drugs, 737 unique cancer cell lines, and 23 target pathways.

The analysis progressed through data ingestion and quality assessment, followed by univariate and multivariate profiling of key variables, drug sensitivity characterisation, cancer-type-level comparisons, and an investigation into how genomic features — specifically copy number alterations, gene expression, methylation status, and microsatellite instability — relate to observed drug response. All significant preprocessing decisions, including how duplicate records were handled and how binary columns were standardised, are documented in detail throughout the report.

## 2. Background and Objectives

The GDSC project systematically screens cancer cell lines against a broad panel of anti-cancer compounds and records pharmacological sensitivity metrics for each combination. Understanding which drugs are most potent, which cancer types are most resistant, and how molecular features modulate drug response is fundamental to the development of precision oncology strategies.

This analysis was undertaken with the following primary objectives in mind. First, to characterise the structure and quality of the GDSC dataset before any modelling work begins. Second, to understand the distribution of key pharmacological metrics, namely LN_IC50, AUC, and Z-Score. Third, to rank drugs by both potency and response variability. Fourth, to compare drug sensitivity profiles across cancer types using the TCGA classification scheme. Finally, to assess whether the availability of genomic feature data — specifically CNA, Gene Expression, Methylation, and MSI — is systematically associated with differential drug response.

## 3. Dataset Overview

### 3.1 Source

The primary data file is `GDSC.xlsx`, which contains 162,103 rows and 19 columns. Each row represents a single drug–cell-line sensitivity experiment.

### 3.2 Variables

| Variable | Type | Description |
|---|---|---|
| `COSMIC_ID` | Integer (ID) | Unique cell line identifier in the COSMIC catalogue |
| `CELL_LINE_NAME` | String (ID) | Human-readable cell line name |
| `DRUG_ID` | Integer (ID) | Numeric drug identifier |
| `DRUG_NAME` | String (ID) | Drug common name |
| `LN_IC50` | Float | Natural log of the half-maximal inhibitory concentration |
| `AUC` | Float | Area Under the dose-response Curve |
| `Z_SCORE` | Float | Standardised sensitivity score |
| `TARGET` | String | Molecular target of the drug (183 unique values) |
| `TARGET_PATHWAY` | String | Pathway of the drug target (23 unique values) |
| `TCGA_DESC` | String | TCGA cancer-type label |
| `Cancer_Type_(matching_TCGA_label)` | String | Alternative TCGA cancer-type label |
| `GDSC_Tissue_descriptor_1` | String | Broad tissue system (17 unique values) |
| `GDSC_Tissue_descriptor_2` | String | Organ-level descriptor (737 unique values) |
| `Microsatellite_instability_Status_(MSI)` | Categorical | MSI status of the cell line |
| `CNA` | Binary categorical | Copy Number Alteration data availability (Y/N) |
| `Gene_Expression` | Binary categorical | Gene Expression data availability (Y/N) |
| `Methylation` | Binary categorical | Methylation data availability (Y/N) |
| `Screen_Medium` | Categorical | Cell culture medium used |
| `Growth_Properties` | Categorical | Cell growth characteristics |

### 3.3 Key Cardinalities

| Dimension | Count |
|---|---|
| Total entries (raw) | 162,103 |
| Unique cell lines (COSMIC_ID) | 737 |
| Unique drugs | 246 |
| Unique drug targets | 183 |
| Unique target pathways | 23 |
| Unique tissue systems | 17 |
| Maximum theoretical combinations (737 × 246) | 181,302 |
| Actual entries | 162,103 |
| Coverage | ~89.5% of all possible combinations |

The gap between the theoretical maximum of 181,302 combinations and the actual 162,103 entries tells us that not every drug was tested against every cell line. This partial coverage is worth keeping in mind when interpreting global drug rankings, as the absence of certain combinations may not be random.

## 4. Data Preprocessing and Quality Assessment

### 4.1 Missing Values

An initial null check using `gdsc.isnull().sum()` confirmed that there were no missing values across any of the 19 columns. No imputation was needed, which is a fortunate starting point given the size of the dataset.

### 4.2 Column Standardisation

Column names in the raw file contained leading or trailing whitespace and internal spaces, which would cause problems in downstream code. These were resolved early in the `data_exp.ipynb` notebook by stripping whitespace and replacing spaces with underscores:

```python
gdsc.columns = gdsc.columns.str.strip()
gdsc.columns = gdsc.columns.str.replace(" ", "_")
```

### 4.3 Data Type Corrections

The sensitivity columns `LN_IC50` and `AUC` were coerced to `float64` using `pd.to_numeric(..., errors='coerce')` to guard against any non-numeric entries that might have been present. A follow-up null check confirmed that no new missing values were introduced by this conversion.

### 4.4 Binary Column Standardisation

The genomic feature flag columns (`CNA`, `Gene_Expression`, `Methylation`) were found to store values in mixed formats, with some entries reading `Y` or `Yes` and others `N` or `No`. To ensure consistency, all three columns were standardised to uppercase stripped strings before being cast to the `category` dtype:

```python
binary_cols = ['CNA', 'Gene_Expression', 'Methylation']
for col in binary_cols:
    gdsc[col] = gdsc[col].astype(str).str.strip().str.upper().astype('category')
```

### 4.5 Duplicate Detection and Resolution

A check for duplicate combinations of `CELL_LINE_NAME` and `DRUG_NAME` revealed 4,290 repeated pairs. Rather than treating these as data entry errors and dropping them outright, a closer inspection in `notebook.ipynb` showed that the duplicated rows were consistent across all categorical columns, including MSI status, CNA, Gene Expression, and other metadata fields. This consistency suggested that the duplicates arose from the same cell line being screened under different experimental conditions across separate runs, rather than from any clerical mistake.

Given this, the decision was taken to retain the information from these repeated measurements by averaging the continuous sensitivity metrics within each unique cell line and drug combination:

```python
drug_response = (
    gdsc
    .groupby(['COSMIC_ID', 'DRUG_NAME'], as_index=False)
    .agg({'LN_IC50': 'mean', 'AUC': 'mean', 'Z_SCORE': 'mean'})
)
```

Cell line metadata was extracted separately, taking one record per `COSMIC_ID`, and then merged back to produce the clean analytical dataset `gdsc_clean`, in which every cell line and drug combination appears exactly once.

### 4.6 Column Redundancy Check

A comparison between the `TCGA_DESC` and `Cancer_Type_(matching_TCGA_label)` columns revealed that they are not identical across all rows. Both columns were retained for now pending further clarification on how they differ semantically, as using the wrong one as the canonical cancer-type label could affect downstream analyses.

## 5. Exploratory Data Analysis

### 5.1 Univariate Analysis — Continuous Variables

Three continuous variables were profiled: LN_IC50, AUC, and Z_SCORE. For each, a paired visualisation was produced consisting of a histogram with mean (blue dashed line) and median (red solid line) overlaid, and a boxplot positioned directly below to expose the spread and any outliers. A custom `plot_distribution()` function using `GridSpecFromSubplotSpec` was written to render these panels in a responsive, screen-size-aware layout.

LN_IC50 spans a notably wide range, reflecting the large variation in drug potency observed across cell lines. AUC values are bounded between 0 and 1, representing the proportion of cell viability retained across all tested concentrations. Z_SCORE, as expected for a standardised metric, centres around zero.

### 5.2 Univariate Analysis — Categorical Variables

Identifier-type columns such as `COSMIC_ID`, `DRUG_ID`, `DRUG_NAME`, `CELL_LINE_NAME`, and `TARGET` were excluded from categorical profiling, as their high cardinality makes bar chart visualisations uninformative. All remaining discrete variables were visualised using proportional bar charts with percentage labels annotated above each bar.

A custom `plot_proportions()` function was developed for this purpose, using `constrained_layout` to ensure robust multi-panel rendering. Variables profiled in this way included MSI status, CNA availability, Gene Expression availability, Methylation availability, Screen Medium, Growth Properties, TCGA cancer type, and Tissue Descriptor 1.

### 5.3 Multivariate Analysis

#### 5.3.1 Correlation Among Continuous Variables

A Pearson correlation heatmap was computed across LN_IC50, AUC, and Z_SCORE using seaborn's heatmap function with a `coolwarm` diverging palette centred at zero. This gives a first indication of the degree of linear association between the three sensitivity metrics and helps identify whether they are largely capturing the same signal or providing complementary information.

#### 5.3.2 Drug-Level Summary Statistics

A per-drug aggregation was performed to produce a summary table, `drug_stats`, containing the median, standard deviation, and observation count of LN_IC50, AUC, and Z_SCORE for each drug across all tested cell lines. This table forms the analytical basis for the drug sensitivity rankings discussed in Section 6.

## 6. Drug Sensitivity Patterns

### 6.1 Most Potent Drugs

Drugs were ranked by ascending median LN_IC50 across all cell lines. The ten drugs with the lowest median LN_IC50 represent the highest pharmacological potency, meaning they achieve 50% cell growth inhibition at the lowest concentrations. These compounds are of the greatest interest for downstream therapeutic investigation.

The results were displayed as horizontal bar charts with plus-or-minus one standard deviation error bars. The width of the error bars is informative in its own right: narrow bars indicate a drug that is consistently potent across most cell lines, while wide bars suggest potency concentrated within a sensitive subset, which may point towards a targetable biomarker.

### 6.2 Least Potent Drugs

By the same logic, the ten drugs with the highest median LN_IC50 were also identified. These compounds require substantially higher concentrations to achieve 50% inhibition and may represent agents with limited broad-spectrum efficacy, or drugs whose activity is restricted to a narrow set of molecular contexts not well-represented in the screened panel.

### 6.3 Drugs with Highly Variable Response

Beyond average potency, the ten drugs with the highest standard deviation in LN_IC50 across all cell lines were extracted. High variability in response is a biologically meaningful signal, suggesting that a drug's effectiveness is strongly context-dependent and likely linked to specific genomic alterations, cancer subtypes, or pathway activation states. Drugs in this category are strong candidates for biomarker discovery work, where the goal would be to identify the molecular features that predict sensitivity or resistance.

## 7. Cancer Cell Line Analysis

### 7.1 LN_IC50 Distribution by Cancer Type

A boxplot was generated for all TCGA cancer types, sorted by their median LN_IC50 from left (most sensitive) to right (most resistant). A horizontal reference line marking the overall dataset median provides a neutral benchmark against which each cancer type can be compared visually.

It is worth emphasising that this view aggregates drug response across all 246 drugs simultaneously. A cancer type that appears globally resistant may still harbour specific, targetable drug sensitivities that would only emerge from a drug-specific analysis. The global boxplot is therefore best understood as a high-level overview rather than a definitive statement about any cancer type's therapeutic intractability.

### 7.2 Drug–Cancer Type Sensitivity Heatmap

To move toward drug-specific comparisons, a pivot table heatmap was constructed using the median LN_IC50 for each cancer type and drug combination. To keep the visualisation interpretable, only the 20 drugs with the highest variability in LN_IC50 were included. The heatmap uses a diverging colour scale where green indicates sensitivity and red indicates resistance, making it straightforward to identify cancer-type-specific drug vulnerabilities at a glance.

### 7.3 Drug–Cell Line Sensitivity Heatmap

A second heatmap was produced at full cell line resolution, covering all 737 cell lines against the top 10 most variable drugs. This finer-grained view reveals the degree of heterogeneity that exists even within a single TCGA cancer classification and can surface outlier cell lines with unusual sensitivity or resistance profiles that might otherwise be obscured by type-level summaries.

## 8. Genomic, Transcriptomic, and Epigenomic Influence on Drug Response

### 8.1 Variables Assessed

Four genomic feature variables were considered in this section: CNA (Copy Number Alteration availability), Gene Expression (transcriptomic profiling availability), Methylation (DNA methylation data availability), and Microsatellite Instability Status (MSI).

### 8.2 Analysis Approach

A 4-by-3 grid of boxplots was produced, with each row corresponding to one genomic feature flag and each column corresponding to one sensitivity metric (LN_IC50, AUC, or Z_SCORE). Each panel shows the distribution of a sensitivity metric stratified by the binary or categorical value of that genomic feature. The complete figure was saved as `genomic_flags_vs_sensitivity.png` at 150 dpi for archival purposes.

### 8.3 Interpretation

It is important to note that the CNA, Gene Expression, and Methylation columns in this dataset indicate data availability rather than actual molecular measurements. As a result, any observed differences in drug sensitivity across these groups may partly reflect differences in the experimental cohorts or tissue types that tend to have particular data modalities available, rather than a direct biological effect of the molecular feature itself.

MSI status, on the other hand, is a well-established biological biomarker. MSI-high tumours arise from defective DNA mismatch repair and are known to respond differently to both chemotherapy and immunotherapy. The MSI panels therefore offer a more directly interpretable comparison. These visualisations represent a first-pass exploration; formal statistical testing will be needed before any conclusions can be drawn.

## 9. Key Findings and Observations

The analysis yielded several observations worth highlighting going into Stage Two.

The dataset covers approximately 89.5% of all theoretical drug–cell-line combinations, confirming a broad but deliberately incomplete experimental screen. No missing values were present in the raw data, which significantly reduced preprocessing complexity. However, 4,290 duplicate drug–cell-line records were identified and resolved by averaging sensitivity metrics, a decision that preserves information from repeated experiments while producing a clean analytical dataset.

The dataset's pharmacological breadth is substantial: 246 drugs acting across 183 molecular targets and 23 target pathways. Despite this breadth, drug potency as measured by LN_IC50 is highly variable across both compounds and cancer types, underscoring the context-dependent nature of drug efficacy and the importance of identifying predictive biomarkers.

The divergence between `TCGA_DESC` and `Cancer_Type_(matching_TCGA_label)` remains unresolved and should be addressed before either is used as the definitive cancer-type label in modelling. Cancer-type-level heatmaps confirmed that drug sensitivity profiles are heterogeneous even within a single TCGA class, reinforcing the need for cell-line-level analyses. Similarly, drugs with the highest standard deviation in LN_IC50 are strong candidates for biomarker-driven stratification studies, as their variable responses suggest the existence of identifiable sensitive subpopulations.

## 10. Limitations and Future Work

### 10.1 Limitations

While the dataset is large and well-structured, several limitations are worth acknowledging. Averaging repeated measurements to resolve duplicates, although a reasonable approach, may obscure screen-specific technical variation. Ideally, batch identity would be included as a covariate in any subsequent modelling.

All analyses in this phase are descriptive in nature. No formal statistical inference has been applied, meaning the differences and patterns observed in plots have not yet been tested for significance or effect size. The incomplete drug–cell-line matrix introduces a potential source of bias if the missing combinations are non-randomly distributed across cancer types or drug classes. The binary genomic flags, as discussed, reflect data availability rather than biological measurements, and their use as stratifying variables requires careful interpretation. Finally, the discrepancy between the two TCGA label columns has not been fully investigated and could affect any cancer-type-level summaries if the wrong column is selected.

### 10.2 Recommended Next Steps

Beyond statistical testing, integrating actual molecular profiles — including mutation burden scores, continuous CNA values, and expression levels — rather than the binary availability flags would substantially strengthen the biological interpretability of any model. Dimensionality reduction techniques such as PCA or UMAP applied to the drug sensitivity matrix could help identify clusters of co-sensitive cell lines and provide a more intuitive map of the pharmacological landscape. Pathway enrichment analysis on the most potent and most variable drug targets would add biological context to the pharmacological rankings. The `TCGA_DESC` discrepancy should be resolved by consulting the GDSC data dictionary. Finally, building predictive models such as elastic-net regression or random forests to predict LN_IC50 from genomic features is the stated primary goal of Stage Two and will build directly on the foundations laid here.

## 11. Appendix: Variable Glossary

| Term | Definition |
|---|---|
| **IC50** | Half-maximal Inhibitory Concentration: the drug concentration that reduces cell viability by 50% |
| **LN_IC50** | Natural logarithm of IC50; lower values indicate greater potency |
| **AUC** | Area Under the dose-response Curve; lower AUC indicates greater drug effect across all concentrations |
| **Z_SCORE** | Standardised LN_IC50 score adjusted for screen-to-screen variation |
| **COSMIC_ID** | Identifier in the Catalogue of Somatic Mutations in Cancer |
| **TCGA** | The Cancer Genome Atlas, a classification system for cancer types |
| **MSI** | Microsatellite Instability, a marker of defective DNA mismatch repair; MSI-H tumours tend to respond differently to immunotherapy |
| **CNA** | Copy Number Alteration, referring to genomic gains or losses of chromosomal segments |
| **EDA** | Exploratory Data Analysis |
| **GDSC** | Genomics of Drug Sensitivity in Cancer project |

*Report prepared from analysis notebooks `data_exp.ipynb` and `notebook.ipynb` as part of AI for Genomics Internship, Stage One.*


