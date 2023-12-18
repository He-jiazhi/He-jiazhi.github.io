---
layout: page
title: Analysis of      HIV Drug Resistance 
description: STAT 2131 Course Final Project at Pitt
img: assets/img/apply1_1.png
importance: 1
category: application
---

The scientific goal is to determine which mutations of the Human Immunodeficiency Virus Type 1 (HIV-1) are associated with drug resistance. The data set, publicly available from the Stanford HIV Drug Resistance Database [Stanford HIV Drug Resistance Database](http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/), was originally analyzed in (Rhee et al. 2006).

## Preparing the data

The data set consists of measurements for three classes of drugs: protease inhibitors (PIs), nucleoside reverse transcriptase (RT) inhibitors (NRTIs), and nonnucleoside RT inhibitors (NNRTIs). Protease and reverse transcriptase are two enzymes in HIV-1 that are crucial to the function of the virus. This data set seeks associations between mutations in the HIV-1 protease and drug resistance to different PI type drugs, and between mutations in the HIV-1 reverse transcriptase and drug resistance to different NRTI and NNRTI type drugs (The raw data are saved as `gene_df`).

In order to evaluate our results, we compare with the treatment-selected mutation panels created by (Rhee et al. 2005), which can be viewed as the ground true. These panels give lists of HIV-1 mutations appearing more frequently in patients who have previously been treated with PI, NRTI, or NNRTI type drugs, than in patients with no previous exposure to that drug type. Increased frequency of a mutation among patients treated with a certain drug type implies that the mutation confers resistance to that drug type (The raw data are saved as `tsm_df`).

To simplify the analysis, in this project we will confine our attention to the PI drugs.

```{r}
drug_class = 'PI' # Possible drug types are 'PI', 'NRTI', and 'NNRTI'
```

## Fetching and cleaning the data

First, we download the data and read it into data frames.

```{r}
base_url = 'http://hivdb.stanford.edu/_wrapper/pages/published_analysis/genophenoPNAS2006'
gene_url = paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep='/')
tsm_url = paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep='/')

gene_df = read.delim(gene_url, na.string = c('NA', ''), stringsAsFactors = FALSE)
tsm_df = read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
names(tsm_df) = c('Position', 'Mutations')
```

A small sample of the data is shown below.

<table>
<thead>
<tr class="header">
<th align="right">APV</th>
<th align="right">ATV</th>
<th align="right">IDV</th>
<th align="right">LPV</th>
<th align="right">NFV</th>
<th align="right">RTV</th>
<th align="right">SQV</th>
<th align="left">P1</th>
<th align="left">P2</th>
<th align="left">P3</th>
<th align="left">P4</th>
<th align="left">P5</th>
<th align="left">P6</th>
<th align="left">P7</th>
<th align="left">P8</th>
<th align="left">P9</th>
<th align="left">P10</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="right">2.3</td>
<td align="right">NA</td>
<td align="right">32.7</td>
<td align="right">NA</td>
<td align="right">23.4</td>
<td align="right">51.6</td>
<td align="right">37.8</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">I</td>
</tr>
<tr class="even">
<td align="right">76.0</td>
<td align="right">NA</td>
<td align="right">131.0</td>
<td align="right">200.0</td>
<td align="right">50.0</td>
<td align="right">200.0</td>
<td align="right">156.0</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">F</td>
</tr>
<tr class="odd">
<td align="right">2.8</td>
<td align="right">NA</td>
<td align="right">12.0</td>
<td align="right">NA</td>
<td align="right">100.0</td>
<td align="right">41.0</td>
<td align="right">145.6</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
</tr>
<tr class="even">
<td align="right">6.5</td>
<td align="right">9.2</td>
<td align="right">2.1</td>
<td align="right">5.3</td>
<td align="right">5.0</td>
<td align="right">36.0</td>
<td align="right">13.0</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">I</td>
</tr>
<tr class="odd">
<td align="right">8.3</td>
<td align="right">NA</td>
<td align="right">100.0</td>
<td align="right">NA</td>
<td align="right">161.1</td>
<td align="right">170.2</td>
<td align="right">100.0</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">I</td>
</tr>
<tr class="even">
<td align="right">82.0</td>
<td align="right">75.0</td>
<td align="right">400.0</td>
<td align="right">400.0</td>
<td align="right">91.0</td>
<td align="right">400.0</td>
<td align="right">400.0</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">-</td>
<td align="left">I</td>
</tr>
</tbody>
</table>