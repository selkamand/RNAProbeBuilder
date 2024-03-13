
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RNAProbeBuilder

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/RNAProbeBuilder)](https://CRAN.R-project.org/package=RNAProbeBuilder)
[![R-CMD-check](https://github.com/selkamand/RNAProbeBuilder/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/RNAProbeBuilder/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

RNAProbeBuilder constructs probes for targeting RNA mutations. At least
two probes are generated per mutation, one targeted to the mutant, one
targeting wildtype sequences. This workflow will show you how to convert
**a VCF of genomic mutations** -\> **mRNA probes** (e.g. for use with
Visium platform)

Probes are all built in mRNA space, and by default with ‘U’s replaced
with ’T’. You can configure the output based on the \`probe_type\`
argument).

**WARNING**: This package is in *very early development* and not ready
for use yet!!!

## Installation

You can install the development version of RNAProbeBuilder from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("selkamand/RNAProbeBuilder")
```

## Usage

### **1 \| Prepare Inputs:**

RNAProbeBuilder probe designs are based on cDNA HGVS notation. Start by
taking your list of mutations, and annotate it using
[VEP](https://asia.ensembl.org/Tools/VEP) annotate impacts in cDNA space
(importantly this retrieves the impact across all isoforms of a
transcript)

Use the online / offline VEP (latest
[hg38](https://asia.ensembl.org/Tools/VEP) or
[hg37](https://grch37.ensembl.org/Homo_sapiens/Tools/VEP) cache are a
good choice)

You annotations should include cDNA HGVS notation like the following:
NM_004006.2:c.4375C\>T

(transcript:c.positionRef\>Alt)

**Configuring VEP**

![](inst/figs/VepConfig.png)

- Identifiers:
  - Gene Symbol
  - Transcript Version
  - HGVS
- Transcript Database
  - Ensembl<u>/GENCODE</u> transcripts

**Downloading Results**

Download / save all results in TXT format

![](inst/figs/DownloadVepResults.png)

This downloaded text file will serve as the input to

### 2 \| Connect to biomaRt

This will let us pull the mRNA sequences we need to construct probes

``` r
# If varants are GRCh37 / hg19:
ensembl <- load_biomart(GRCh = "37")

# If variants are GRCh38 / hg38
ensembl <- load_biomart(GRCh = "38")
```

### 3 \| Build probes from VCF

``` r
probes_read_vep_txt("path_to_vep_text") |> # Read in Probe data
  probes_construct( # Constructs 1 probe for each mrna mutation
    ensembl = ensembl,
    probe_size = 51,
    probe_type = "mRNA_no_U" # only mRNA supported for now
  )  |>
  probes_collapse_duplicates() |> # Collapse duplicate probes
  probes_write_output(outdir = "outdir", prefix = "myprobes") # Write FASTA / QC files
  
```

### 4 \| Build probes from fusion sequences

Currently we don’t support end-end fusion breakpoint -\> mRNA sequence
creation. However, we provide some utilities to simplify the process.
Start with the following input:

- Fusion_Sequence (mRNA nucleotide sequence with breakpoint indicated by
  ‘\|’)

``` r
probes_construct_fusion_sequence("ACTGACCGAC|TTTCTCTCTACATC", probe_size = 10)
```

## FAQ

**Q: My mutation may affect multiple isoforms? How do i handle this?**

A: By default, VEP describes the impact of mutations for EVERY isoform
in the annotations Ensembl/Gencode. `probes_read_vep_txt()` parses
isoform level impacts, and `probes_construct()` creates one WT/Mut pair
of probes per isoform hit, NOT 1 probe per genomic mutation.

Where probes for multiple mutant isoforms would have identical
sequences, `probes_collapse_duplicates()` will collapse these into
single probe pairs. If you know in advanced what isoforms you care
about - simply filter the output of `probes_read_vep_txt()` to include
only those isoforms

**Q: What variant types are supported**

A: Currently just SNVs and Dups. No Indels, fusions or large deletions
supported

**Q: What type of sequences are output? mRNA? cDNA? etc?**

A: By default 5’ -\> 3’ mRNA sequences with Uracils replaced by Thymines
are supplied (mRNA_no_U)

![](inst/figs/Sequence_Output_Type)

**Q: Why don’t you support cDNA output?**

A: cDNA can have different definitions:

1.  A synthetic sequence transcribed from an mRNA molecule, this ‘cDNA’
    is complementary to the mRNA
2.  The mRNA molecule itself with Uracils replaced by Thymines
    (identical to **mRNA_no_U**)

Due to this ambiguity we don’t support an ‘explicit’ cDNA output. If you
want ‘cDNA’, we advise you to convert directly from **mRNA_no_U** based
on which type of ‘cDNA’ you’re after.
