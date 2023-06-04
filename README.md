
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RNAProbeBuilder

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/RNAProbeBuilder)](https://CRAN.R-project.org/package=RNAProbeBuilder)
[![R-CMD-check](https://github.com/selkamand/RNAProbeBuilder/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/RNAProbeBuilder/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

RNAProbeBuilder lets you input mutations in mRNA space and produces
mutant and wildtype mRNA sequences that could be used to probe the
sequences in supported experimental systems

**WARNING**: This package is in *very early development* and not ready
for use yet!!!

## Installation

You can install the development version of RNAProbeBuilder from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("selkamand/RNAProbeBuilder")
```

## Quick Start

Your input should be mutations in protein coding mRNA (mRNA HGVS
notation). For example

NM_004006.2:c.4375C\>T

(transcript:c.positionRef\>Alt)

``` r
mRNA_mutations <- c("NM_004006.2:c.4375C>T")

mRNA_mutations |> 
  probes_construct( # Constructs 1 probe for each mrna mutation
    bp_upstream = 20, 
    bp_downstream = 20,
    probe_type = "mRNA_no_U" # build cDNA probes setting to 'cDNA'
  ) |> 
  probes_collapse_duplicates() |> # Collapse identical probes (e.g. because 2 isoforms are different, they should be collapsed 
  probes_write_fastq('probes.fasta')
  
```

## FAQ

**Q: My mutations are in genomic coordinates. How do I get mutations in
mRNA HGVS notation**

A: Run VEP (see working with genomic inputs vignette)

**Q: My mutation may affect multiple isoforms? How do i handle this?**

A: The `probes_prepare_input_from_vep()` function will, by default,
produce a vector with the impact for EVERY isoform in the annotation
used by VEP. So when you pipe this into `probes_construct()` you’ll get
one probe per isoform hit, NOT 1 probe per genomic mutation.

Where probes for multiple mutant isoforms would have identical
sequences, `probes_collapse_duplicates()` will collapse these into
single probes. If you know in advanced what isoforms you care about -
simply filter the output of probes_prepare_input_from_vep() to include
only those isoforms

**Q: What variant types are supported**

A: Currently just SNVs and Indels. No fusions or large deletions
supported
