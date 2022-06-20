---
sidebar_position: 2
---

# Variant calling, phasing and imputation

This morning we focussed on quality control, aligning, and inspecting some sequence data for a single sample.

In this afternoon's session we're going to focus on using a set of sequenced samples to do a few things:

1. identify genetic variants in a region (and work out the sample genotypes at those variants)
2. *phase* the genotypes so we can see how they align along haplotypes
3. and use this to impute variants from a second dataset (for which we only have microarray data).

Not every possible genetic variant is a real one and to make this work well we will need to do some
more **quality control**. This time we'll take care to filter the set of variants based on sensible
metrics before we use the data.

All the data in this practical comes from the [IGSR](https://www.internationalgenome.org), which
lists a huge number of open-access datasets that you can use in your analysis (including the 1000
Genomes Project data).

**Note.** Before doing anything else, please make sure you have [downloaded the data](Prerequisites.md).
Then come back here.

## Steps in the practical

The practical has four main steps:

1. The first step is to [generate a VCF file of variant calls](Variant_calling.md) from some 1000 Genomes Project samples.

2. We'll then perform [quality control on the initial variant calls](Variant_quality_control.md).

3. Next we will [phase the calls](phasing.md) so we know how they stack up on haplotypes.

4. And then we'll [use these haplotypes to impute some microarray data](imputation.md)

Finally we have - you guessed it - some [challenge question](Challenge_questions.md).  Good luck!

