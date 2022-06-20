---
sidebar_position: 5
---

# Phasing the variant calls

If you got this far you should have a robust, well-QCd set of variants. But we're not done! Humans
are diploid and so we'd like to know how the variants stack up on each chromosome - that is how
they form **haplotypes**.

We don't directly observe haplotypes but they can be inferred through a process known as
**phasing**. This is done using a statistical algorithm that iteratively updates haplotype phase
(starting at random) by modelling each individual's haplotypes as a mosaic of segments of other
haplotypes in the data. (Even though this starts with random phase, people with homozygous
genotypes have stretches of genuine haplotype so there is some information there to start with).
'Switches' between haplotypes generally occur at recombination hotspots, so this process also uses
an estimate of the underlying recombination rate. This underlying model is known as the Li and
Stephens model after [the famous paper in which it was
introduced](https://pubmed.ncbi.nlm.nih.gov/14704198/), although it has seen many updates in
practice. Recent versions for large samples are based on the [Positional Burroughs-Wheeler
Transform](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3998136/) to quickly match haplotypes.

We'll use a recent implementation called EAGLE to do this.

One limitation is that EAGLE expectes variants to be bi-allelic. To make this true we will
represent all variants as multi-allelic in the VCF file. (This doesn't change the underlying calls,
but we end up with multiple rows for each multi-allelic variant):

```
bcftools norm -m -Oz -o GWD_30x_calls.filtered_split_.vcf.gz GWD_30x_calls.filtered.vcf.gz
```

**Question.** Compare these files (using `bcftools view [file] - less -S`.)  Can you see what this has done?

Now we are ready to phase.  The command is:

```
eagle \
  --geneticMapFile genetic_map/genetic_map_hg38_withX.txt.gz
  --vcf GWD_30x_calls.filtered_split_.vcf.gz \
  --numThreads 2 \
  --outPrefix GWD_30x_calls.phased
```

**Note.** This will take a few minutes to run. As always, watch the screen output for a wealth of
information about what is happening.


### What has phasing done?


