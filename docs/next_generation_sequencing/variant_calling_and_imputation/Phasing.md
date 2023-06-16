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

**Question.** Compare these files (using `bcftools view [file] | less -S`.)  Can you see what this has done?

Now we are ready to phase.  The command is:

```
eagle \
  --geneticMapFile genetic_map/genetic_map_hg38_withX.txt.gz
  --vcf GWD_30x_calls.filtered_split_.vcf.gz \
  --numThreads 2 \
  --outPrefix GWD_30x_calls.phased
```

**NOte.** The backslash (`\`) is a *line continuation character*.  The above is run as a single command.

**Note.** This will take a few minutes to run. As always, watch the screen output for a wealth of
information about what is happening.

### What has phasing done?

Let's take a look at the first few variants and the first three samples.  Before phasing:
```sh
bcftools view -H GWD_30x_calls.filtered.vcf.gz | head -20 | cut -f 10-12
```

And after:
```sh
bcftools view -H GWD_30x_calls.phased.vcf.gz | head -20 | cut -f 10-12
```

You should see something like this:
```
1/1:255,111,0	1/1:255,108,0	1/1:255,105,0
0/0:0,99,255	0/0:0,102,255	0/0:0,111,255
0/0:0,108,255	0/0:0,117,255	0/0:0,102,255
0/0:0,102,255	0/0:0,99,255	0/0:0,96,255
0/0:0,114,255	0/0:0,87,255	0/0:0,111,255
0/0:0,120,255	0/0:0,99,255	0/0:0,93,255
0/0:0,120,255	0/1:181,0,250	0/0:0,117,255
0/0:0,75,255	0/0:0,123,255	0/0:0,120,255
0/0:0,93,255	0/0:0,108,255	0/0:0,141,255
0/0:0,84,255	0/0:0,93,255	0/0:0,102,255
0/0:0,114,255	0/1:219,0,239	0/0:0,111,255
0/1:254,0,203	0/0:0,154,255	0/0:0,111,255
0/0:0,108,255	0/0:0,151,255	0/0:0,111,255
1/1:255,81,0	0/1:243,0,183	1/1:255,102,0
1/1:255,99,0	0/0:0,138,255	0/1:255,0,221
0/0:0,105,255	0/0:0,141,255	0/0:0,141,255
0/0:0,96,255	0/0:0,111,255	0/0:0,105,255
0/0:0,87,255	0/0:0,111,255	0/0:0,111,255
0/0:0,84,255	0/1:255,0,187	0/0:0,114,255
1/1:255,99,0	0/1:158,0,255	1/1:255,126,0
```

and
```
1|1:255,111,0	1|1:255,108,0	1|1:255,105,0
0|0:0,99,255	0|0:0,102,255	0|0:0,111,255
0|0:0,108,255	0|0:0,117,255	0|0:0,102,255
0|0:0,102,255	0|0:0,99,255	0|0:0,96,255
0|0:0,114,255	0|0:0,87,255	0|0:0,111,255
0|0:0,120,255	0|0:0,99,255	0|0:0,93,255
0|0:0,120,255	1|0:181,0,250	0|0:0,117,255
0|0:0,75,255	0|0:0,123,255	0|0:0,120,255
0|0:0,93,255	0|0:0,108,255	0|0:0,141,255
0|0:0,84,255	0|0:0,93,255	0|0:0,102,255
0|0:0,114,255	1|0:219,0,239	0|0:0,111,255
1|0:254,0,203	0|0:0,154,255	0|0:0,111,255
0|0:0,108,255	0|0:0,151,255	0|0:0,111,255
1|1:255,81,0	0|1:243,0,183	1|1:255,102,0
1|1:255,99,0	0|0:0,138,255	1|0:255,0,221
0|0:0,105,255	0|0:0,141,255	0|0:0,141,255
0|0:0,96,255	0|0:0,111,255	0|0:0,105,255
0|0:0,87,255	0|0:0,111,255	0|0:0,111,255
0|0:0,84,255	1|0:255,0,187	0|0:0,114,255
1|1:255,99,0	0|1:158,0,255	1|1:255,126,
```

Phasing has **not changed the genotype likelihoods** (the numbers, in the PL field). But it has
turned forward slashes (`/`) into vertical bars (`|`)!  Great!

If you look more closely though, you'll see that Eagle has done something interesting to the
heterozygous variants. For example, look at sample 2. It has several heterozygous genotypes that
look like `0/1`. But phasing has turned these either into `0|1` or `1|0` - it has *phased* them
relative to each other.

You can interpret this as follows: for each sample (column), the number before the `|` reflects the
first haplotype in the sample (it might for example be the maternally-inherited chromosome,
although we don't know that from this file). The number after the `|` reflects the second haplotype
in the sample (this might be the paternally inherited chromosome, say). And before a `0` means the
haplotype carries the reference allele and a `1` means the non-reference allele.

### Next steps

When you're read [go back to the practical](README.md#steps-in-the-practical) and go on to **imputation**.

