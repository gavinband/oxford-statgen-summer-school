---
sidebar_position: 6
---

# Imputing a set of microarray data

Ok, so we've [called genetic variants](Variant_calling.md), we've [QCd
them](Variant_quality_controls.md), and we've [phased them into haplotypes](Phasing.md). What can
we do with that?

One of the important ways that genetic variation reference panels is used is to **impute**
genotypes in other sample sets that (for reasons of cost or practicality) haven't been sequenced.
This is the basic paradigm, for example, for most analyses of the
[genetic resource in the UK Biobank](https://www.nature.com/articles/s41586-018-0579-z).

Imputation will be covered in more detail in another session, but the basic idea is that - even if
we haven't genotyped a particular variant in our sample set - we have genotyped variants on the
same haplotype. So by statistically matching haplotypes between our genotyped set and the reference
panel we've created using sequencing we might be able to infer their genotypes, even at untyped
markers. (This process relies on the interesting structure that human haplotypes have, in which
variants along the genome become correlated due to genetic drift, but these patterns are broken
down by recombination.)

Imputation from very large reference panels can now be done easily online using for example the
[Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html)
or the [Sanger imputation service](https://imputation.sanger.ac.uk).  But here we want to demonstrate how it's done.

To this end we've placed a dataset from microarray-based genotyping of a second set of Gambian
samples in the file: `GGVP/omni2.5M/GGVP/omni2.5M/GGVP-illumina_omni2.5M.phased.vcf.gz`. (They come
from the
[Gambian Genome Variation Project](https://www.internationalgenome.org/gambian-genome-variation-project/)
which is why we've called them `GGVP`.)

**Note.** We have **already QCd and phased these genotypes** for you - in a real analysis you might
have to do that yourself.

Let's use our phased sequence data to see if we can impute a particular SNP - the *secretor status* SNP rs601338. This is a G>A mutation at `chr19:48,703,417` which determines
whether individuals secrete the 'Lewis antigen' - and hence A or B blood type antigens - into mucus. You can check that it is not genotype in the input data:

```
bcftools view \
-H \
-i 'POS=48703417' \
GGVP/omni2.5M/GGVP-illumina_omni2.5M.phased.vcf.gz \
| wc -l
```

**Note.** recall that `wc -l` counts the number of lines in its output.

To use imputation we will use the program `MINIMAC`.  The first step is to convert our reference panel to the MINIMAC format:

```
Minimac3 \
  --refHaps GWD_30x_calls.phased.vcf.gz \
  --processReference \
  --prefix GWD_30x_calls.phased \
  --chr "chr19"
```
This should output a transformed file `GWD_30x_calls.phased.m3vcf.gz`.

Now let's use MINIMAC to impute:
```
minimac4 \
  --refHaps "GWD_30x_calls.phased.m3vcf.gz" \
  --mapFile genetic_map/genetic_map_hg38_withX.txt.gz \
  --haps GGVP/omni2.5M/GGVP-illumina_omni2.5M.phased.vcf.gz \
  --ignoreDuplicates \
  --format GT \
  --prefix GGVP-illumina_omni2.5M.imputed
```

### What did imputation do?

What's happened there?  Well let's look at how many variants were in the files.
```
#Number of variants in the microarray data:
bcftools view -H GGVP/omni2.5M/GGVP-illumina_omni2.5M.phased.vcf.gz | wc -l

# Number of variants in the reference panel
bcftools view -H GWD_30x_calls.phased.vcf.gz | wc -l

# Number of variants in the imputed data
bcftools view -H GGVP-illumina_omni2.5M.imputed.dose.vcf.gz | wc -l
```

**Question.** What do these numbers mean?

### Inspecting the secretor status

Let's have a look at the secretor status SNP in the output:
```
bcftools view -H -i 'POS=48703417' GGVP-illumina_omni2.5M.imputed.dose.vcf.gz
```

You should see something like this:
```
chr19	48703417	chr19:48703417:G:A	G	A	.	PASS	AF=0.53395;MAF=0.46605;R2=0.99822;IMPUTED	GT	0|0	1|0	1|0	1|1	0|1	1|0	0|0	1|1	0|0	0|1	0|1	1|0	0|1	0|1	0|0	1|1	1|1	1|0	0|1	1|1	0|0	0|0	1|0	1|1	0|0	0|1	0|0	1|1	0|0	1|0	0|0	1|1	0|1	0|1	0|0	1|1	1|0	1|1	1|1	0|0	1|0	0|0	0|1	1|0	1|0	1|1	0|0	0|0	1|0	0|1	0|0	0|0	0|1	0|0	0|1	0|1	1|0	1|1	0|1	1|1	1|1	1|0	1|1	0|0	0|1	1|0	0|0	1|0	1|1	1|0	0|0	1|0	0|1	0|1	1|1	1|0	1|1	0|0	0|1	0|0	0|1	0|1	1|0	1|1	0|0	1|0	0|1	0|1	1|0	1|1	1|0	0|0	1|0	1|0	0|0	1|1	1|1	1|0	1|0	1|0	0|1	1|1	1|1	1|0	0|0	0|1	0|1	0|0	1|1	0|0	1|1	1|1	0|0	1|0	1|0	0|1	1|0	1|1	1|1	0|1	1|0	0|1	1|1	1|1	1|1	1|1	1|0	0|0	1|0	0|1	1|0	1|0	1|0	0|1	1|1	0|1	1|1	0|1	0|1	0|0	0|1	1|0	0|0	1|0	1|1	0|0	1|1	0|0	0|1	1|1	0|1	0|1	0|1	1|1	0|1	0|0	0|1	1|0	0|1	0|0	1|1	1|0	0|0	1|0	1|0	1|0	1|1	1|1	1|0	1|1	1|0	1|1	1|1	1|0	1|1	1|0	1|0	1|0	1|1	1|1	0|0	0|1	0|1	0|1	1|1	0|1	1|1	1|1	0|1	0|0	0|0	0|0	1|1	0|1	0|0	1|1	1|1	1|0	1|0	1|1	0|1	0|1	1|1	1|0	1|0	1|1	1|1	1|0	0|0	0|0	1|0	1|0	0|1	1|1	0|1	0|1	1|1	1|0	1|1	0|0	1|1	0|0	1|1	1|1	1|0	0|0	1|1	1|0	1|1	0|0	1|1	0|1	0|1	1|0	1|0	0|0	0|1	0|0	1|0	0|0	0|0	1|1	0|1	1|1	0|1	0|1	1|0	0|1	0|1	1|0	0|1	1|0	0|0	0|1	0|0	1|1	1|0	1|0	1|1	0|0	1|1	0|1	1|0	0|0	1|1	0|1	1|0	1|0	1|0	0|1	0|0	0|1	0|1	1|1	1|1	0|1	0|0	1|0	1|1	1|0	1|1	1|0	1|0	1|0	0|0	1|0	1|1	1|1	1|0	0|0	1|1	0|0	1|0	1|1	1|0	1|1	1|0	0|1	0|0	1|0	0|1	0|1	0|1	0|0	1|1	1|1	1|1	1|0	0|0	1|1	1|1	1|0	1|0	1|0	0|1	1|1	1|1	1|0	0|0	1|1	1|0	0|1	0|1	0|1	1|1	1|0	1|0	1|0	0|1	0|1	1|0	0|1	1|0	0|1	1|0	1|1	1|0	1|0	0|1	1|1	1|0	0|1	1|0	0|1	1|1	1|0	1|1	0|0	0|1	1|1	0|0	0|0	1|1	1|0	0|0	0|1	0|1	0|0	1|0	1|1	0|0	0|1	0|0	0|1	1|0	0|1	1|0	0|1	1|0	0|0	1|0	0|0	1|1	0|1	1|1	0|1	1|0	0|0	1|1	0|1	1|1	0|0	1|1	1|0	1|1	0|0	1|0	1|0	1|0	0|1	1|0	1|1	1|0	0|0	1|1	1|0	0|1	0|0	0|1	1|1	0|1	0|0	1|1	0|0	1|0	1|1	0|1
```

It worked!  Imputation has generated best-guess genotypes for `rs601338` (at `chr19:48703417`) for us.

:::tip Question
Inspect the VCF INFO fields for this SNP.  What frequency does it have?  What is the imputation R2 (a measure of predicted imputation accuracy)?
Use a similar `bcftools view` command to look up the SNP in our phased reference panel - is the frequency similar?
:::

### Next steps

You have successfully used a set of sequence data to identify genetic variants,
quality control and phase them. And you have used that to impute the important secretor status (and
other variants) into another dataset.  Congratulations!

To finish the practical, go back and try the [challenge questions](README.md#steps-in-the-practical).
