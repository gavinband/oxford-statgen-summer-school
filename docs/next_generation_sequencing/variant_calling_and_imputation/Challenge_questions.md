---
sidebar_position: 10
---

# Variant calling challenge questions

## Hidden dangers of filtering

Our variant filtering approach actually had a unintended **hidden side-effect** (also known as a serious bug.)
Can you figure out what it is?

Hint: look at different types of variant in the `info` dataframe in your python code.  Or look at the output file:
```
bcftools view -H GWD_30x_calls.filtered.vcf.gz | less -S
```

One type of variant just isn't there isn't that output file!  Can you figure out why?

:::tip 

This illustrates the need, in real work, to carefully inspect and inspect outputs at all stages - it's very easy for unexpected things to go wrong.

:::

## Assessing the transition / transversion ratio.

We learned in the lectures that Ts/Tv can be used to assess the quality of a set of variant calls,
and that this number should be 2 or greater. **What is the Ts/Tv for this data before and after
filtering?**

**Note.** The `bcftools stats` command can calculate Ts/TV (as TSTV) statistics - search for 'TSTV' in the output.

## Inspecting failed variants

A useful test of the filtering strategy is to pick at some variants that failed the filters, and look at what the pileups look like.
This might make it clear what (if anything) is wrong with the variant.

You could of course find these variants from your python code, but here's an alternative approach using `bcftools`.

We'll first generate a new version of the filtered VCF that keeps the failed variants in but marks them as 'FAIL' (as the original version removed them entirely). You can produce
it by adding the `-s` option to the `bcftools filter` command, to 'soft filter' the file:

```
bcftools filter -Oz \
-i 'DP>=1000 & DP<=10000 & MQ>=40 & MQ0F<=0.05 & AN>=300 & AC>=2 && QUAL>=50 & RPB>=0.0001 & MQB>=0.0001 & BQB>=0.0001 & MQSB>=0.0001' \
-o 'GWD_30x_calls.soft-filtered.vcf.gz' \
-s FAIL \
calls/GWD_30x_calls.vcf.gz
```

:::tip Note
Use the same filters as you used in the [Variant quality control](./Variant_quality_control.md) section, if you changed them.
:::

Let's also *index* the file, which helps `bcftools` select specific variants or regions:
```
bcftools index GWD_30x_calls.soft-filtered.vcf.gz
```

Now the SNPs in the data can be viewed like this:
```
bcftools view \
-H \
-r chr19:48693822-48708100 \
--types snps \
GWD_30x_calls.soft-filtered.vcf.gz \
| less -S
```

:::tip Note
We used the `-r` option above to restrict to the region that our subsetted read data in `reads/` corresponds to.
:::

Scroll around in this file a bit and pick a SNP that says `FAIL` in the `FILTER` column - for example, in my filtering, there's a SNP at chr19:48700152 that failed the filter.

We can now use `bcftools query` to pull out a list of samples that appeared to carry the alternate alleles at this variant. The syntax is a bit arcane, but very powerful - it's
described on the [bcftools manual page](https://samtools.github.io/bcftools/bcftools.html#query).  For our purposes it can be done like this:

```
bcftools query \
--regions chr19:48700152-48700152 \
--include 'TYPE="snp" & GT="alt"' \
--format 'Variant:\n%CHROM:%POS %REF>%ALT %FILTER\nINFO=%INFO\n\nNon-ref samples:\n[%SAMPLE %GT\n]' \
GWD_30x_calls.soft-filtered.vcf.gz \
| head -n 30
```

You should see the INFO fields for the variant, followed by a list of the first few samples that were called with non-reference genotypes.

Can you see why it failed?

You can also `samtools tview`, as in the [earlier practical](../basic_sequence_data_analysis/Viewing_alignments.md) to inspect the pileup for these samples:

```
samtools tview --reference GRCh38_chr19.fa.gz -p chr19:48700152 reads/HG02588
```

Better yet - load up the data in IGV and see if you agree with the filtering.
