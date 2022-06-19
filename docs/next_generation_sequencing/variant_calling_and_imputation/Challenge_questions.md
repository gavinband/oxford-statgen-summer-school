---
sidebar_position: 10
---

# Challenge questions

## Questions on variant QC.

### Viewing variants in `tview`

Pick some variants that failed and visualise the BAM files with `samtools tview`. (Pick samples
that have the alternative alleles.)

**Note.** You will have to take variants in the range `chr19:48693971-48707951` to be able to see them in
the read files supplied in `reads`.

### Assessing the transition / transversion ratio.

We learned in the lectures that Ts/Tv can be used to assess the quality of a set of variant calls,
and that this number should be 2 or greater. **What is the Ts/Tv for this data before and after
filtering?**

**Note.** The `bcftools stats` command can calculate Ts/TV (as TSTV) statistics.

## Post-imputation questions.

### Imputation $r^2$

The R2 field is a measure of how well the imputation software THINKS it is doing - it is an interal
estimate of what it thinks the correlation between the imputed and true genotype should be.

**Question.** What is the distribution of R2 across all variants? How does it vary with MAF?

### Imputation ER2

The ER2 field measures the actual "leave-one-out" R2 for this variant (i.e. the correlation between
what the model would have predicted for this variant and what it actually observed).

**Question.** How do R2 and ER2 compare? Is the software's guess at its accuracy (R2) close to the
emperical accuracy (ER2) well?  How does this vary with MAF?

### Calling secretor status

The SNP rs601338 at `chr19:48703417` in *FUT2* determines whether an individual secrets the Lewis
antigen. An individual who is homozygous AA doesn't secret Lewis antigen ("secretors"), and
everything else does ("non-secretors"). We have provided a file
(`secretor/GGVP_secretor_status.txt`) that includes the real secretor status for these individuals.

**Question**. How accurate is the secretor status called from imputation, in comparison to the real data provided?

