---
sidebar_position: 21
---

# Appendix 3: read trimming

Read trimming refers to the process of removing trailing parts of the reads before analysis for various reasons including

- the presence of adapter sequences in the reads (trim everything after the adapter)
- the presence of low-quality bases toward the end of the read (trim everything after the low-quality bases)
- the presence of other contaminant or artifactual sequences (for example poly-G sequences).

There are several programs that can do this, including
[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic),
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) and
[fastp](https://github.com/OpenGene/fastp). If you have lots of adapter contamination you will may
to explore using these.

However, there are some reasons not to worry about this:

- many datasets (such as the one we are working on here) have low levels of contamination by
  adapters or other sequences.
  
- base qualities are also generally high, and they are taken into account by many downstream tools
  such as variant callers

- most importantly, modern aligners allow for subregions of reads that do not align. (They will
  generally be 'clipped' off). Thus, these sections might not affect analysis much anyway.
  
For the purposes of this practical we will skip the read trimming step. However, if I were
analysing a dataset with high levels of contamination - like [this
one](https://www.well.ox.ac.uk/~gav/projects/oxford_statgen_summer_school/sequence_data_analysis/fastqc_examples/human/HV31-illumina_novaseq_2_fastqc.html) or worse - I would definitely be tempted to implement trimming.

