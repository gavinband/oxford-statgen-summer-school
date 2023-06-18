---
sidebar_position: 2
---

# Practical outline

### Overview 

In this tutorial we will demonstrate a basic pipeline for analysing paired-end short-read genomic
sequencing data. We will start with raw data in a FASTQ file, inspect quality control metrics,
align the data, and then use it to look for genetic variation.

#### Prerequisites

If you got here you should hopefully have already downloaded the practical data - **if not**, please
follow the instructions for doing that on the [prerequisites page](Prerequisites.md), and then come
back here.

#### A look at the data

You should now have a folder called `sequence_data_analysis` filled with a number of data
files. Now would be a good point to explore what's in there.  The folder contains

- sequence data reads from malaria parasites (under the `malaria/` folder). These reads are in
  *gzipped fastq* format.

- A malaria reference genome assembly (`Pf3D7_v3.fa`) in *FASTA* format.

- A similar set of sequence data reads and reference genome for a human sample - in the `human/` folder.

**Note.** I placed online a
[set of solutions files for steps in the practical](https://www.well.ox.ac.uk/~gav/projects/oxford_statgen_summer_school/sequence_data_analysis/solutions).
Feel free to check these as you go along.

During the practical we'll bring one or more of these datasets to an analysis-ready state.

To get started, start a terminal window and change directory into that folder:
```
cd sequence_data_analysis
```

### The practical in a nutshell

This practical works as follows: for each step there's a page giving you some information about how
to run the step; the page then links back to this one so you can see the next step. 

**Note.** Most of these example work with the malaria data in:
```
malaria/QG0033-C_Illumina-HiSeq_read1.fastq.gz
malaria/QG0033-C_Illumina-HiSeq_read2.fastq.gz
```

If you want to go off-piste, feel free to work with any (or all) of the others (that's one of the
[Challenge questions](Challenge_questions.md)). Just remember that the human data will align to
human genome and the malaria data to the malaria genome.

#### Now you are ready to start - go!

The steps are as follows. *Can you work out the answers to the questions below?*

1. First **[have a look at the FASTQ files](Inspecting_the_fastqs.md)**.  

Questions: *How many read pairs are in the file?  What is the read length?*

2. **[Perform quality control (QC) on the sequence reads](Quality_control.md)**.  

Questions: *What is the GC content in the reads? What is the fragment duplication rate? Are there
any sequencing artifacts?*

3. **[Align the reads](Aligning_reads.md)**.

Questions: *How are the reads represented in the aligned output file?  How many reads were aligned?  How many were not nmapped?*

4. **[Inspect read pileups and looking for variation](Viewing_alignments.md)**.

Questions: *Can you find a SNP?  An insertion or deletion?  A structural variant?*

Feel free to work through at your own pace, and to ask an instructor or your neighbour for help if you get stuck or
have a question. We plan to spend approximately one hour on this part, so spending up to ~30mins on steps 1 & 2, and
~30mins on steps 3 & 4 would be fine. 

:::tip Note
This might not be enough time to explore all the material fully so don’t worry if
you don’t finish in the time today – you can always re-visit these pages later.
:::

## Challenge questions

If you get this far, congratulations - you're an expert!

To test your mettle, here are some [challenge questions](Challenge_questions.md).  Good luck!

