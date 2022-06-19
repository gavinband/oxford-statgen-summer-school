---
sidebar_position: 2
---

# Practical outline

### Overview 

In this tutorial we'll implement a basic pipeline for inspecting and analysing short-read sequence
data. We'll use the results to look for variation in the sample compared to the reference genome
assembly.

#### Getting started

You should hopefully have already downloaded the practical data - if not please follow the
instructions for doing that on the [prerequisites page](Prerequisites.md), and then come back here.

#### What's in the data?

You should hopefully have a folder called `sequence_data_analysis` filled with a number of data
files. Now would be a good point to explore what's in there.  The folder contains

- sequence data reads from malaria parasites (under the `malaria/` folder). These reads are in
  *gzipped fastq* format.

- A malaria reference genome assembly (`Pf3D7_v3.fa`) in *FASTA* format.

- A similar set of sequence data reads and reference genome for a human sample - in the `human/` folder.

It also contains a set of solutions files - these are the expected output from this practical. Feel
free to check these out as you go along.

The aim of this practical is to bring one or more of these datasets through a basic pipeline that
brings them to an analysis-ready state. We've include multiple datasets so you can pick the one you
want to work with and see some differences between them. (Or, if you're quick, you can process them
all.)

To get started, start a terminal window and change directory into that folder:
```
cd sequence_data_analysis
```

### The practical in a nutshell

This practical works as follows: for each step there's a page giving you some information about how
to run the step; the page then links back to this one so you can see the next step. 

**Note.** The examples will generally work with the malaria data - specifically the files:

* `malaria/QG0033-C_Illumina-HiSeq_read1.fastq.gz` and
* `malaria/QG0033-C_Illumina-HiSeq_read2.fastq.gz`

But if you want to go off-piste, feel free to work with any (or all) of the others. Just remember
that the human data will align to human genome and the malaria data to the malaria genome.

The steps are as follows. *Can you work out the answers to the questions below? Please make a note of the answers for the consolidation session.*

1. First **[have a look at the FASTQ files](Inspecting_the_fastqs.md)**.  

Questions: *How many read pairs are in the file?  What is the read length?*

2. **[Perform quality control (QC) on the sequence reads](Quality_control.md)**.  

Questions: *What is the GC content in the reads? What is the fragment duplication rate? Are there
any sequencing artifacts?*

3. **[Align the reads](Aligning_reads.md)**.

Questions: *How are the reads represented in the aligned output file?  How many reads were aligned?  How many were not nmapped?*

4. **[Inspect read pileups and looking for variation](Viewing_alignments.md)**.

Questions: *Can you find a SNP?  An insertion or deletion?  A structural variant?*


## Challenge questions

If you get this far, congratulations - you are now expert!

To test your mettle, here are some [challenge questions](Challenge_questions.md).  Good luck!

