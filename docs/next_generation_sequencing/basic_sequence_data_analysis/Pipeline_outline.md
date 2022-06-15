---
sidebar_position: 0
---

# An outline of the pipeline

In this tutorial we'll implement a basic pipeline for analysing short-read sequence data. And we
will use the results to look for variation in the sample compared to the reference genome assembly.

To make this as simple as possible to run, we've organised this tutorial into several steps with
commands detailed below. Then for each step we have a seperate page that explains more about how to
interpret the results.

**Note.** Before you start, you should have obtained the sequence reads (FASTQ files) and reference
sequence (FASTA file) as described on the [Prerequisiutes page](Prerequisistes.md). Most commands
will be run in a Terminal window, and you should ensure you are in the correct directory by typing:
```
cd ~/sequence_data_tutorial
````

first.  (The command `cwd` can be used to see where you are in the filesystem.)

## The pipeline in brief

The pipeline goes like this:

0. First [have a look at the FASTQ files](Inspecting_the_fastqs.md).  *How many read pairs are in the file?  What is the read length?*

1. [Perform quality control (QC) on the sequence reads](Quality_control.md).  *What is the GC content in the reads?  What is the duplication rate?*

**Note.** This step reads in FASTQ files (sequence reads) and outputs an HTML report with lots of details of
the sequencing. It can be run like this:

```
mkdir fastqc_output
fastqc -o fastqc_output ERR*.fastq.gz
```

2. [Align the reads](Aligning_reads.md).  *How are the reads represented in that file?  How many were aligned?  How many were not mapped?*

**Note.** This step reads in the FASTQ files and the genome assembly FASTA file. It outputs a list
of read alignents in SAM format.  It can be run in the terminal by running `bwa mem`:
```
bwa mem -o ERR377582-aligned.sam Pf3D7_v3.fa.gz ERR377582_1.fastq.gz ERR377582_2.fastq.gz
```

This step will take 10 minutes to run. When you have looked at the output, consider reading about
[short read theory](Short_read_theory.md) or [estimate the sequencing error
rate](De_novo_error_rate_estimation.md).

3. [Postprocess the alignments](Aligning_reads.md) to make them analysis-ready. *What steps are taken here?  How does this change the alignments?*

This step can be run like this:
```
mkdir tmp
samtools fixmate -m ERR377582-aligned.sam ERR377582-fixmate.bam
samtools sort -T tmp -o ERR377582-sorted.bam ERR377582-fixmate.bam
samtools markdup ERR377582-sorted.bam ERR377582-final.bam

# rename output file and index
mv ERR377582-final.bam ERR377582.bam
samtools index ERR377582.bam
```

4. [Inspect alignment statistics](alignment_statistics.md). 

This can be run like this:
```
samtools stats ERR377582.bam
```

5. [

## Challenges

So you 
