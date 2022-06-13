# Introduction

To run this tutorial you will need a few things:

## Start in a new, empty folder

We recommend first creating a new empty folder to run the practical in.  You could do that like so:
```
mkdir sequence_data_tutorial
cd sequence_data_tutorial
```

## Getting the read data

To run the tutorial you first need a pair of fastq files representing paired-end short-read sequencing.

To make the tutorial run in a reasonable time, this needs to be not too big, but to make it useful,
it needs to be big enough. I've prepared a set of files from sequencing of *Plasmodium falciparum*
which you can find [here](https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Next_Generation_Sequencing/practicals/ngs_processing_pipeline/data/reads/subsampled/).

Pick one of these samples and download its data - note you need both the `[something]_1.fastq.gz`
and `[something]_2.fastq.gz` files (because these contain, respectively, the first and second read
in each pair).

**Note** if you are working in a remote environment (such as a JupyterHub instance) then use the `wget` command to download it, for example:
```
wget https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Next_Generation_Sequencing/practicals/ngs_processing_pipeline/data/reads/subsampled/ERR377582_1.fastq.gz
wget https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Next_Generation_Sequencing/practicals/ngs_processing_pipeline/data/reads/subsampled/ERR377582_2.fastq.gz
```

Feel free to choose a different sample (or download more of them) as well.

## Getting the reference sequence

For *P.falciparum* malaria the standard reference is known as `Pf3D7_v3` (version 3 of the assembly
built from the 3D7 isolate. 3D7 was originally isolated from an individual who lived near Schipol
airport. This person had never left The Netherlands but nevertheless became sick with malaria,
presumably because an infected mosquito hitched a ride on an incoming flight. The parasite is
thought to be of west African origin.

You can download the reference sequence
[here](https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Next_Generation_Sequencing/prac
ticals/ngs_processing_pipeline/data/reference/).

**Note** if you are working in a remote environment (such as a JupyterHub instance) then use the
`wget` command in a terminal to download it:

```
wget https://www.well.ox.ac.uk/~gav/projects/gms/statistics-course/Next_Generation_Sequencing/practicals/ngs_processing_pipeline/data/reference/Pf3D7_v3.fa.gz
```

## Getting the software

We will be using `jellyfish`, `fastqc`, `samtools`, and `bwa` for this analysis. Later, we will also use
`IGV`.

To check that you have the right software, open a bash terminal and type the program names followed
by &lt;enter&gt;. (Don't type the dollar signs, these are there to show you the bash prompt). You
should see something like:

```bash
$ jellyfish --version
jellyfish 2.3.0

$ fastqc --version
FastQC v0.11.9

$ samtools --version
samtools 1.8
Using htslib 1.8
Copyright (C) 2018 Genome Research Ltd.

$ bwa --version

Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1188
```

That's it!  When you're ready, go to the page on [Inspecting the fastqs](Inspecting_the_fastqs.md).
