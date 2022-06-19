---
sidebar_position: 1
---
# Getting the prerequisites

To run this tutorial you will need a few things.  This page will get you set up.

#### Using  JupyterHub

If you are running this as part of the Oxford Statistical Genomics Summer School, you'll be using
JupyterHub to run the practicals, hosted at
[statgenschool.bmrc.ox.ac.uk](https://statgenschool.bmrc.ox.ac.uk). You should hopefully have
already logged in to this site. Please make sure you are working in the 'Monday' image for this
work.  (See [how to switch your JupyterHub image](../switching_images.md).)

For this practical we'll be using both a terminal window and an R session.

#### Getting the data

To get the data, open a terminal and run the following commands:

```sh
# download the data
wget https://www.well.ox.ac.uk/~gav/projects/oxford_statgen_summer_school/day_one_morning.tgz
# extract it
tar -xzf day_one_morning.tgz
```

You should now have a folder called `sequence_data_analysis`.  

**Note.** for the practical we'll generally assume you are working in that folder. In a terminal
you do this by typing `cd sequence_data_analysis`; in R, you could do it by typing
`setwd("sequence_data_analysis")`.  

**Note.** If something's not working for you, please check which folder you are in first, by typing
`cwd` in the terminal or `getwd()` in R.

#### Getting the software

If you are running this as part of the Oxford Statistical Genomics Summer School, most of the
software should already be installed for you (remember to use the `Monday` image - )

We will be using `jellyfish`, `fastqc`, `samtools`, and `bwa` for this analysis. Later, we will also use
`IGV`.

To check that you have the right software, open a bash terminal and type the program names followed
by &lt;enter&gt;. (Don't type the dollar signs, these are there to show you the bash prompt ). You
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

That's it!  When you're ready, go back to the [practical](Pipeline_outline.md).
