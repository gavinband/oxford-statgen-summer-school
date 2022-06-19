---
sidebar_position: 1
---
# Getting the prerequisites

To run this tutorial you will need a few things.  This page will get you set up.

**Short version** This page basically tells you to:

* Log into [statgenschool.bmrc.ox.ac.uk](https://statgenschool.bmrc.ox.ac.uk) and select the Monday image.
* Download and extract [the data](https://www.well.ox.ac.uk/~gav/projects/oxford_statgen_summer_school/day_one_morning.tgz) into your jupyterhub image using `wget` and `tar`.  (The commands are detailed below).

Once you have done that feel free to **go ahead [to the practical](Pipeline_outline.md)**.

#### Using  JupyterHub

The Oxford Statistical Genomics Summer School practicals are being run on JupyterHub hosted at
[statgenschool.bmrc.ox.ac.uk](https://statgenschool.bmrc.ox.ac.uk). You should have already logged
in to this site - please make sure you are working in the 'Monday' image for this work. (If you're not already in there, see
[how to access and switching your JupyterHub image](../../switching_images.md) and then come back here.)

For this practical we'll be using both a terminal window and an R session.

#### Getting the data

To get the data, open a terminal and run the following commands:

Download the data: 
```sh
wget https://www.well.ox.ac.uk/~gav/projects/oxford_statgen_summer_school/day_one_morning.tgz
````

Extract it :
```
tar -xzf day_one_morning.tgz
```

You should now have a folder called `sequence_data_analysis`. Let's first delete the tarball and
then change into that directory:

```
rm day_one_morning.tgz
cd sequence_data_analysis
```

#### Getting the software

If you are running this as part of the Oxford Statistical Genomics Summer School, most of the
software should already be installed for you (remember to use the `Monday` image.)

The exception is the [IGV desktop app](https://igv.org) which you will need to download and install
on your laptop (**not** on jupyterhub). That's only needed for the [challenge
question](Challenge_questions.md) though, so don't worry if you don't have this.

When you're ready, go on to the [practical](Pipeline_outline.md).

#### Checking the software (optional)

We will be using `jellyfish`, `fastqc`, `samtools`, and `bwa` for this analysis. (A challenge
question also uses the [IGV desktop app](https://igv.org) which you will need installed on your
laptop if you choose that challenge.)

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

