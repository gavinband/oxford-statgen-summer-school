# Performing quality control

## Running `fastqc`
There are a set of analyses that can be carried out directly on a FASTQ file to assess the quality
of sequencing. These include looking for read duplication, assessing GC content, looking for
over-represented sequences and the presence of sequence adapters in reads.

Fortunately there are good software packages to do this for us. The one we'll try is called
`fastqc` and it can be run as follows.  First make a directory for the output:
```
mkdir fastqc_output
```

```
fastqc -o fastqc_output ERR_*.fastq.gz
```

This will process both reads and will take a couple of minutes to run.

If you look in the folder now you'll see there are some new files:
```
$ ls fastq_output
ERR377582_1_fastqc.zip	ERR377582_2_fastqc.html
ERR377582_1_fastqc.html	ERR377582_2_fastqc.zip
...
```

`fastqc` has created some `.html` files, as well as some `.zip` files. Let's look by opening the
`.html` files in a web browser.  In case you can't see your versions, a copy of these files can be found
[here](https://www.well.ox.ac.uk/~gav/projects/oxford_statgen_summer_school/ngs/qc/fastqc/malaria/).  (I've
renamed these so that you can see the sample and machine, and to emphasise there's a file for the first and second read.)

## Interpreting fastqc output

Look at this file now.  Here are a few features:

### Basic statistics

The 'Basic statistics' section tells you some stuff you already know - the number of reads, read
length and so on. It also tells you the average GC content (the proportion of read bases that are
Gs or Cs). This is around 23%, which is about right because the *P.falciparum* genome [has such
high A/T content](https://doi.org/10.1038/nature01097).

### Per base sequence quality

The 'Per base sequence quality' section is interesting. This tells you that the *sequencer's
estimate of base quality* appears lower at the ends of the read. 

**This tailoff is generally seen at the end of the read and to a lesser extent at the start**.

### Per tile sequence quality

The 'Per tile sequence quality' section shows the average per-base quality for each tile on the Illumina
sequencer. In Illumina terminology, a **tile** is a physical region of a sequencer's flowcell lane
that is imaged to detect bases. The chemistry is arranged so that - with luck - that a single DNA
fragment has annealed and amplified within each tile, so that the results from one tile represent a
single 'read'. In this output you can see some tiles that have performed badly (red indicating
below average base quality) which could mean that multiple molecules ended up amplified in the
tile, or that something else prevented good imaging.

### Per-base sequence content

This shows the average proportion of each base in each location in the reads. In this file you see
a fairly typical picture, as follows:

* In general there are the same number of G as C bases, and the same number of A and T bases.
  **This is as it should be**, because either strand of DNA ought to be equally likely to be
  sequenced.
  
* There is nevertheless a skew of base content toward the start of the read. This could be caused
  by a number of things; for example there are thought to be [biases in the sequence content due to
  the chemistry of fragment generation processes](https://doi.org/10.1186/s12859-016-0976-y). Apart
  from this things look pretty flat - as expected.
  
### Per-sequence GC content

Here things get interesting. This is measuring the distribution of G and C bases ('GC' content) in
all reads. P.falciparum has around 24% GC content on average, but it varies a bit around the
genome. So we should expect to see a 'bump' around 24% with some spread around it.

That is indeed what we see some of the other samples in this folder -
like [ERR417621](https://www.well.ox.ac.uk/~gav/projects/oxford_statgen_summer_school/ngs/qc/fastqc/malaria/ERR417621_Illumina-HiSeq_read1_fastqc.html).

But this sample `ERR377582` has a second bump of reads at around 60% GC. (If you look at the [read 2
fastqc output file](https://www.well.ox.ac.uk/~gav/projects/oxford_statgen_summer_school/ngs/qc/fastqc/malaria/ERR377582_Illumina-HiSeq_read2_fastqc.html), you'll see it's there as well.

**Question.** What could this be?

### Sequence duplication levels

In general it's hoped that each 'tile' in the sequencer lane gets a unique fragment, and that the
process of generating these fragments is random, so that we are effectively sequencing millions of
randomly-placed fragments in the genome. If we did that we would generally not see the exact same
sequence twice. However there are a few reasons sequences can end appear duplicated:

* At high depth (which is not the case here - we have about 17-fold coverage) you might get
  multiple reads on the same location just by chance.

* Genomes have repetitive regions (for example centromeres and telomeres tend to be highly
  repetitive). So two reads from different locations might end up identical.

* We are studying whole-genome sequencing here.  But if you had 

In some sense the above are 'everything is ok' reasons.

However in general there are two reasons that 






