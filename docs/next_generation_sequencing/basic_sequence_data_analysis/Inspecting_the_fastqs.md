---
sidebar_position: 3
---

# What's in the fastq files?

### Using UNIX to inspect the fastqs

The simplest thing to do is look at the data in the terminal. The command `less` will do this for
us - actually we'll use `zless` because the data is compressed).  Try it now:

```
zless -S malaria/QG0033-C_Illumina-HiSeq_read1.fastq.gz
```

**Note**. You ought to be able to copy/paste the above into your terminal, or else type it
directly. (Use &lt;tab&gt; to auto-complete filenames.) The `-S` option tells `zless` to stretch
long lines off to the right of the screen. You can navigate using the arrow keys. You can also quit
any time by pressing `q`.

Look at the structure of the fastq file. **Each read occupies four lines in the file**.  For example the
first read in the above file looks like:

```
@ERR377582.7615542 HS23_10792:2:2307:6524:31920#15/1
AAAAATCCATTTATATCTTTTATGGTTAGTATTATTTATACCT...
+
B@DECEFEEGFEEHFGGEFFFFFEFFEEFEGHHGGGFFFAEFF...
...
```

In the above:

* The first line is a header line that identifies the read.
* The second line is the read itself
* the third line is just a `+`
* the fourth line contains base qualities in "PHRED" encoding.

(More details on what this all means are given below, but this is the basic idea.)

### Warmup questions

Here are some basic questions to ask about the FASTQ files.  Can you answer them? (Fear not - hints are below):

* **Question 1** How does the header differ between the *first read in the read1 file* and the
*first read in the read2 file*?
* **Question 2**: how many reads are in the file?
* **Question 3**: how long are the reads?

The answers to Questions 2-3 should let you calculate an important quantity:

* **Question 4**: if the *P.falciparum* genome is about 23 Mb long, what sequencing depth do you expect to get?

#### Questions hints

There are several ways to answer these questions, but the quickest and easiest involve using basic
UNIX command line tools - notably `zcat` (which decompressed the file), `wc`, which counts lines or
characters in its input, and `head` and `tail` which isolate the top or bottom lines. Here we go:

* To look at the header of the first read you can use `head` (we also need to decompress the file using `zcat`).  So:

```
zcat malaria/QG0033-C_Illumina-HiSeq_read1.fastq.gz | head -n 1
zcat malaria/QG0033-C_Illumina-HiSeq_read2.fastq.gz | head -n 1
```

Spot the difference?

* To count the number of reads you could do:

```
zcat malaria/QG0033-C_Illumina-HiSeq_read1.fastq.gz | wc -l
```

This might take a minute or so to run - it is decompressing the whole file of course. The number
output is the number of lines in the file - so what is the number of reads?

* To count the read length you could inspect the first read - which is on the 2nd line of the file:
```
zcat malaria/QG0033-C_Illumina-HiSeq_read1.fastq.gz | head -n 2 | tail -n 1 | wc -c
```
The number output is the number of characters in the second line of the file, i.e. the read length.

**Note.** If you don't yet have facility with these types of command - don't worry, you will gain
it.  However, don't worry because we'll move to using more sophisticated graphical packages in the next step.

When you're ready to move on, [continue the practical](Pipeline_outline.md#the-practical-in-a-nutshell).

## More details on the FASTQ format
### Structure

As described above, each FASTQ record spans four lines: the **header**, the **read sequence**, a
**separator**, and the encoded **base qualities**.

### The read header

The header row typically contains a whole bunch of information.  For example:

```
@ERR377582.7615542 HS23_10792:2:2307:6524:31920#15/1
```

This tells us the sample ID (`ERR377582`) and the read identifier (`7615542`). And this is
followed by information identifying the instrument that generated the reads (`HS23_10792`), the
flowcell lane and tile number in the lane (`2:2307`), the coordinates of the
[cluster](https://www.broadinstitute.org/files/shared/illuminavids/clusterGenSlides.pdf) within the
tile (`6524`, `31920`), a number identifying the index of the sample within a multiplexed set of
samples (i.e. all run at the same time; `#15`). And finally it identifies whether it's read 1 or 2
(`/1` or `/2`).

(If you look in the second file in the pair, you should see exactly the same IDs in the same order,
except that the end `/2` will be present indicating the second read.)

**Note.** The format of this information is not standard across platforms, and it changes depending
on your data provider. Some other examples can be found [on
wikipedia](https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers) or on the
[GATK read groups page](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups).

### The read sequence

The second row of each record contains the read bases themselves:
```
AAAAATCCATTTATATCTTTTATGGTTAGTATTATTTATACCT...
```
These are the DNA bases as called by the sequencer, in the order they were sequenced.  Simple!

Sounds obvious but the length of this line shows you the *read length*.  So you can compute the read length by:
```
zcat file.fastq.gz | head -n 2 | tail -n 1 | wc -c
```

The number output is the number of characters in the line, that is, the read length.

#### The base qualities

For each of the above bases, the sequencer also emits an estimate of the *base quality*. These look
a bit confusing at first:

```
B@DECEFEEGFEEHFGGEFFFFFEFFEEFEGHHGGGFFFAEFFFGGEFFH...
```

The computation goes:

* start with an estimate $x$ of the *base error*, that is, the probability that the base is incorrectly called.

* transform it to PHRED scale (that is, take $-10 * log_{10} (x)$). In principle this brings it
  into the range from zero (very likely to be an error) to infinity (no chance at all of being an
  error). However, in practice it generally maxes out at 41 (except for PacBio Hifi which goes up to ~100).

* Finally add 33 and take the corresponding [ASCII character](https://en.wikipedia.org/wiki/ASCII).

The point of this is that ASCII characters from 33 onwards are visible / printable, so this
generates a sequence of characters that encode the estimated quality of each base. (**Note.** There
is also a 'PHRED64' encoding which adds 64 instead. It's generally only used on older sequencers
and you aren't likely to run across it for new data.)

It is important to realise that these qualities are only the values *estimated* by the sequencer.
They are generally [pretty good](https://lh3.github.io/2017/07/24/on-nonvaseq-base-quality) but can
be [improved by re-estimation](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-) (but that's beyond the scope of this tutorial.)

## Next steps

When you're ready to move on, [continue the practical](Pipeline_outline.md#the-practical-in-a-nutshell) (go on to step 2).
