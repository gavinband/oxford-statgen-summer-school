# Inspecting the fastqs

## First look at the fastqs

Let's take a look at the data using the UNIX command `less` (actually we'll use `zless`, as this
will uncompress the file on the fly):

```
zless -S ERR377582_1.fastq.gz
```

**Note**. the `-S` tells `zless` to stretch long lines off to the right of the screen.  You can navigate using the arrow keys.

Look at the structure of the fastq file. Each read occupies four lines in the file, for example the
first read in the above file looks like:

```
@ERR377582.7615542 HS23_10792:2:2307:6524:31920#15/1
AAAAATCCATTTATATCTTTTATGGTTAGTATTATTTATACCT...
+
B@DECEFEEGFEEHFGGEFFFFFEFFEEFEGHHGGGFFFAEFF...
```

**Note**. you can press `q` at any time to exit `zless` - but don't do that yet.  

### Warmup questions

Here are some basic questions to ask about the FASTQ files.  Can you answer them? (Hints below):

* **Question 1**: how many reads are in the file?
* **Question 2**: how long are the reads?

So this should let you calculate:

* **Question 3**: if the *P.falciparum* genome is about 23 Mb long, what sequencing depth do you expect to get?

### Warmup question hints

There are several ways to answer these questions, but the quickest and easiest involve
using basic UNIX command line tools - notably `wc`, which counts lines or characters in its input,
and `head` and `tail` which isolate the top or bottom lines. So to count the number of reads you
could do:

```
zcat ERR377582_1.fastq.gz | wc
```

This might take a minute or so to run - it is decompressing the whole file of course. The number
output is the number of lines in the file.  (So what is the number of reads?)

And to count the read length you could do:
```
zcat ERR377582_1.fastq.gz | head -n 2 | tail -n 1 | wc -c
```
(The number output is the number of characters in the second line of the file, i.e. the read length.)

### Next steps

More details on what's in the FASTQ file are below. When you're satisfied you know what is in the
files, move on to [read quality control](quality_control.md). Or

## More details
### The fastq read header line

The first row of each fastq record contains a whole bunch of information:

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

### The read bases line

The second row of each record contains the read bases themselves:
```
AAAAATCCATTTATATCTTTTATGGTTAGTATTATTTATACCT...
```
These are the DNA bases as called by the sequencer, in the order they were sequenced.  Simple!

Sounds obvious but the length of this line shows you the *read length*.  So you can compute the read length by:
```
zcat ERR377582_1.fastq.gz | head -n 2 | tail -n 1 | wc
```

The third number output is the number of characters, that is, the read length.

#### What's in the base quality line?

For each of the above bases, the sequencer also emits an estimate of the *base quality*.  These look a bit confusing at first:
```
B@DECEFEEGFEEHFGGEFFFFFEFFEEFEGHHGGGFFFAEFFFGGEFFH...
```

The computation goes:

* start with an estimate $x$ of the *base error*, that is, the probability that the base is incorrectly called.

* transform it to PHRED scale (that is, take $-10 * log_{10} (x)$). This brings it into the range
  from zero (very likely to be an error) to infinity (no chance at all of being an error), although
  in practice it generally maxes out at 41 (except for PacBio Hifi).

* Finally add 33 and take the corresponding [ASCII character](https://en.wikipedia.org/wiki/ASCII).

The point of this is that ASCII characters from 33 onwards are visible / printable, so this
generates a sequence of characters that encode the estimated quality of each base. (**Note.** There
is also a 'PHRED64' encoding which adds 64 instead. It's generally only used on older sequencers
and you aren't likely to run across it for new data.)

It is important to realise that these qualities are only the values *estimated* by the sequencer.
They are generally [pretty good](https://lh3.github.io/2017/07/24/on-nonvaseq-base-quality) but can
be [improved by re-estimation](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-) (but that's beyond the scope of this tutorial.)

## Next steps

Now move on to [read quality control](quality_control.md).
