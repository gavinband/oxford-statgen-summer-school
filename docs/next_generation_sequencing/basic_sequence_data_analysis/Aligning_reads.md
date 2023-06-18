---
sidebar_position: 5
---

# Step 3: Aligning reads

Most analyses of short-read sequence data are based on the basic paradigm made up of two steps:

* first, *align reads to an established reference genome assembly*.
* then *use those alignments to identify and call genetic variation for analysis*.

The first step in this paradigm is to get reads aligned. On this page we'll implement the basic
read alignment pipeline that starts with the fastqs and ends with alignments for downstream
analysis. 

## The basic pipeline

The basic pipeline is:

1. align reads to the reference genome assembly
2. sort them so they are in genome position order
3. 'mark' reads that look like duplicates (since these are [likely to be artifacts](Quality_control.md).)
4. Finally we'll gather some statistics about the alignments.

In theory that's it - in practice it gets a bit more fiddly because of the need for various format
conversion and indexing steps. (Some data might also need a [read trimming step](Read_trimming.md)
but we're skipping that.)

Without further ado let's implement these steps using a popular aligner `bwa mem` and the swiss
army knife-like `samtools` to do the work. As before I'll work with the data for the `QG0033-C`
malaria sample (which was ascertained in Congo in 2012), but if you know what you are doing you are
free to work with any of the other samples in the folder - or the chromosome 19 human data in the
`human/` folder if you prefer. (Just remember to use the human reference sequence `GRCh38_chr19.fa`
in that case.)

To get started, move to a terminal window and make sure you are in the
`sequence_data_analysis/malaria` folder:

```
cd ~/sequence_data_analysis/malaria
```

The pipeline can be run like this:
```sh
# Make a temp dir to hold intermediate files
mkdir -p tmp

# Create a FM index for the reference
bwa index Pf3D7_v3.fa

# align the reads (using 2 threads - don't use more please!)
bwa mem -t 4 -o tmp/QG0033-C-aligned.sam Pf3D7_v3.fa QG0033-C_Illumina-HiSeq_read1.fastq.gz QG0033-C_Illumina-HiSeq_read2.fastq.gz

# convert to BAM
samtools view -b -o tmp/QG0033-C-aligned.bam tmp/QG0033-C-aligned.sam

# fix mate-pair information and convert to BAM
samtools fixmate -m tmp/QG0033-C-aligned.bam tmp/QG0033-C-fixmate.bam

# sort reads by genomic position
samtools sort -T tmp -o tmp/QG0033-C-sorted.bam tmp/QG0033-C-fixmate.bam

# Identify probable duplicate reads
samtools markdup -s tmp/QG0033-C-sorted.bam tmp/QG0033-C-markdup.bam

# Rename to this folder, and index
mv tmp/QG0033-C-markdup.bam QG0033-C.bam
samtools index QG0033-C.bam

```

(**Note.** You can see what this has created by running `ls`.)

If you run all that together it might take around 10 minutes to complete (the longest steps being the `bwa mem`
indexing and alignment steps.) It's worth running the steps one by one to make sure they complete without errors, and
watching the output to see what is being done. While you're waiting for the alignment to complete, consider reading the
information below about what's in the aligned data, [reading more about paired-end sequencing theory](Short_read_theory.md),
or more about the SAM/BAM format and using samtools.

For example here is a [useful overview of SAM format](https://davetang.org/wiki/tiki-index.php?page=SAM), or a page
with [more tips on using samtools](https://github.com/davetang/learning_bam_file).


### What did the pipeline generate?

If you followed the above you'll now have:

- a file `tmp/QG0033-C-aligned.sam` which contains the alignments in the [SAM text file format](https://en.wikipedia.org/wiki/SAM_(file_format)).

- a binary version `tmp/QG0033-C-aligned.bam` of the same thing, in [BAM format](https://en.wikipedia.org/wiki/Binary_Alignment_Map).

- a coordinate-sorted version of that, in `tmp/QG0033-C-sorted.bam`

- a final version of the alignments with duplicate reads 'marked' - called `QG0033-C.bam`

- and a corresponding index file `QG0033-C.bam.bai`.

**Congratulations!** You should now have an aligned set of reads ready for analysis. 

The rest of this page gives details about what is in the alignment files you've just created. When
you're satisfied you know enough, [go back to the practical](Pipeline_outline.md#the-practical-in-a-nutshell).

## Inspecting the alignment output

### The contents of a SAM file

If you look at the unsorted SAM file (e.g. by typing `less -S tmp/QG0033-C-aligned.sam` - press `q` when you want to quit) you'll see something like this:

```
@SQ     SN:Pf3D7_08_v3  LN:1472805
@SQ     SN:Pf_M76611    LN:5967
@SQ     SN:Pf3D7_01_v3  LN:640851
@SQ     SN:Pf3D7_09_v3  LN:1541735
@SQ     SN:Pf3D7_03_v3  LN:1067971
@SQ     SN:Pf3D7_05_v3  LN:1343557
@SQ     SN:Pf3D7_02_v3  LN:947102
@SQ     SN:Pf3D7_14_v3  LN:3291936
@SQ     SN:Pf3D7_11_v3  LN:2038340
@SQ     SN:Pf3D7_10_v3  LN:1687656
@SQ     SN:Pf3D7_04_v3  LN:1200490
@SQ     SN:Pf3D7_12_v3  LN:2271494
@SQ     SN:Pf3D7_13_v3  LN:2925236
@SQ     SN:Pf3D7_06_v3  LN:1418242
@SQ     SN:Pf3D7_07_v3  LN:1445207
@SQ     SN:PF_apicoplast_genome_1       LN:29430
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem -t 2 -o QG0033-C-aligned.sam Pf3D7_v3.fa QG0033-C_Illumina-HiSeq_read1.fast
ERR377582.7615542       99      Pf3D7_09_v3     37181   21      100M    =       37244   163     AAAAATCCATTTATATCTTTTATGGTTAGTATTATTTATACCTTCCTGTTTTATTTATCCACGTTTTATAAAATTGACGTATTGTTAATAAGTAGTTGTA    B@DECEFEEGFEEHFGGEFFFFFEFFEEFEGHHGGGFFFAEFFFGGEFFHFGEEFFEDEGDDEF@GFGEGEEGHEEDGGGEFEEFFEEEDEFDGDECEEF    NM:i:1  MD:Z:58T41      MC:Z:100M       AS:i:95 XS:i:91
ERR377582.7615542       147     Pf3D7_09_v3     37244   21      100M    =       37181   -163    TTTTATAAAATTGACGTATTGTTAATAAGTAGTTGTATATTATAAATTATAAGGTAATATAATTGCTATAATTATAATAATCATTATCATAATTATTATC    DHEDGBHEDGFEDAAGEDEDGFHGGGDEFCFGFGFFEEEIBCFCIFDEEDFFGGBEFDGEEGEFFEFHIFDEFEHHEFHEDDFFEEFDFEEFDEDECA@9    NM:i:2  MD:Z:43G49A6    MC:Z:100M       AS:i:90 XS:i:90
ERR377582.19970658      77      *       0       0       *       *       0       0       TGATCGGTGCGGGCCGGACTGAATTGCTACGCCTGATTTTCGGTGCCGACCTGGCCGACAGTGGCACGGTGGCCTTGGGGTCGCCAGCGCAGGTGGTGAG    *:?D<CF3-D@EE@EDDGFF?DDEFACEFEBHCGD/EFEE:?E;F?GE=2FE>FFC:EHD??H@DFFB89CD@D-ECC8;EF.A:;6EEDDFA,D>'D@F    AS:i:0  XS:i:0
...
```

The file consists of some *metadata* (lines starting with `@`, and recording information about the sequences in the reference
assembly and the program that was run to generate the alignments) followed by the alignments themselves.

The alignment columns rows consist of

* The read ID (straight from the fastq file)
* Some *flags* encoded as an integer.  The best way to understand these is through the ["Explain SAM flags" webpage](https://broadinstitute.github.io/picard/explain-flags.html).  We'll go through some examples below.
* The chromosome and position at which the read aligns
* The *mapping quality*, detailed more below.
* The *CIGAR* string which gives some detail about how the read aligns
* The chromosome and position of the other read in the pair *if it's aligned*. (The
chromosome is often `=`, which means the same chromosome as this read, which you'd expect since they both came from the same fragment. A `*` means the other read was not aligned).
* The 'template length' (how long was the span of the read pair on the reference contig?)
* The read bases and mapping qualities.  These are straight from the fastq but *reverse-complemented if necessary* so they are in the same order as the reference bases.
* Finally there are some **optional fields**.  Among these are **NM** (number of mismatches between the read and the reference), and **MD** (a kind of counterpoint to the CIGAR string in column 6).

If you look at the read IDs and flags above (via
[this page](https://broadinstitute.github.io/picard/explain-flags.html)), you'll see that the read pairs
are represented as pairs of alignments with the same identifiers but with one flagged as 'first in
pair' and the other as 'second in pair'.

### How many alignments are there?

We can use the UNIX `wc` command to count the alignments:
```
samtools view tmp/QG0033-C-aligned.sam | wc -l
```

**Note.** The BAM file is exactly the same, it's just encoded in a binary format.  So this also works and is slightly quicker:
```
samtools view tmp/QG0033-C-aligned.bam | wc -l
```

If you run this with the above data it says there are 4,057,972 alignments. But wait - there were only 2 million read pairs!

This illustrates an important point: **there are more alignments than reads in the data**.  How can this be?

### Primary and supplementary alignments

To figure out why, let's look at an example: the read with ID `ERR377582.20226793`. You can extract
this using `grep` like this:

```
grep 'ERR377582.20226793' tmp/QG0033-C-aligned.sam
```

If you look at the alignment flags (second columns) via [this page](https://broadinstitute.github.io/picard/explain-flags.html), you'll see that:

* The first alignment (flag=81) represents the first read in the pair.  It is aligned on the reverse strand of the reference.
* The second alignment (flag=2113) *also* represents the first read in the pair.  It is flagged as a **supplementary alignment**.
* The third alignment (flag=161) represents the second read in the pair.

**Terminology**: A *supplementary alignment* occurs when a read aligns in more than one part. (The
additional alignments aren't adjacent to the primary one, otherwise they would be part of it.)

You can see how this works for this read by looking at the *CIGAR string*. For the first alignment
the CIGAR is `60M40S`. This means that the first 60 bases of the 100bp read were aligned, and the
following 40 were 'soft clipped' - not part of the alignment. For the supplementary alignment, the
CIGAR is `54M46H`. If you look at the read bases you will see they are the reverse complement of
the bases from the original read (which are fully represented in the first alignment). The upshot
is that `bwa` thinks the first ~60 bases align to chromosome `Pf3D7_10_v3`, but the last ~54 bases
align in reverse orientation to chromosome `Pf3D7_09_v3`.

### Interpreting mapping qualities

Taken at face value you might think that this means this sample has a translocation between
chromosome 9 and 10. But hold on. If you look at the *mapping quality* column (column 5) you'll see
that all these alignments actually have **zero mapping quality**. This means that confidence in
these alignments is very low!

Like base qualities, mapping qualities are expressed on the PHRED scale. They are the aligner's
estimate of the quantity: the *probability that the read actually aligns elsewhere*. On the PHRED
scale this is expressed as:

$$
\text{mapping quality} = -10 \log_{10} \left( P( \text{the read actually maps elsewhere}) \right)
$$

which gives this relationship:

| mapping quality | Estimated probability alignment is wrong |
| --------------- | ---------------------------------------- |
| 0               | 100%                                     |
| 10              |  10%                                     |
| 20              |   1%                                     |
| 30              |   0.1%                                   |
| 40              |   0.001%                                 |

A mapping quality of zero translates to a 100% probability that the read aligns elsewhere, that is,
that `bwa` is not confident at all in this alignment. 

**Note.** Also like base qualities, this is only an estimate of the probability.  It may not be well calibrated.

(What's really going on is that these reads original from the telomeres. Telomeres in malaria, as
in humans and other organisms, are highly repetitive and are very hard to analyse using short-read
data.)

### What types of alignment are there?

To get a better sense of the types of alignment let's use `samtools flagstat`.  
This counts reads according to the `flags` column. Run it like this:

```
samtools flagstat tmp/QG0033-C-aligned.sam 
```

From the output you should see: 

* Of 4 million reads ("4000000 paired in sequencing"), only 3,804,730 were actually aligned by `bwa`.  That's nearly 200,000 that weren't aligned at all.
* Most of these (3731534) were were aligned in pairs.
* Most of *these* (3537816) were 'properly paired' - that is, both reads mapped to the same chromosome, in the right orientation (i.e. [facing each other](Short_read_theory.md)) and about the right distance apart.
* But a small subset (264863 in total) had unmapped mate, mate on a different chromosome, or not close or in the wrong orientation to the original read.

**Note.** A much more detailed view of reads can be generated by running `samtools stats`. See more
on the [`samtools stats` documentation page](http://www.htslib.org/doc/samtools-stats.html).

### What did `markdup` do?

**Question** what **did** the `samtools markdup` step do?

**Hint 1**. Try viewing the first reads in the coordinate-sorted and final (post-markdup BAM).  E.g.:

```
samtools view tmp/QG0033-C-sorted.bam  | head
samtools view QG0033-C.bam  | head
```

Can you spot the difference

**Hint**. Try running `samtools flagstat` on these two files.  What's different?

### Important!  Cleaning up

Our pipeline has left behind several copies of the original data that we don't need.  See them by looking in the temp directory:

```
ls -lh tmp
du -ch tmp
```

That's 2.4Gb of space essentially wasted!  Now is a good time to get rid of these:
```
rm tmp/*
rmdir tmp
```

**Note.** **Don't** delete the main output file `QG0033-C.bam` or its index - we'll look at that in a moment.

## Next steps

Now that we have an aligned set of reads, move on to [viewing your alignments](./Viewing_alignments.md).
