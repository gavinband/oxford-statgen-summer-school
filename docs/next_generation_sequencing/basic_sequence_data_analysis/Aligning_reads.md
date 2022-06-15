# A basic read alignment pipeline.

In this section we will run a basic read alignment pipeline.  The pipeline will:

1. use `bwa mem` to align each read pair to a reference genome assembly
2. use `samtools` to identify probable duplicate read pairs
3. perform various sorting and indexing steps to make the data useful.
4. Use `samtools stats` to gather some statistics about the alignment

**Note.** A full pipeline might also implement additional steps - such as **read trimming** to
remove trailing low-quality bases, adapter sequence at the 3' end of reads, or perhaps trailing Gs
reflecting read-through the fragment end. Tools like
[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic),
[cutadapt](https://cutadapt.readthedocs.io/en/stable/), or
[fastp](https://github.com/OpenGene/fastp) can be used to do this.  I'll come back to this below.

Let's get started.

## Using bwa to align short-read data

Here are the steps we have to take:

* Create an index of the reference assembly for the aligner to work with.
* Align reads to the reference assembly.
* Convert the output file to the binary [BAM format](https://en.wikipedia.org/wiki/Binary_Alignment_Map).
* Convert the SAM to a binary BAM format (because this uses less space)
* Sort the file by genomic position
* 'Mark' duplicate reads.
* Finally we index the resulting file.  

Without further ado let's implement these steps. We are working in the bash terminal.) I recommend
running these commands one at a time. In particular, I recommend watching the output of the `bwa
mem` command as it does its work:

```sh
# Create a FM index for the reference
bwa index Pf3D7_v3.fa.gz

# align the reads
bwa mem -o ERR377582-aligned.sam Pf3D7_v3.fa.gz ERR377582_1.fastq.gz ERR377582_2.fastq.gz

# convert to BAM
samtools view -b -o ERR377582-aligned.bam ERR377582-aligned.sam

# sort reads by genomic position
mkdir tmp
samtools sort -T tmp -o ERR377582-sorted.bam ERR377582-aligned.bam

# Mark probable duplicate reads
samtools markdup ERR377582-sorted.bam ERR377582-markdup.bam

# Rename and inex
mv ERR377582-markdup.bam ERR377582.bam
samtools index ERR377582.bam

```

### Inspecting alignments

Did it work?  We can use `samtools tview` to have a look at some reads.  Let's look at a particular gene, 

```

```