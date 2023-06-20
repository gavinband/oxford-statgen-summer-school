---
sidebar_position: 5
---

# Dataset E: ATAC-seq

**Dataset E URL**:
```
https://tinyurl.com/3chwns49/ATAC-monocytes.bam
```

Dataset E is [ATAC-seq](https://en.wikipedia.org/wiki/ATAC-seq) from monocytes.  ATAC-seq is a way of measuring regions of open chromatin.
It uses a transposase to introduce sequencing primers specifically into DNA in open regions, which are then sequenced.

Open chromatin is found upstream of transcribed genes, whereas silenced genes tend to have closed chromatin.

:::caution Warning

A reminder that all data in this sightseeing tour is included strictly for training purposes - it is **not**
publicly-available data. Please do not share outside this course.
Contact me (Gavin Band) if you have any queries about this.

:::

## Questions

* Look at the ATAC, RNA, and long-read 5mC data around some genes.  Do the patterns make sense?

You might roughly expect to see:

- open chromatin (ATAC peaks) near the start of genes, if they are transcribed.
- if the DNA is open, it is probably unmethylated (blue)
- if the gene is expressed, it should have reasonable coverage of RNA-seq reads over the exons.

:::caution warning

The ATAC-seq and RNA-seq here come from monocytes (or more exactly - they come from cells isolated due to expressing
CD14 on the cell surface. What does the gene CD14 look like in these data?)

However the genomic sequencing is from peripheral blood mononuclear cells, i.e. PBMCs. These are a mixture of cell
types including a majority of T cells, B cells, as well as monocytes. Because of this the patterns of methylation, open
chromatin, and expression may not perfectly correspond.

:::
