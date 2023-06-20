---
sidebar_position: 8
---

# Dataset E: 10X linked-read sequencing

**Dataset link**:
```
https://tinyurl.com/3chwns49/10x.bam
```

Yet another type of data is linked-read sequencing. This is Illumina short-read sequencing, but including molecular 'barcodes' that allows the technology to link the reads together.

:::tip IGV hint

To see the barcode linking, you need to turn on 'Linked read view (BX)' in the context menu. (To make this easy to see, you may also then need to choose 'Group alignments by -> None' and then 'Collapsed' from
the same menu.)

At this point you should be able to see multiple reads linked together by thin lines - the linked reads all come from the same DNA fragment.
:::
