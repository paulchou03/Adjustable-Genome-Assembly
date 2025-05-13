# Adjustable-Genome-Assembly using de Bruijn Graphs

This project implements a genome assembler that reconstructs the underlying genome sequence from a set of short sequencing reads using a **de Bruijn graph** approach. It features adjustable k-mer frequency thresholds and mismatch tolerance for flexible assembly accuracy and read alignment sensitivity.

---

## Project Overview

Given a FASTA file of short DNA reads, the assembler:

1. Extracts k-mers and filters them based on customizable frequency thresholds.
2. Constructs a de Bruijn graph from the filtered k-mers.
3. Traverses the graph to assemble the genome by identifying an **Eulerian path**.

This method supports control over:

- **k-mer inclusion thresholds**
- **Repeat resolution**
- **Read-to-genome mismatch tolerance**

---

## Input

The assembler expects an input FASTA file called `reads.fasta`, formatted as:

```fasta
>read_0
ACGTACGTGACTAGCT...
>read_1
TACGTGACTGGTAGCA...
...
```

## Parameters

The following parameters can be configured in the script to control the behavior of the genome assembler:

- `readLength` *(int)*  
  Represents the number of base pairs per read (e.g., 50, 150, 300, etc.).  

- `k` *(int)*  
  Length of k-mers used to construct the de Bruijn graph. Affects assembly sensitivity and specificity.  

- `max_mismatches` *(int)*  
  Maximum number of mismatches allowed when aligning a read back to the assembled genome.  

- `low_thresh` *(int)*  
  Minimum frequency a k-mer must have to be included in the graph. Filters out sequencing errors and low-confidence k-mers.  

- `high_thresh` *(int)*  
  Frequency above which a k-mer is considered repetitive. Such k-mers are downweighted using the `step` value.  

- `step` *(int)*  
  Determines how much to increase the count weight for overrepresented k-mers (used to adjust path traversal frequency).  

##  Output

After assembling the genome and aligning reads with mismatch tolerance, the script generates a file called `output.txt`.

### `output.txt`

This file contains the subset of input reads that successfully align to the assembled genome (with up to `max_mismatches` differences).

Each aligned read is output in FASTA format:

```fasta
>read_3
ACGTACGTGACTAGCTAGTTACGGAGGACCT...
>read_17
TACGTGACTGGTAGCACAGGTTAGCGTCAG...
...
