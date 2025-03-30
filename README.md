# Adjustable-Genome-Assembly using de Bruijn Graphs

This project implements a genome assembler that reconstructs the underlying genome sequence from a set of short sequencing reads using a **de Bruijn graph** approach. It features adjustable k-mer frequency thresholds and mismatch tolerance for flexible assembly accuracy and read alignment sensitivity.

---

## Project Overview

Given a FASTA file of short DNA reads, the assembler:

1. Extracts k-mers and filters them based on customizable frequency thresholds.
2. Constructs a de Bruijn graph from the filtered k-mers.
3. Traverses the graph to assemble the genome by identifying an **Eulerian path**.
4. Aligns the assembled genome back to the reads, allowing for mismatches.

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
