# Comb aligner

Comb aligner is a prototype implementation of a graphical local alignment algorithm, calculating a set of high-scoring paths between two nucleotide string graphs. The aligner accepts FASTA, FASTQ, and [GFA](https://github.com/GFA-spec/GFA-spec) format for the input (query and reference) sequences and dumps alignments in the [SAM](https://github.com/samtools/hts-specs) or a [Graphical Pairwise Alignment (GPA)](https://github.com/ocxtal/gpa) format.


## Overview of the algorithm

The algorithmic design basically follows the seed-and-extend method. Input files are parsed with [libfna](https://github.com/ocxtal/libfna) into individual sequence segments and links. The set of parsed sequences and links are stored in [gref](https://github.com/ocxtal/libgref) objects, then indexed by k-mer hashing. [Libggsea](https://github.com/ocxtal/libggsea), which takes two gref objects, iterates all the k-mer over one object and matches them on the other one. The matched seeds are extended upwards and downwards with [libgaba](https://github.com/ocxtal/libfna) in an adaptive banded way, where the graphs are dynamically expanded into trees with their root at the seed positions and traversed in a breadth-first order. The alignment paths are generated between pairs of maximum score positions in the trees and dumped into the SAM or GPA format with [libaw](https://github.com/ocxtal/libaw).


## Build

Python (>= 2.7 or 3.3) is required to run the build script dependent on the waf build framework. The programs, entirely written in C99, can be compiled with gcc-compatible compilers passing the CC option to the waf configure argument.

```
./waf configure CC=clang clean build
sudo ./waf install
```

## Usage

```
comb [options] <ref. seq.> <query seq.> <output>
```

### Options

#### Indexing options

* **-k\<int\>** k-mer and seed length

#### Scoring options

A penalty for a gap with length k is represented in a Gi + k*Ge form.

* **-a\<int\>** match award M in a positive integer
* **-b\<int\>** mismatch penalty X in a positive integer
* **-p\<int\>** gap open penalty Gi in a positive integer or zero
* **-q\<int\>** gap extension penalty Ge in a positive integer
* **-x\<int\>** xdrop threshold in a positive integer
* **-m\<int\>** minimum score for reporting

#### Others

* **-h** print help
* **-v** print version info
* **-t** number of threads (multi-thread mode is broken for now)


## Examples

### Linear-to-linear

Linear-to-linear alignment is the normal, conventional alignment concept that is implemented in many programs like the BLAST, BWA and so on... The comb aligner can handle the similar tasks, taking FASTA reference and FASTA/Q reads then generating SAM output.

```
comb example/linear1.fa example/linear2.fa test.sam
```

### Linear-to-graph

The GFA format is acceptable as an input reference sequence object. The sequence segments in the GFA file are indexed in the same way as the linear references. The alignments are reported in the GPA format by default.

```
comb example/graph1.fa example/linear2.fa test.gpa
```

### Graph-to-graph

The comb aligner can find a set of high-score paths between two graphical sequence objects.

```
comb example/graph1.fa example/graph2.fa test.gpa
```

## Notes

### Issues

* The software is not stable. It may report wrong results and segfaults in an unexpected way.
* The libfna parser library and libgref graph indexing library cannot handle links with overlaps (links with non-0M cigar in the GFA format).
* The affine-gap penalty alignment routine in the gaba library has a bug, reporting wrong alignment paths.
* Overlapping results are not filterd.
* Multi-thread mode segfaults.

### TODO

* Fix known bugs listed above.
* Add VCF parser to enable SNP and short indel modifications. It also requires implementing two functions, `append_snp` and `split_segment` in the gref library.
* Add seed filtering to improve performance.
* Add matrix merge to reduce computational complexity.


## License

Apache v2.

Copyright (c) 2016, Hajime Suzuki
