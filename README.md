# FMLRC
## Update - Version 2
FMLRC has been succeeded by `fmlrc2` ([https://github.com/HudsonAlpha/rust-fmlrc](https://github.com/HudsonAlpha/rust-fmlrc)). 
Preliminary results show `fmlrc2` has near identical results, but runs in <50% of the time. 
It is also implemented in Rust, leveraging the cargo ecosystem. 
Assuming the results are satisfactory for your use case, we recommend switching to `fmlrc2` for the reduced run-time and to receive future updates to the algorithm.

## Introduction
FMLRC, or FM-index Long Read Corrector, is a tool for performing hybrid correction of long read sequencing using the BWT and FM-index of short-read sequencing data.
Given a BWT of the short-read sequencing data, FMLRC will build an FM-index and use that as an implicit de Bruijn graph.
Each long read is then corrected independently by identifying low frequency k-mers in the long read and replacing them with the closest matching high frequency k-mers in the implicit de Bruijn graph.
In contrast to other de Bruijn graph based implementations, FMLRC is not restricted to a particular k-mer size and instead uses a two pass method with both a short "k-mer" and a longer "K-mer".
This allows FMLRC to correct through low complexity regions that are computational difficult for short k-mers.

Included in this package are two implementations of the FM-index component of FMLRC.
The default implementation is requires less CPU time but uses a higher sampled FM-index that requires more memory.
The second implementation is more similar to a traditional sampled FM-index that requires less memory, but at the cost of longer computation times.
Both implementation handle parallelization by distributing the reads across all available threads.

## Quick-start
A full example is available in the `example` subfolder.  Please refer to the [README](https://github.com/holtjma/fmlrc/tree/master/example) for directions.

## Installation and Setup
First, download the latest version of FMLRC and unzip it.  Then simply make the program and run it with the "-h" option to verify it installed.

    cd fmlrc
    make
    ./fmlrc -h

## Building the short-read BWT
Prior to running FMLRC, a BWT of the short-read sequencing data needs to be constructed.
Currently, the implementation expects it to be in the Run-Length Encoded (RLE) format of the [*msbwt*](https://github.com/holtjma/msbwt) python package.
We recommend building the BWT using [*ropebwt2*](https://github.com/lh3/ropebwt2) by following the instructions on [Converting to the fmlrc RLE-BWT format](https://github.com/holtjma/fmlrc/wiki/Converting-to-the-fmlrc-RLE-BWT-format).
Alternatively, the *msbwt* package can directly build these BWTs ([Constructing the BWT wiki](https://github.com/holtjma/msbwt/wiki/Constructing-the-MSBWT)), but it may be slower and less memory efficient.

## Running FMLRC
Once a short-read BWT is constructed, the execution of FMLRC is relatively simple:

    ./fmlrc [options] <comp_msbwt.npy> <long_reads.fa> <corrected_reads.fa>

Here is a partial list of the more useful options of FMLRC:

* -k - sets the length for the short k-mer pass (default: 21)
* -K - sets the length for the long K-mer pass (default: 59)
* -p - sets the number of threads allowed for correction (default: 1)

## Reference

[Wang, Jeremy R. and Holt, James and McMillan, Leonard and Jones, Corbin D. FMLRC: Hybrid long read error correction using an FM-index. BMC Bioinformatics, 2018. 19 (1) 50.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2051-3)
