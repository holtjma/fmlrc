# FMLRC
## Introduction
FMLRC, or FM-index Long Read Corrector, is a tool for performing hybrid correction of long read sequencing using the BWT and FM-index of short-read sequencing data.  Given a BWT of the short-read sequencing data, FMLRC will build an FM-index and use that as an implicit de Bruijn graph.  Each long read is then corrected independently by identifying low frequency k-mers in the long read and replacing them with the closest matching high frequency k-mers in the implicit de Bruijn graph.  In contrast to other de Bruijn graph based implementations, FMLRC is not restricted to a particular k-mer size and instead uses a two pass method with both a short "k-mer" and a longer "K-mer".  This allows FMLRC to correct through low complexity regions that are computational difficult for short k-mers.

Included in this package are two implementation of the FM-index component of FMLRC.  The default implementation is requires less CPU time but uses a higher sampled FM-index that requires more memory.  The second implementation is more similar to a traditional sampled FM-index that requires less memory, but at the cost of longer computation times.  Both implementation handle parallelization by distributing the reads across to available threads.

## Installation and Setup
First, download the latest version of FMLRC and unzip it.  Then simply make the program and run it with the "-h" option to verify it installed.

    cd fmlrc
    make
    ./fmlrc -h

## Building the short-read BWT
Prior to running FMLRC, a BWT of the short-read sequencing data needs to be constructed.  Currently, the implementation expects it to be in the Run-Length Encoded (RLE) format of the [*msbwt*](https://github.com/holtjma/msbwt) python package.  The *msbwt* package can directly build these BWTs ([Constructing the BWT wiki](https://github.com/holtjma/msbwt/wiki/Constructing-the-MSBWT)) or they can be [Converted to Run-Length Encoded (RLE) format](https://github.com/holtjma/msbwt/wiki/Converting-to-msbwt's-RLE-format) from faster tools like [*ropebwt2*](https://github.com/lh3/ropebwt2).

## Running FMLRC
Once a short-read BWT is constructed, the execution of FMLRC is relatively simple:

    ./fmlrc [options] <comp_msbwt.npy> <long_reads.fa> <corrected_reads.fa>

Here is a partial list of the more useful options of FMLRC:

* -k - sets the length for the short k-mer pass (default: 21)
* -K - sets the length for the long K-mer pass (default: 59)
* -p - sets the number of threads allowed for correction (default: 1)