PathEncode - Statistical reference-based compression for short reads

About Path Encoding
===================

Path encoding is a technique for compressing short-read sequence files. It uses
a reference (any gzipped mulit-FASTA file) to build a statistical model of the
sequences, which is adaptatively updated during compression.

The path encoding software is written in Go, and is open source.

If you use this software, please cite:

    Carl Kingsford and Rob Patro. Compression of short-read sequences using
    path encoding. Under review (2014).  
    

Installation & Requirements
===========================

Binaries for Mac OS X and Linux are available below. If neither of them work on
your system, you can easily build the software from the sources.


To install using the binary:
----------------------------

* Download the latest version of kpath

* Decompress the tarball: tar xzf kpath-0.6.1.tar.gz

* Copy the kpath-X.X.X-OS binary to a location in your path (for easy access)


To install using the source:
----------------------------

* Install Go, version 1.2 or 1.3 (1.3 recommended) [this is easy].

* Download the latest version of kpath

* Decompress the tarball: tar xzf kpath-0.6.1.tar.gz
  and copy the "src" directory and its subdirectories into the "src"
  subdirectory of your GOPATH workspace. 

  Alternatively, make the directory that tar created the root of your Go
  workspace:

	export GOPATH=/path/to/kpath-0.6.1/

* Compile: 
	cd kpath-0.6.1/src/kingsford/kpath 
	go build

* Copy the kpath executable that is created to a location in your path (for
  easy access)


Usage
=====

To compress:
------------

    kpath encode -ref=REF -reads=IN.fastq -out=OUT

where REF is the path to a gzipped multi-fasta file containing your reference
sequences (i.e. a set of transcripts, or genomes, or chromosomes); IN.fastq is
the fastq file you want to compress; OUT is the prefix of the output files
where compressed version are stored.  kpath will create OUT.enc, OUT.bittree,
OUT.counts, OUT.flipped, and OUT.ns. The first three files (.enc, .bittree,
.counts) are needed to decompress the sequences if you don't care about Ns the
orientation of the reads. You can delete one or both of .flipped and .ns.


To decompress:
--------------

    kpath decode -ref=REF -reads=OUT -out=RECOVERED.fasta

where OUT is the basename for the file to decompress (same as the OUT used when
encoding). REF is the same reference used when encoding. This will write the
sequences in FASTA format to RECOVERED.fasta. If the OUT.ns file is present the
Ns in the original reads will be recovered. If the OUT.flipped file is present,
the reads will be put in their original orientation. If either of these files
is missing, the corresponding step will be skipped.  The reads will NOT be in
the same order as in the original file.


Other Options:
--------------

      -k=16: length of k

Change the value of the context length used. Smaller k and larger k generally
result in worse compression, but smaller k can use less resources.

      -fasta=true: If false, output seqs, one per line

Use "-fasta=false" to write out the reads without fasta headers.

      -p=10: The maximum number of threads to use

Allow kpath to use more or fewer threads.

      -flip=true: if true, reverse complement reads as needed

Use -flip=false to skip writing out the file that records which reads were
reverse complemented. 


Special options:
----------------

The following options are primarily for debugging and likely need not be used
during normal operation.

      -dups=true: if true, record dups specially
      -update=true: if true, update the reference dynamically

