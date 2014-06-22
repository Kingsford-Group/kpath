/*
   kpath - Compression of short-read sequence data
   Copyright (C) 2014  Carl Kingsford & Rob Patro

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Contact: carlk@cs.cmu.edu
*/

package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"log"
	"os"

	"kingsford/kpath/bitio"
)

/*
Read in a sorted list of kmers
Write out the kmer code using a pre-order traversal:
    1 means next character is present
    0 means next character is not present
*/

// ALPHA is the alphabet we are working over. Some code assumes it is ACGT.
const ALPHA string = "ACGT"

// children() computes the children of the given kmer "node" in the kmer list.
func children(kmers []string, start, end, depth int) [len(ALPHA)][2]int {
	var p [len(ALPHA)][2]int

	var a int
	b := start
	for ic, c := range ALPHA {
		a = b
		for b < end && kmers[b][depth] == byte(c) {
			b++
		}
		p[ic] = [...]int{a, b} // range for c is [a,b)
	}

	return p
}

// traverseToBitTree() takes a sorted list of kmers and treats it as a trie
// that it traverses in DFS-order, outputting to the bits channel a 1 whenever
// an edge exists and 0 when it does not.
func traverseToBitTree(kmers []string, bits chan<- byte) {

	count := 0

	// stack starts with all sequences on it at depth 0
	stack := make([][3]int, 0)
	start := 0
	end := len(kmers)
	depth := 0
	stack = append(stack, [...]int{start, end, depth})

	// while there are unprocessed nodes
	for len(stack) > 0 {

		// pop a node off the stack
		SED := stack[len(stack)-1]
		start, end, depth = SED[0], SED[1], SED[2]
		stack = stack[:len(stack)-1]

		if start == end {
			bits <- 0
		} else {
			bits <- 1

			// construct the node's children
			if depth < len(kmers[0]) {
				C := children(kmers, start, end, depth)
				for i := range ALPHA {
					a, b := C[i][0], C[i][1]
					stack = append(stack, [...]int{a, b, depth + 1})
				}
			} else {
				count++
			}
		}
	}
	log.Printf("Wrote %v kmers\n", count)
	if count != len(kmers) {
		panic(fmt.Errorf("Should have written %d kmers, but wrote %d!", len(kmers), count))
	}
	close(bits)
}

// decodeBitTree() reads bits from the given channel and outputs kmers on the
// output channel that were stored in the bittree. The output kmers are in no
// particular order.
func decodeBitTree(bits <-chan byte, k int, out chan<- string) {
	// stack starts with the root string
	stack := make([]string, 0)
	stack = append(stack, "")

	bitsread := 0
	for len(stack) > 0 {
		// pop node off stack
		cur := stack[len(stack)-1]
		stack = stack[:len(stack)-1]

		bitsread++
		if <-bits != 0 {
			if len(cur) == k {
				out <- cur
			} else {
				for _, c := range ALPHA {
					stack = append(stack, cur+string(c))
				}
			}
		}
	}
	log.Printf("Processed %v bits", bitsread)
	close(out)
}

// given a list of kmers, encode them to a file using the bittree scheme. The
// kmers must be sorted and they must be unique.
func encodeKmersToFile(kmers []string, out *bitio.Writer) {
	log.Printf("Encoding %v kmers to bittree file...", len(kmers))
	bits := make(chan byte, 1000000)
	go traverseToBitTree(kmers, bits)

	count := 0
	for c := range bits {
		out.WriteBit(c)
		count++
	}
	log.Printf("done. Wrote %v bits", count)
}

// readBits() creates a bit channel from a bitio.Reader().
func readBits(in *bitio.Reader, bits chan<- byte) {
	count := 0
	for {
		b, err := in.ReadBit()
		count++
		if err != nil {
			log.Printf("Stopping after %v bits", count)
			close(bits)
			return
		}
		bits <- b
	}
}

// decodeKmersFromFile() opens the given gzipped bittree file and extracts the
// stored kmers.
func decodeKmersFromFile(filename string, k int) []string {
	log.Printf("Decoding kmer buckets from %v", filename)
	// open the file and wrap a bit reader around it
	bittree, err := os.Open(filename)
	DIE_ON_ERR(err, "Couldn't open bitree file %s", filename)
	defer bittree.Close()

	bittreeZ, err := gzip.NewReader(bittree)
	DIE_ON_ERR(err, "Couldn't create gzipper")
	defer bittreeZ.Close()

	in := bitio.NewReader(bufio.NewReader(bittreeZ))
	defer in.Close()

	// start a routine to produce the bits
	bits := make(chan byte, 1000000)
	go readBits(in, bits)

	// make a channel to get the output
	out := make(chan string, 1000000)

	// decode and pass the input to the decoded output
	go decodeBitTree(bits, k, out)

	kmers := make([]string, 0)
	for s := range out {
		kmers = append(kmers, s)
	}
	log.Printf("done; found %v kmers", len(kmers))
	return kmers
}
