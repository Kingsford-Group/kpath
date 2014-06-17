package main

/* Version June 16, 2014 */

/* TODO:

4. conserve memory with a DNAString type (?)
5. add some more unit tests
*/

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strings"
    "strconv"
    "time"
    "math"

	"kingsford/arithc"
	"kingsford/bitio"
)

//===================================================================
// Kmer types
//===================================================================

// A Kmer represents a kmer of size <= 16.
type Kmer uint32

type KmerCount uint16
const MAX_OBSERVATION uint32 = (1 << 16) - 1

// A KmerInfo contains the information about a given kmer context.
type KmerInfo struct {
	next [len(ALPHA)]KmerCount
}

// A KmerHash contains the context database.
type KmerHash map[Kmer]KmerInfo

// A WeightXformFcn represents a function that can transform distribution
// counts.
type WeightXformFcn func(int, [len(ALPHA)]KmerCount) uint64

//===================================================================
// Globals
//===================================================================

var (
	encodeFlags   *flag.FlagSet
	outFile       string
	refFile       string
	readFile      string
	globalK       int
	shiftKmerMask Kmer

	defaultInterval   [len(ALPHA)]uint32 = [...]uint32{2, 2, 2, 2}
    defaultIntervalSum uint64 = 4*2

	defaultUsed   int
	contextExists int
	smoothed      int
	flipped       int
)

var (
	flipReadsOption    bool = true
    dupsOption         bool = true
    writeNsOption      bool = true
    writeFlippedOption bool = true
    updateReference    bool = true
    maxThreads         int = 10

	cpuProfile         string = ""  // set to nonempty to write profile to this file
    writeQualOption    bool = false // NYI completely
)

const (
	pseudoCount       uint64 = 1
	observationWeight uint64 = 10
	seenThreshold     KmerCount = 2 // before this threshold, increment 1 and treat as unseen
	observationInc    KmerCount = 2 // once above seenThreshold, increment by this on each observation

	//smoothOption    bool = false
)

//===================================================================
// Int <-> String Kmer representations
//===================================================================

// acgt() takes a letter and returns the index in 0,1,2,3 to which it is
// mapped. 'N's become 'A's and any other letter induces a panic.
func acgt(a byte) byte {
	switch a {
	case 'A':
		return 0
	case 'N':
		return 0
	case 'C':
		return 1
	case 'G':
		return 2
	case 'T':
		return 3
	}
	panic(fmt.Errorf("Bad character: %s!", string(a)))
}

// baseFromBits() returns the ASCII letter for the given 2-bit encoding.
func baseFromBits(a byte) byte {
	return "ACGT"[a]
}

// stringToKmer() converts a string to a 2-bit kmer representation.
func stringToKmer(kmer string) Kmer {
	var x uint64
	for _, c := range kmer {
		x = (x << 2) | uint64(acgt(byte(c)))
	}
	return Kmer(x)
}

// isACGT() returns true iff the given character is one of A,C,G, or T.
func isACGT(c rune) bool {
	return c == 'A' || c == 'C' || c == 'G' || c == 'T'
}

// kmerToString() unpacks a 2-bit encoded kmer into a string.
func kmerToString(kmer Kmer, k int) string {
	s := make([]byte, k)
	for i := 0; i < k; i++ {
		s[k-i-1] = baseFromBits(byte(kmer & 0x3))
		kmer >>= 2
	}
	return string(s)
}

// setShiftKmerMask() initializes the kmer mask. This must be called anytime
// globalK changes.
func setShiftKmerMask() {
	for i := 0; i < globalK; i++ {
		shiftKmerMask = (shiftKmerMask << 2) | 3
	}
}

// shiftKmer() creates a new kmer by shifting the given one over one base to
// the left and adding the given next character at the right.
func shiftKmer(kmer Kmer, next byte) Kmer {
	return ((kmer << 2) | Kmer(next)) & shiftKmerMask
}

// RC computes the reverse complement of a single given nucleotide. Ns become
// Ts as if they were As. Any other character induces a panic.
func RC(c byte) byte {
	switch c {
	case 'A':
		return 'T'
	case 'N':
		return 'N'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'T':
		return 'A'
	}
	panic(fmt.Errorf("Bad character: %s!", string(c)))
}

// reverseComplement() returns the reverse complement of the given string
func reverseComplement(r string) string {
	s := make([]byte, len(r))
	for i := 0; i < len(r); i++ {
		s[len(r)-i-1] = RC(r[i])
	}
	return string(s)
}

func AbsInt(x int) int {
    if x < 0 { return -x }
    return x
}

//===================================================================

// newKmerHash() constructs a new, empty kmer hash.
func newKmerHash() KmerHash {
	return make(KmerHash, 100000000)
}

// readReferenceFile() reads the sequences in the gzipped multifasta file with
// the given name and returns them as a slice of strings.
func readReferenceFile(fastaFile string) []string {
	// open the .gz fasta file that is the references
	log.Println("Reading Reference File...")
	inFasta, err := os.Open(fastaFile)
	DIE_ON_ERR(err, "Couldn't open fasta file %s", fastaFile)
	defer inFasta.Close()

	// wrap the gzip reader around it
	in, err := gzip.NewReader(inFasta)
	DIE_ON_ERR(err, "Couldn't open gzipped file %s", fastaFile)
	defer in.Close()

	out := make([]string, 0, 10000000)
	cur := make([]string, 0, 100)

	scanner := bufio.NewScanner(in)
	for scanner.Scan() {
		line := strings.TrimSpace(strings.ToUpper(scanner.Text()))
        if len(line) == 0 { continue }

		if line[0] == byte('>') {
			if len(cur) > 0 {
				out = append(out, strings.Join(cur, ""))
				cur = make([]string, 0, 100)
			}
		} else {
			cur = append(cur, line)
		}
	}
	DIE_ON_ERR(scanner.Err(), "Couldn't finish reading reference")
	return out
}

// countKmersInReference() reads the given reference file (gzipped multifasta)
// and constructs a kmer hash for it that mapps kmers to distributions of next
// characters.
func countKmersInReference(k int, fastaFile string) KmerHash {
	seqs := readReferenceFile(fastaFile)
	hash := newKmerHash()

	log.Printf("Counting %v-mer transitions in reference file...\n", k)
	for _, s := range seqs {
		if len(s) <= k {
			continue
		}
		contextMer := stringToKmer(s[:k])
		for i := 0; i < len(s)-k; i++ {
			info := hash[contextMer]
			next := acgt(s[i+k])
            //info.next[next] += observationInc XXX
            info.next[next] = observationInc
            hash[contextMer] = info

			contextMer = shiftKmer(contextMer, next)
		}
	}
	return hash
}

// capTransitionCounts() postprocesses the kmer hash to make all transition
// counts at most the given max.
func capTransitionCounts(hash KmerHash, max int) {
	M := KmerCount(max)
	for _, v := range hash {
		for i, count := range v.next {
			if count > M {
				v.next[i] = M
			}
		}
	}
}

//===================================================================
// Encoding
//===================================================================

// contextWeight() is a weight transformation function that will change the
// distribution weights according to the function for real contexts. If the
// count is too small, it returns the pseudocount; if the count is big enough
// it returns observationWeight * the distribution value.
func contextWeight(charIdx int, dist [len(ALPHA)]KmerCount) (w uint64) {
	if dist[charIdx] >= seenThreshold {
		w = observationWeight * uint64(dist[charIdx]) / uint64(observationInc)
		return
	}
	w = pseudoCount
	return
}

// defaultWeight() is a weight transformation function for the default
// distribution. It returns the weight unchanged.
func defaultWeight(charIdx int, dist [len(ALPHA)]KmerCount) uint64 {
	return uint64(dist[charIdx])
}

// intervalFor() returns the interval for the given character (represented as a
// 2-bit encoded base) according to the given distribution (transformed by the
// given weight transformation function).
func intervalFor(
        letter byte, 
        dist [len(ALPHA)]KmerCount, 
        weightOf WeightXformFcn,
) (a uint64, b uint64, total uint64) {

	letterIdx := int(letter)
	for i := 0; i < len(dist); i++ {
		w := weightOf(i, dist)

		total += w
		if i <= letterIdx {
			b += w
			if i < letterIdx {
				a += w
			}
		}
	}
	return
}

// intervalForDefault() computes the interval for the given character using the
// default interval
func intervalForDefault(letter byte) (a uint64, b uint64, total uint64) {
    letterIdx := int(letter)
    for i := 0; i < len(defaultInterval); i++ {
        w := uint64(defaultInterval[i])
        total += w
        if i <= letterIdx {
            b += w
            if i < letterIdx {
                a += w
            }
        }
    }
    return
}
    

// nextInterval() computes the interval for the given context and updates the
// default distribution and context distributions as required.
func nextInterval(
    hash KmerHash, 
    contextMer Kmer, 
    kidx byte,
) (a uint64, b uint64, total uint64) {
	info, ok := hash[contextMer]
	// if the context exists, use that distribution
	if ok {
		contextExists++
		a, b, total = intervalFor(kidx, info.next, contextWeight)
        if updateReference {
            prevVal := info.next[kidx]
            if prevVal >= seenThreshold { 
                if uint32(prevVal) + uint32(observationInc) < MAX_OBSERVATION {
                    info.next[kidx] += observationInc
                }
            } else {
                info.next[kidx]++
            }
            hash[contextMer] = info
        }
	} else {
		// if the context doesnt exist, use a simple default interval
		defaultUsed++
		//a, b, total = intervalFor(kidx, defaultInterval, defaultWeight)
        a, b, total = intervalForDefault(kidx) //XXX
		defaultInterval[kidx]++
        defaultIntervalSum++

        if updateReference {
            // add this to the context now
            info.next[kidx]++
            hash[contextMer] = info
        }
	}
	return
}

// countMatchingObservations() counts the number of observaions of kmers in the
// read.
func countMatchingObservations(hash KmerHash, r string) (n KmerCount) {
	contextMer := stringToKmer(r[:globalK])
	for i := globalK; i < len(r); i++ {
		symb := acgt(r[i])
		if H, ok := hash[contextMer]; ok {
			n += H.next[symb]
		}
		contextMer = shiftKmer(contextMer, symb)
	}
	return
}

// support sorting the fastq list lexicographically
type Lexicographically []*FastQ
func (a Lexicographically) Len() int { return len(a) }
func (a Lexicographically) Swap(i, j int) { a[i], a[j] = a[j], a[i] }

func (a Lexicographically) Less(i, j int) bool { 
    for i, c := range a[i].Seq[:globalK] {
        d := a[j].Seq[i]
        if c < d { return true }
        if c > d { return false }
    }
    return false
}
        
// flipRange() flips the reads in the given slice if the reverse complement
// matches the reference better.
func flipRange(block []*FastQ, hash KmerHash) int {
    flip := 0
    for _, fq := range block {
        n1 := countMatchingObservations(hash, string(fq.Seq))
        rcr := reverseComplement(string(fq.Seq))
        n2 := countMatchingObservations(hash, rcr)

        // if they are tied, take the lexigographically smaller one
        if n2 > n1 || (n2 == n1 && string(rcr) < string(fq.Seq)) {
            fq.SetReverseComplement(rcr)
            flip++
        }
    }
    return flip
}

// readAndFlipReads() reads the reads and reverse complements them if the
// reverse complement matches the hash better (according to a countMatching*
// function above). It returns a slice of the reads. "N"s are treated as "A"s.
// No other characters are transformed and will eventually lead to a panic.
func readAndFlipReads(
    readFile string, 
    hash KmerHash, 
    flipReadsOption bool,
) []*FastQ {
    // read the reads from the file into memory
    log.Printf("Reading reads...")
    readStart := time.Now()
    fq := make(chan*FastQ, 10000000)
    go ReadFastQ(readFile, fq)
    reads := make([]*FastQ, 0, 10000000)
    for rec := range fq {
        reads = append(reads, rec)
    }
    readEnd := time.Now()
    log.Printf("Time: read %v reads; spent %v seconds.", len(reads), readEnd.Sub(readStart).Seconds())

    // if enabled, start several threads to flip the reads
    if flipReadsOption {
        // start maxThreads-1 workers to flip the read ranges
        wait := make([]chan int, maxThreads-1)
        for i := range wait {
            wait[i] = make(chan int)
        }
        blockSize := 1 + len(reads) / len(wait)
        log.Printf("Have %v read flippers, each working on %v reads", len(wait), blockSize)
        for i, c := range wait {
            go func(i int, c chan int) {
                end := (i+1)*blockSize
                if end > len(reads) { end = len(reads) }
                log.Printf("Worker %v flipping [%d, %d)...", i, i*blockSize, end)
                count := flipRange(reads[i*blockSize : end], hash)
                c <- count
                close(c)
                log.Printf("Worker %v flipped %d reads", i, count)
                runtime.Goexit()
                return
            }(i, c)
        }

        // wait for all the workers to finish and sum up their 
        for _, c := range wait { 
            for f := range c {
                flipped += f
            }
        }
    }
    flipEnd := time.Now()
    log.Printf("Time: flipping: %v seconds.", flipEnd.Sub(readEnd).Seconds())

    // sort the records by sequence
    sort.Sort(Lexicographically(reads))
    readSort := time.Now()
    log.Printf("Time: sorting reads: %v seconds.", readSort.Sub(flipEnd).Seconds())

	log.Printf("Read %v reads; flipped %v of them.", len(reads), flipped)
    return reads

}
/*
// readAndFlipReads() reads the reads and reverse complements them if the
// reverse complement matches the hash better (according to a countMatching*
// function above). It returns a slide of the reads. "N"s are treated as "A"s.
// No other characters are transformed and will eventually lead to a panic.
func readAndFlipReads(readFile string, hash KmerHash, flipReadsOption bool) []*FastQ {
    // start the reading routine
    log.Printf("Reading reads...")
    readStart := time.Now()
    fq := make(chan *FastQ, 10000000)
    go ReadFastQ(readFile, fq)

    reads := make([]*FastQ, 0, 10000000)

    // for every record
    for rec := range fq {
        // possibly flip it
        if flipReadsOption {
            n1 := countMatchingObservations(hash, string(rec.Seq))
            rcr := reverseComplement(string(rec.Seq))
            n2 := countMatchingObservations(hash, rcr)

			// if they are tied, take the lexigographically smaller one
			if n2 > n1 || (n2 == n1 && string(rcr) < string(rec.Seq)) {
				rec.SetReverseComplement(rcr)
				flipped++
			}
        }
        // save it in our read list
        reads = append(reads, rec)
    }
    readEnd := time.Now()

    // sort the records by sequence
    sort.Sort(Lexicographically(reads))
    readSort := time.Now()

	log.Printf("Read %v reads; flipped %v of them.", len(reads), flipped)
    log.Printf("Time: reading and flipping: %v seconds.", readEnd.Sub(readStart).Seconds())
    log.Printf("Time: sorting reads: %v seconds.", readSort.Sub(readEnd).Seconds())
    return reads
}
*/


// listBuckets() processes the reads and creates the bucket list and the list
// of the bucket sizes and returns them.
func listBuckets(reads []*FastQ) ([]string, []int) {
	curBucket := ""
    prevRead := ""
    allSame := false
	buckets := make([]string, 0, 1000000)
	counts := make([]int, 0, 1000000)

	for _, rec := range reads {
        r := string(rec.Seq)
		if r[:globalK] != curBucket {
            // if all the reads in a bucket are the same, record this
            // by negating the bucket count
            if dupsOption && allSame && counts[len(counts)-1] > 1 {
                counts[len(counts)-1] = -counts[len(counts)-1]
            }

			curBucket = r[:globalK]
            prevRead = r
			buckets = append(buckets, curBucket)
			counts = append(counts, 1)
            allSame = true
		} else {
            allSame = allSame && (r == prevRead)
            prevRead = r
			counts[len(counts)-1]++
		}
	}
	return buckets, counts
}

// writeCounts() writes the counts list out to the given writer.
func writeCounts(f io.Writer, readlen int, counts []int) {
	log.Printf("Writing counts...")
	fmt.Fprintf(f, "%d ", readlen)
	for _, c := range counts {
		fmt.Fprintf(f, "%d ", c)
	}
	log.Printf("Done; write %d counts.", len(counts))
}


// writeNLocations() writes out the locations of the translated Ns in the file.
func writeNLocations(f io.Writer, reads []*FastQ) {
    log.Printf("Writing location of Ns...")
    // every read's locations are written as a space separated list of ascii
    // integers
    c := 0
    for _, fq := range reads {
        for i, p := range fq.NLocations {
            fmt.Fprintf(f, "%d", p)
            c++
            if i != len(fq.NLocations)-1 {
                fmt.Fprintf(f, " ")
            }
        }
        fmt.Fprintf(f, "\n")
    }
    log.Printf("Done; wrote %d Ns.", c)
}

// writeFlipped() writes out a stream of bits that says whether or not the
// reads were flipped.
func writeFlipped(out *bitio.Writer, reads []*FastQ) {
    for _, fq := range reads {
        if fq.IsFlipped {
            out.WriteBit(1)
        } else {
            out.WriteBit(0)
        }
    }
}

// encodeSingleReadWithBucket() encodes a single read: uses a bucketing scheme
// for initial part, and arithmetic encoding for the rest.
func encodeSingleReadWithBucket(contextMer Kmer, r string, hash KmerHash, coder *arithc.Encoder) {
	// encode rest using the reference probs
	for i := globalK; i < len(r); i++ {
		char := acgt(r[i])
		a, b, total := nextInterval(hash, contextMer, char)
		err := coder.Encode(a, b, total)
		DIE_ON_ERR(err, "Error encoding read: %s", r)
		contextMer = shiftKmer(contextMer, byte(char))
	}
}


// encodeWithBuckets() reads the reads, creates the buckets, saves the buckets
// and their counts, and then encodes each read.
func encodeWithBuckets(
    readFile, 
    outBaseName string, 
    hash KmerHash, 
    coder *arithc.Encoder,
) (n int) {
	// read the reads and flip as needed
	reads := readAndFlipReads(readFile, hash, flipReadsOption)

    readLength := len(reads[0].Seq)

    log.Printf("Estimated 2-bit encoding size: %d", 
        uint64(math.Ceil(float64(2*len(reads)*readLength) / 8.0)))

    // if the user wants the qualities written out 
    waitForFlipped := make(chan struct{})
    if writeFlippedOption {
        outFlipped, err := os.Create(outBaseName + ".flipped")
        DIE_ON_ERR(err, "Couldn't create flipped file: %s", outBaseName + ".flipped")
        defer outFlipped.Close()

        outFlippedZ, err := gzip.NewWriterLevel(outFlipped, gzip.BestCompression)
        DIE_ON_ERR(err, "Couldn't create gzipper for flipped file.")
        defer outFlippedZ.Close()

        flippedBits := bitio.NewWriter(outFlippedZ)
        defer flippedBits.Close()

        go func() {
            writeFlipped(flippedBits, reads)
            close(waitForFlipped)
            runtime.Goexit()
            return
        }()
    } else {
        close(waitForFlipped)
    }

    // if the user wants to write out the N positions, write them out
    waitForNs := make(chan struct{})
    if writeNsOption {
        outNs, err := os.Create(outBaseName + ".ns")
        DIE_ON_ERR(err, "Couldn't create N location file: %s", outBaseName + ".ns")
        defer outNs.Close()

        outNsZ, err := gzip.NewWriterLevel(outNs, gzip.BestCompression)
        DIE_ON_ERR(err, "Couldn't create gzipper for N location file.")
        defer outNsZ.Close()

        go func() {
            writeNLocations(outNsZ, reads)
            close(waitForNs)
            runtime.Goexit()
            return
        }()
    } else {
        close(waitForNs)
    }
     
	// create the buckets and counts
	buckets, counts := listBuckets(reads)

	// write the bittree for the bucket out to a file
	outBT, err := os.Create(outBaseName + ".bittree")
	DIE_ON_ERR(err, "Couldn't create bucket file: %s", outBaseName+".bittree")
	defer outBT.Close()

	// compress the file with gzip as we are writing it
	outBZ, err := gzip.NewWriterLevel(outBT, gzip.BestCompression)
	DIE_ON_ERR(err, "Couldn't create gzipper for bucket file")
	defer outBZ.Close()

	// create a writer that lets us write bits
	writer := bitio.NewWriter(outBZ)
	defer writer.Close()

	/*** The main work to encode the bucket names ***/
	waitForBuckets := make(chan struct{})
	go func() {
		encodeKmersToFile(buckets, writer)
		close(waitForBuckets)
        runtime.Goexit()
        return
	}()

	// write out the counts
	countF, err := os.Create(outBaseName + ".counts")
	DIE_ON_ERR(err, "Couldn't create counts file: %s", outBaseName+".counts")
	defer countF.Close()

	// compress it as we are writing it
	countZ, err := gzip.NewWriterLevel(countF, gzip.BestCompression)
	DIE_ON_ERR(err, "Couldn't create gzipper for count file")
	defer countZ.Close()

	/*** The main work to encode the bucket counts ***/
	waitForCounts := make(chan struct{})
	go func() {
		writeCounts(countZ, readLength, counts)
		close(waitForCounts)
        runtime.Goexit()
        return
	}()
	// Wait for each of the coders to finish
	<-waitForBuckets
	<-waitForCounts
    <-waitForNs
    <-waitForFlipped

    fmt.Printf("Currently have %v Go routines...", runtime.NumGoroutine())
	/*** The main work to encode the read tails ***/
    encodeStart := time.Now()
	log.Printf("Encoding reads, each of length %d ...", readLength)
	//waitForReads := make(chan struct{})

    curRead := 0
    for i, c := range counts {
        bucketMer := stringToKmer(buckets[i])
        if c > 0 {
            // write out the given number of reads
            for j := 0; j < c; j++ {
                encodeSingleReadWithBucket(bucketMer, string(reads[curRead].Seq), hash, coder)
                curRead++
                n++
            }
        } else {
            // all the reads in this bucket are the same, so just write one
            // and skip past the rest.
            encodeSingleReadWithBucket(bucketMer, string(reads[curRead].Seq), hash, coder)
            curRead += AbsInt(c)
            n++
        }
    }

	//<-waitForReads
	log.Printf("done. Took %v seconds to encode the tails.", time.Now().Sub(encodeStart).Seconds())
	return
}

//===============================================================================
// DECODING
//===============================================================================

// readBucketCounts() opens the file with the given name and parses it to
// extract a list of bucket sizes that were written by the encoding. The given
// file must have been written by the coder --- it is assumed to be a gzipped
// list of space-separated ASCII numbers.
func readBucketCounts(countsFN string) ([]int, int) {
	log.Printf("Reading bucket counts from %v", countsFN)

	// open the count file
	c1, err := os.Open(countsFN)
	DIE_ON_ERR(err, "Couldn't open count file: %s", countsFN)
	defer c1.Close()

	// the count file is compressed with gzip; uncompress it as we read it
	c, err := gzip.NewReader(c1)
	DIE_ON_ERR(err, "Couldn't create gzip reader: %v")
	defer c.Close()

	var n, readlen int
	_, err = fmt.Fscanf(c, "%d", &readlen)
	DIE_ON_ERR(err, "Couldn't read read length from counts file")

	counts := make([]int, 0)
	err = nil
	for x := 1; err == nil && x > 0; {
		x, err = fmt.Fscanf(c, "%d", &n)
		if x > 0 && err == nil {
			counts = append(counts, n)
		}
	}
	log.Printf("done; read %d counts", len(counts))
	return counts, readlen
}

// readFlipped() reads the compressed bitstream that indicates whether a read
// was flipped or not. If the file does not exist, returns nil.
func readFlipped(flippedFN string) []bool {
    // open the file; return empty if nothing there
    flippedIn, err := os.Open(flippedFN)
    if err == nil {
        log.Printf("Reading flipped bits from %s", flippedFN)
        defer flippedIn.Close()

        flippedZ, err := gzip.NewReader(flippedIn)
        DIE_ON_ERR(err, "Couldn't create unzipper for flipped file")
        defer flippedZ.Close()

        flippedBits := bitio.NewReader(bufio.NewReader(flippedZ))
        defer flippedBits.Close()
        
        flipped := make([]bool, 0, 1000000)
        for {
            b, err := flippedBits.ReadBit()
            if err != nil { break }
            if b > 0 {
                flipped = append(flipped, true)
            } else {
                flipped = append(flipped, false)
            }
        }
        log.Printf("Read %d bits indicating whether reads were flipped.", len(flipped))
        return flipped
    } else {
        log.Printf("No flipped bit file (%s) found; ignoring.", flippedFN)
        return nil
    }
}

// readNLocations() reads the compressed N location file and returns a slice of
// slices that contain the positions of the Ns. An optimization is made that if
// there are no Ns in a read, then out[r] will be nil rather than an empty
// list.  If the file is not found, will return nil
func readNLocations(nLocFN string) [][]byte {
    // open the file; return empty if nothing there
    inNs, err := os.Open(nLocFN)
    if err == nil {
        log.Printf("Reading locations of Ns from %s", nLocFN)
        defer inNs.Close()
        inZ, err := gzip.NewReader(inNs)
        DIE_ON_ERR(err, "Couldn't create gzipper for N locations")
        defer inZ.Close()

        locs := make([][]byte, 0, 10000000)
        ncount := 0
        
        // for every line in the input file
        scanner := bufio.NewScanner(inZ)
        for scanner.Scan() {
            // split into the list of integers (as strings)
            posns := strings.Split(strings.TrimSpace(scanner.Text()), " ")

            // if there are any Ns in this read
            if len(posns) > 0  && posns[0] != "" {
                // create a new slice to hold them, and convert them to integers
                locs = append(locs, make([]byte, 0))
                for _, v := range posns {
                    p, err := strconv.Atoi(v)
                    DIE_ON_ERR(err, "Badly formatted N location file!")
                    locs[len(locs)-1] = append(locs[len(locs)-1], byte(p))
                }
                ncount += len(posns)
            } else {
                // otherwise, for reads with no Ns, the slice is just nil
                locs = append(locs, nil)
            }
        }
        DIE_ON_ERR(scanner.Err(), "Couldn't finish reading N locations")
        log.Printf("Read locations for %d Ns.", ncount)
        return locs
    } else {
        log.Printf("No file with N locations (%s) was found; ignoring.", nLocFN)
        return nil
    }
}

// dart() finds the interval in the given distribution that contains the given
// target, after transformming the distribution using the given weightOf
// function. This is called by lookup() during decode.
func dart(
    dist [len(ALPHA)]KmerCount, 
    target uint32, 
    weightOf WeightXformFcn,
) (uint64, uint64, uint64) {
	sum := uint32(0)
	for i := range dist {
		w := uint32(weightOf(i, dist))
		sum += w
		if target < sum {
			return uint64(sum - w), uint64(sum), uint64(i)
		}
	}
	panic(fmt.Errorf("Couldn't find range for target %d", target))
}

func dartDefault(target uint64) (uint64, uint64, uint64) {
    sum := uint64(0)
    for i, w := range defaultInterval {
        sum += uint64(w)
        if target < sum {
            return sum - uint64(w), sum, uint64(i)
        }
    }
	panic(fmt.Errorf("Couldn't find range for target %d", target))
}

// lookup() is called by arithc.Decoder to find an interval that contains the
// given value t.
func lookup(hash KmerHash, context Kmer, t uint64) (uint64, uint64, uint64) {
	if info, ok := hash[context]; ok {
		return dart(info.next, uint32(t), contextWeight)
	} else {
		//return dart(defaultInterval, uint32(t), defaultWeight)
		return dartDefault(t)
	}
}

// sumDist() computes the sum of the items in the given distribution after
// first transforming them via the given weightOf function.
func sumDist(d [len(ALPHA)]KmerCount, weightOf WeightXformFcn) (total uint64) {
	for i := range d {
		total += uint64(weightOf(i, d))
	}
	return
}

// contextTotal() returns the total sum of the appropriate distribution: the
// distribution of the given context (if found) or the default distribution
// (otherwise).
func contextTotal(hash KmerHash, context Kmer) uint64 {
	if info, ok := hash[context]; ok {
		return sumDist(info.next, contextWeight)
	} else {
        return defaultIntervalSum
		//return sumDist(defaultInterval, defaultWeight)
	}
}

// decodeSingleRead() does the work of decoding a single read. 
func decodeSingleRead(
    contextMer Kmer, 
    hash KmerHash, 
    tailLen int, 
    decoder *arithc.Decoder, 
    out []byte,
) {
    // function called by Decode
	lu := func(t uint64) (uint64, uint64, uint64) {
		return lookup(hash, contextMer, t)
	}

    for i := 0; i < tailLen; i++ {
        // decode next symbol
        symb, err := decoder.Decode(contextTotal(hash, contextMer), lu)
        DIE_ON_ERR(err, "Fatal error decoding!")
        b := byte(symb)

        // write it out
        out[i] = baseFromBits(b)

        // update hash counts (throws away the computed interval; just
        // called for side effects.)
        nextInterval(hash, contextMer, b)

        // update the new context
        contextMer = shiftKmer(contextMer, b)
    }
}

// replace the letters at the given position by Ns
func putbackNs(s string, p []byte) string {
    b := []byte(s)
    for _, v := range p {
        b[v] = 'N'
    }
    return string(b)
}


// decodeReads() decodes the file wrapped by the given Decoder, using the
// kmers, counts, and hash table provided. It writes its output to the given
// io.Writer.
func decodeReads(
    kmers []string, 
    counts []int, 
    isFlipped []bool,
    nLocations [][]byte,
    hash KmerHash, 
    readLen int, 
    out io.Writer, 
    decoder *arithc.Decoder,
) {
	log.Printf("Decoding reads...")

	n := 0
    ncount := 0
	buf := bufio.NewWriter(out)

    patchAndWriteRead := func(head, tail string) {
        // put the head & tail together
        s := fmt.Sprintf("%s%s", head, tail)
        // put back the ns if we have them
        if nLocations != nil {
            s = putbackNs(s, nLocations[n])
            ncount += len(nLocations[n])
        }
        // unflip the reads if we have them
        if isFlipped != nil && isFlipped[n] {
            s = reverseComplement(s)
            flipped++
        }
        // write it out
        buf.Write([]byte(s))
        buf.WriteByte('\n')
        return
    }

    // tailBuf is a buffer for read tails returned by decodeSingleRead
    tailLen := readLen-len(kmers[0])
    tailBuf := make([]byte, tailLen)

    // for every bucket
    for curBucket, c := range counts {
        contextMer := stringToKmer(kmers[curBucket])

        // if bucket is a uniform bucket, write out |c| copies of the decoded
        // string
        if c < 0 {
            decodeSingleRead(contextMer, hash, tailLen, decoder, tailBuf)
            for j := 0; j < AbsInt(c); j++ {
                patchAndWriteRead(kmers[curBucket], string(tailBuf))
                n++
            }
        } else {
            // otherwise, decode a read for each string in the bucket
            for j := 0; j < c; j++ {
                decodeSingleRead(contextMer, hash, tailLen, decoder, tailBuf)
                patchAndWriteRead(kmers[curBucket], string(tailBuf))
                n++
            }
        }
    }
	buf.Flush()
	log.Printf("Added back %d Ns to the reads.", ncount)
	log.Printf("done. Wrote %v reads; %d were flipped", n, flipped)
}

//===================================================================
// Command line and main driver
//===================================================================

// DIE_ON_ERR() logs a fatal error to the standard logger if err != nil and
// exits the program. It also prints the given informative message.
func DIE_ON_ERR(err error, msg string, args ...interface{}) {
	if err != nil {
		log.Printf("Error: "+msg, args...)
		log.Fatalln(err)
	}
}

// init() is called automatically on program start up. Here, it creates the
// command line parser.
func init() {
	encodeFlags = flag.NewFlagSet("encode", flag.ContinueOnError)
	encodeFlags.StringVar(&refFile, "ref", "", "reference fasta filename")
	encodeFlags.StringVar(&outFile, "out", "", "output filename")
	encodeFlags.StringVar(&readFile, "reads", "", "reads filename")
	encodeFlags.IntVar(&globalK, "k", 16, "length of k")
    encodeFlags.BoolVar(&flipReadsOption, "flip", true, "if true, reverse complement reads as needed") 
    encodeFlags.BoolVar(&dupsOption, "dups", true, "if true, record dups specially")
    encodeFlags.BoolVar(&updateReference, "update", true, "if true, update the reference dynamically")

    encodeFlags.StringVar(&cpuProfile, "cpuProfile", "", "if nonempty, write pprof profile to given file.")
}

// writeGlobalOptions() writes out the global variables that can affect the
// encoding / decoding. Files encoded with one set of options can only be
// decoded using the same set of options.
func writeGlobalOptions() {
	log.Printf("Option: psudeoCount = %d", pseudoCount)
	log.Printf("Option: observationWeight = %d", observationWeight)
	log.Printf("Option: seenThreshold = %d", seenThreshold)
	log.Printf("Option: observationInc = %d", observationInc)
	//log.Printf("Option: smoothOption = %v", smoothOption)
	log.Printf("Option: flipReadsOption = %v", flipReadsOption)
    log.Printf("Option: dupsOption = %v", dupsOption)
    log.Printf("Option: updateReference = %v", updateReference)
}

// main() encodes or decodes a set of reads based on the first command line
// argument (which is either encode or decode).
func main() {
	log.Println("Starting kpath version 6-16-14")

    startTime := time.Now()

	log.Printf("Maximum threads = %v", maxThreads)
	runtime.GOMAXPROCS(maxThreads)

	// parse the command line
	const (
		ENCODE int = 1
		DECODE int = 2
	)
	if len(os.Args) < 2 {
		encodeFlags.PrintDefaults()
		os.Exit(1)
	}
	var mode int
	if os.Args[1][0] == 'e' {
		mode = ENCODE
		log.SetPrefix("kpath (encode): ")
	} else {
		mode = DECODE
		log.SetPrefix("kpath (decode): ")
	}
	encodeFlags.Parse(os.Args[2:])
	if globalK <= 0 || globalK > 16 {
		log.Fatalf("K must be specified as a small positive integer with -k")
	}
	log.Printf("Using kmer size = %d", globalK)
	setShiftKmerMask()

	if cpuProfile != "" {
		log.Printf("Writing CPU profile to %s", cpuProfile)
		cpuF, err := os.Create(cpuProfile)
		DIE_ON_ERR(err, "Couldn't create CPU profile file %s", cpuProfile)
		pprof.StartCPUProfile(cpuF)
		defer pprof.StopCPUProfile()
	}

	// count the kmers in the reference
	var hash KmerHash
	waitForReference := make(chan struct{})
	go func() {
        refStart := time.Now()
		hash = countKmersInReference(globalK, refFile)
		log.Printf("There are %v unique %v-mers in the reference\n",
			len(hash), globalK)
        log.Printf("Time: Took %v seconds to read reference.", time.Now().Sub(refStart).Seconds())
		close(waitForReference)
        return
	}()

	writeGlobalOptions()

	if mode == ENCODE {
		/* encode -k -ref -reads=FOO.seq -out=OUT
		   will encode into OUT.{enc,bittree,counts} */
		log.Printf("Reading from %s", readFile)
		log.Printf("Writing to %s, %s, %s", 
            outFile+".enc", outFile+".bittree", outFile+".counts")

        /*
		var err error
		if smoothOption {
			smoothFile, err = os.Create("smoothed.txt")
			DIE_ON_ERR(err, "Couldn't create smoothed file 'smoothed.txt'")
			defer smoothFile.Close()
		}
        */

		// create the output file
		outF, err := os.Create(outFile + ".enc")
		DIE_ON_ERR(err, "Couldn't create output file %s", outFile)
		defer outF.Close()

        //outBuf := bufio.NewWriter(outF)
        //defer outBuf.Flush()

		writer := bitio.NewWriter(outF)
		defer writer.Close()

		// create encoder
		encoder := arithc.NewEncoder(writer)
		defer encoder.Finish()

		// encode reads
		<-waitForReference
		n := encodeWithBuckets(readFile, outFile, hash, encoder)
		//log.Printf("Smoothed (changed) %v characters", smoothed)
		log.Printf("Reads Flipped: %v", flipped)
		log.Printf("Encoded %v reads.", n)

	} else {
		/* decode -k -ref -reads=FOO -out=OUT.seq
		   will look for FOO.enc, FOO.bittree, FOO.counts and decode into OUT.seq */

		tailsFN := readFile + ".enc"
		headsFN := readFile + ".bittree"
		countsFN := readFile + ".counts"

		log.Printf("Reading from %s, %s, and %s", tailsFN, headsFN, countsFN)

		// read the bucket names
		var kmers []string
		waitForBuckets := make(chan struct{})
		go func() {
			kmers = decodeKmersFromFile(headsFN, globalK)
			sort.Strings(kmers)
			close(waitForBuckets)
            runtime.Goexit()
            return
		}()

		// read the bucket counts
		var counts []int
		var readlen int
		waitForCounts := make(chan struct{})
		go func() {
			counts, readlen = readBucketCounts(countsFN)
			close(waitForCounts)
            runtime.Goexit()
            return
		}()

        // read the flipped bits --- flipped by be 0-length if no file could be
        // found; this indicates that either nothing was flipped or we don't
        // care about orientation
        var flipped []bool
        waitForFlipped := make(chan struct{})
        go func() {
            flipped = readFlipped(readFile + ".flipped")
            close(waitForFlipped)
            runtime.Goexit()
            return
        }()

        // read the NLocations, which might be 0-length if no file could be
        // found; this indicates that the Ns were recorded some other way.
        var NLocations [][]byte
        waitForNLocations := make(chan struct{})
        go func() {
            NLocations = readNLocations(readFile + ".ns")
            close(waitForNLocations)
            runtime.Goexit()
            return
        }()

		// open encoded read file
		encIn, err := os.Open(tailsFN)
		DIE_ON_ERR(err, "Can't open encoded read file %s", tailsFN)
		defer encIn.Close()

		// create a bit reader wrapper around it
		reader := bitio.NewReader(bufio.NewReader(encIn))
		defer reader.Close()

		// create a decoder around it
		decoder, err := arithc.NewDecoder(reader)
		DIE_ON_ERR(err, "Couldn't create decoder!")

		// create the output file
		log.Printf("Writing to %s", outFile)
		outF, err := os.Create(outFile)
		DIE_ON_ERR(err, "Couldn't create output file %s", outFile)
		defer outF.Close()

		<-waitForReference
		<-waitForBuckets
		<-waitForCounts
        <-waitForFlipped
        <-waitForNLocations
		log.Printf("Read length = %d", readlen)
		decodeReads(kmers, counts, flipped, NLocations, hash, readlen, outF, decoder)
	}
	log.Printf("Default interval used %v times and context used %v times",
		defaultUsed, contextExists)

    endTime := time.Now()
    log.Printf("kpath took %v to run.", endTime.Sub(startTime).Seconds())

    /* UNCOMMENT TO DEBUG GARBAGE COLLECTION WITH GO 1.2
    var stats debug.GCStats
    stats.PauseQuantiles = make([]time.Duration, 5)
    debug.ReadGCStats(&stats)
    log.Printf("Last GC=%v\nNum GC=%v\nPause for GC=%v\nPauseHistory=%v",
        stats.LastGC, stats.NumGC, stats.PauseTotal.Seconds(), stats.Pause)
    */
}

