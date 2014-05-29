package main

/* TODO:
x0. change pseudocount to be max (or else change the seen threshold)
x1. fix filename handling (including outputting filename information as a log)
x why isn't TTTT...TTT flipped? b/c it's tied
x2. compress counts and bittree directly from here (and uncompress bittree when reading it)
x3. move unused code over to unused.no file
x4. update error messages and panic messages to be more consistent
    // DIE_ON_ERROR(err, "Couldn't create bucket file: %v", err)
x4.1 write out all the global options when encoding / decoding
x5.2. add comments
x6. profile to speed up
x7. parallelize
x5.0. read / write READLEN someplace

5.1. test out on 3 more files

5. conserve memory with a DNAString type (?)
8. add some more unit tests
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

	"kingsford/arithc"
	"kingsford/bitio"
)

//===================================================================
// Kmer types
//===================================================================

// A Kmer represents a kmer of size < 32.
type Kmer uint32

// A KmerInfo contains the information about a given kmer context.
type KmerInfo struct {
	next [len(ALPHA)]uint32
}

// A KmerHash contains the context database.
type KmerHash map[Kmer]*KmerInfo

// A WeightXformFcn represents a function that can transform distribution
// counts.
type WeightXformFcn func(int, [len(ALPHA)]uint32) uint64

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
	readStartInterval [len(ALPHA)]uint32 = [...]uint32{2, 2, 2, 2}

	defaultUsed   int
	contextExists int
	smoothed      int
	flipped       int
)

const (
	pseudoCount       uint64 = 1
	observationWeight uint64 = 10
	seenThreshold     uint32 = 2 // before this threshold, increment 1 and treat as unseen
	observationInc    uint32 = 2 // once above seenThreshold, increment by this on each observation

	smoothOption    bool = false
	flipReadsOption bool = true

	cpuProfile string = "" //"encode.pprof" // set to nonempty to write profile to this file
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

// shiftKmer() creates a new kmer by shifting the given one over one base to the left
// and adding the given next character at the right.
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
		return 'T'
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
		if len(line) > 0 && line[0] == '>' {
			if len(cur) > 0 {
				out = append(out, strings.Join(cur, ""))
				cur = make([]string, 0, 100)
			}
		} else if len(line) > 0 {
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
			if _, ok := hash[contextMer]; !ok {
				hash[contextMer] = &KmerInfo{}
			}
			next := acgt(s[i+k])
			hash[contextMer].next[next] += observationInc
			contextMer = shiftKmer(contextMer, next)
		}
	}
	return hash
}

// capTransitionCounts() postprocesses the kmer hash to make all transition
// counts at most the given max.
func capTransitionCounts(hash KmerHash, max int) {
	M := uint32(max)
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
func contextWeight(charIdx int, dist [len(ALPHA)]uint32) (w uint64) {
	if dist[charIdx] >= seenThreshold {
		w = observationWeight * uint64(dist[charIdx]) / uint64(observationInc)
		return
	}
	w = pseudoCount
	return
}

// defaultWeight() is a weight transformation function for the default
// distribution. It returns the weight unchanged.
func defaultWeight(charIdx int, dist [len(ALPHA)]uint32) uint64 {
	return uint64(dist[charIdx])
}

// intervalFor() returns the interval for the given character (represented as a
// 2-bit encoded base) according to the given distribution (transformed by the
// given weight transformation function).
func intervalFor(letter byte, dist [len(ALPHA)]uint32, weightOf WeightXformFcn) (a uint64, b uint64, total uint64) {
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

// nextInterval() computes the interval for the given context and updates the
// default distribution and context distributions as required.
func nextInterval(hash KmerHash, contextMer Kmer, kidx byte) (a uint64, b uint64, total uint64) {
	info, ok := hash[contextMer]
	// if the context exists, use that distribution
	if ok {
		contextExists++
		a, b, total = intervalFor(kidx, info.next, contextWeight)
		if info.next[kidx] >= seenThreshold { // increment double if in the transcriptome
			info.next[kidx] += observationInc
		} else {
			info.next[kidx]++
		}
	} else {
		// if the context doesnt exist, use a simple default interval
		defaultUsed++
		a, b, total = intervalFor(kidx, defaultInterval, defaultWeight)
		defaultInterval[kidx]++

		// add this to the context now
		hash[contextMer] = &KmerInfo{}
		hash[contextMer].next[kidx]++
	}
	return
}

// countMatchingObservations() counts the number of observaions of kmers in the
// read.
func countMatchingObservations(hash KmerHash, r string) (n uint32) {
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

// readAndFlipReads() reads the reads and reverse complements them if the
// reverse complement matches the hash better (according to a countMatching*
// function above). It returns a slide of the reads. "N"s are treated as "A"s.
// No other characters are transformed and will eventually lead to a panic.
func readAndFlipReads(readFile string, hash KmerHash, flipReadsOption bool) []string {
	// open the read file
	log.Println("Reading and flipping reads...")
	in, err := os.Open(readFile)
	DIE_ON_ERR(err, "Couldn't open read file %s", readFile)
	defer in.Close()

	// put the reads into a global array, flipped if needed
	reads := make([]string, 0, 1000000)
	scanner := bufio.NewScanner(in)
	for scanner.Scan() {
		// remove spaces and convert on-ACGT to 'A'
		r := strings.Replace(strings.TrimSpace(strings.ToUpper(scanner.Text())), "N", "A", -1)
		if flipReadsOption {
			n1 := countMatchingObservations(hash, r)
			rcr := reverseComplement(r)
			n2 := countMatchingObservations(hash, rcr)
			// if they are tied, take the lexigographically smaller one
			if n2 > n1 || (n2 == n1 && rcr < r) {
				r = rcr
				flipped++
			}
		}
		reads = append(reads, r)
	}
	DIE_ON_ERR(scanner.Err(), "Couldn't read reads file to completion")

	// sort the strings and return
	sort.Strings(reads)
	log.Printf("Read %v reads; flipped %v of them.\n", len(reads), flipped)
	return reads
}

// listBuckets() processes the reads and creates the bucket list and the list
// of the bucket sizes and returns them.
func listBuckets(reads []string) ([]string, []int) {
	curBucket := ""
	buckets := make([]string, 0)
	counts := make([]int, 0)

	for _, r := range reads {
		if r[:globalK] != curBucket {
			curBucket = r[:globalK]
			buckets = append(buckets, curBucket)
			counts = append(counts, 1)
		} else {
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

// encodeSingleReadWithBucket() encodes a single read: uses a bucketing scheme
// for initial part, and arithmetic encoding for the rest.
func encodeSingleReadWithBucket(r string, hash KmerHash, coder *arithc.Encoder) {
	// encode rest using the reference probs
	context := r[:globalK]
	contextMer := stringToKmer(context)

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
func encodeWithBuckets(readFile, outBaseName string, hash KmerHash, coder *arithc.Encoder) int {
	// read the reads and flip as needed
	reads := readAndFlipReads(readFile, hash, flipReadsOption)

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
	waitForBuckets := make(chan bool)
	go func() {
		encodeKmersToFile(buckets, writer)
		close(waitForBuckets)
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
	waitForCounts := make(chan bool)
	go func() {
		writeCounts(countZ, len(reads[0]), counts)
		close(waitForCounts)
	}()

	/*** The main work to encode the read tails ***/
	log.Printf("Encoding reads, each of length %d ...", len(reads[0]))
	waitForReads := make(chan bool)
	go func() {
		for _, r := range reads {
			encodeSingleReadWithBucket(r, hash, coder)
		}
		close(waitForReads)
	}()

	// Wait for each of the coders to finish
	<-waitForBuckets
	<-waitForCounts
	<-waitForReads
	log.Printf("done.")
	return len(reads)
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

	var n int
	_, err := fmt.Fscanf(c, "%d", &readlen)
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

// dart() finds the interval in the given distribution that contains the given
// target, after transformming the distribution using the given weightOf
// function. This is called by lookup() during decode.
func dart(dist [len(ALPHA)]uint32, target uint32, weightOf WeightXformFcn) (uint64, uint64, uint64) {
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

// lookup() is called by arithc.Decoder to find an interval that contains the given value t.
func lookup(hash KmerHash, context Kmer, t uint64) (uint64, uint64, uint64) {
	if info, ok := hash[context]; ok {
		return dart(info.next, uint32(t), contextWeight)
	} else {
		return dart(defaultInterval, uint32(t), defaultWeight)
	}
}

// sumDist() computes the sum of the items in the given distribution after
// first transforming them via the given weightOf function.
func sumDist(d [len(ALPHA)]uint32, weightOf WeightXformFcn) (total uint64) {
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
		return sumDist(defaultInterval, defaultWeight)
	}
}

// decodeReads() decodes the file wrapped by the given Decoder, using the
// kmers, counts, and hash table provided. It writes its output to the given
// io.Writer.
func decodeReads(kmers []string, counts []int, hash KmerHash, readLen int, out io.Writer, decoder *arithc.Decoder) {
	log.Printf("Decoding reads...")

	buf := bufio.NewWriter(out)

	var contextMer Kmer
	lu := func(t uint64) (uint64, uint64, uint64) {
		return lookup(hash, contextMer, t)
	}

	n := 0
	curBucket := 0
	bucketCount := 0
	for curBucket < len(kmers) {
		// write the bucket
		buf.Write([]byte(kmers[curBucket]))

		contextMer = stringToKmer(kmers[curBucket])
		// write the reads
		for i := 0; i < readLen-len(kmers[0]); i++ {
			// decode next symbol
			symb, err := decoder.Decode(contextTotal(hash, contextMer), lu)
			DIE_ON_ERR(err, "Fatal error decoding!")
			b := byte(symb)

			// write it out
			buf.WriteByte(baseFromBits(b))

			// update hash counts (throws away the computed interval; just
			// called for side effects.)
			nextInterval(hash, contextMer, b)

			// update the new context
			contextMer = shiftKmer(contextMer, b)
		}

		// at the end of the read; write a newline and increment the # of reads
		// from this bucket that we've written
		buf.WriteByte('\n')
		bucketCount++
		n++

		if bucketCount == counts[curBucket] {
			curBucket++
			bucketCount = 0
		}
	}
	buf.Flush()
	log.Printf("done. Wrote %v reads.", n)
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
}

// writeGlobalOptions() writes out the global variables that can affect the
// encoding / decoding. Files encoded with one set of options can only be
// decoded using the same set of options.
func writeGlobalOptions() {
	log.Printf("Option: psudeoCount = %d", pseudoCount)
	log.Printf("Option: observationWeight = %d", observationWeight)
	log.Printf("Option: seenThreshold = %d", seenThreshold)
	log.Printf("Option: observationInc = %d", observationInc)
	log.Printf("Option: smoothOption = %v", smoothOption)
	log.Printf("Option: flipReadsOption = %v", flipReadsOption)
}

// main() encodes or decodes a set of reads based on the first command line
// argument (which is either encode or decode).
func main() {
	log.Println("Starting kpath version 5-28-14")

	log.Println("Maximum threads = 4")
	runtime.GOMAXPROCS(4)

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
	if globalK <= 0 {
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
	waitForReference := make(chan bool)
	go func() {
		hash = countKmersInReference(globalK, refFile)
		log.Printf("There are %v unique %v-mers in the reference\n",
			len(hash), globalK)
		capTransitionCounts(hash, 2)
		close(waitForReference)
	}()

	writeGlobalOptions()

	if mode == ENCODE {
		/* encode -k -ref -reads=FOO.seq -out=OUT
		   will encode into OUT.{enc,bittree,counts} */
		log.Printf("Reading from %s", readFile)
		log.Printf("Writing to %s, %s, %s", outFile+".enc", outFile+".bittree", outFile+".counts")

		var err error

		if smoothOption {
			smoothFile, err = os.Create("smoothed.txt")
			DIE_ON_ERR(err, "Couldn't create smoothed file 'smoothed.txt'")
			defer smoothFile.Close()
		}

		// create the output file
		outF, err := os.Create(outFile + ".enc")
		DIE_ON_ERR(err, "Couldn't create output file %s", outFile)
		defer outF.Close()

		writer := bitio.NewWriter(outF)
		defer writer.Close()

		// create encoder
		encoder := arithc.NewEncoder(writer)
		defer encoder.Finish()

		// encode reads
		<-waitForReference
		n := encodeWithBuckets(readFile, outFile, hash, encoder)
		log.Printf("Smoothed (changed) %v characters", smoothed)
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
		waitForBuckets := make(chan bool)
		go func() {
			kmers = decodeKmersFromFile(headsFN, globalK)
			sort.Strings(kmers)
			close(waitForBuckets)
		}()

		// read the bucket counts
		var counts []int
		var readlen int
		waitForCounts := make(chan bool)
		go func() {
			counts, readlen = readBucketCounts(countsFN)
			close(waitForCounts)
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
		log.Printf("Read length = %d", readlen)
		decodeReads(kmers, counts, hash, readlen, outF, decoder)
	}
	log.Printf("Default interval used %v times and context used %v times",
		defaultUsed, contextExists)
}

/*
// require seeing an item twice before counting it
x transcript sequences are assigned 2
* when we use a character in a context, increment by 1
x w = observationWeight * min64(uint64(dist[i]-1), 1)
x Smooth anytime count < 2
*/

/*
1. When encoding a read, check a few kmers randomly and
choose the strand that matches the most (done)

2. If you encounter a kmer you haven't seen in this context,
    record that you saw it
    output the maximum probability kmer
    If you see it a second time, change prob to match reference prob

3. When encoding the first k letters, we use a special prob distribution
(done)

4. update special prob distribution (done)

10  smoothed    13479801
12  smoothed    12106330
14  smoothed    6952281

head file:
14  smoothed  6702274
14  smoothed no pseudo   6688925
14  smoothed no pseudo, update-in-xcript, 4741447
14  no-smooth, pseudo=0, update-in-xscript, 5814616
    if context exists & interval is non-zero, we use it and increment the use
        we increment by 1, and devide the total counts by 2 so that we have to see a new
        transition twice before it gets a non-zero prob
    if context doesn't exist, or the interval is zero, we use default and increment it
14  no-smooth, pseudo=0, update by 1 until seen twice, then by 2 (transcriptomic transitiosn start at 2)
    size = 5783643

20  same as above, size = 6135591
18  same as above, size = 5667019
17  same as above, size = 5453441
16  same as above, size = 5295103
15  same as above, size = 5453441
14  same as above, size = 5783643

next idea: introduce new contexts when we need to:
    if context is not in hash table
    use default encoding
    then add context with 0-weight distribution and increment the item for the next guy

16  no-smoot, pseudo=0, update by 1 until seen twice, then by 2 (transcroptomic transistiosn start at 2)
        new contexts are added as above
    size = 5123801

16  all same as previous, but now with smoothing, size = 4511492

16  all same as previous, but with smoothing, size = 4511492 and bziped2 ascii smoothed file = 1972894
    total size = 6484386

16  same as previous, with smoothing, bzipped binary smoothed file; size = 6458275 = 4511492 + 1946783

16  same as prev, no smoothing, except increment default counters whenever used; size = 5124373

16  same as prev, no smooth, increment counters, use default weights direct with defaultWeight size = 5124373

16  same as prev, no smoot, added "start of read interval"; size = 5124413

next step:
    to really test benefit of smoothing, you have to record, in a separate file:
        delta since last smooth
        2 bit encoding of new character
        (or just delta since last A smooth, T smooth, C smooth, G smooth)
*/
