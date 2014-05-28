package main

/* TODO:
x0. change pseudocount to be max (or else change the seen threshold)
x1. fix filename handling (including outputting filename information as a log)
x why isn't TTTT...TTT flipped? b/c it's tied

x2. compress counts and bittree directly from here (and uncompress bittree when reading it)
x3. move unused code over to unused.no file
x4. update error messages and panic messages to be more consistent
    // DIE_ON_ERROR(err, "Couldn't create bucket file: %v", err)

5. conserve memory with a DNAString type (?)
6. profile to speed up
7. parallelize
8. add some more unit tests
*/

import (
    "log"
    "os"
    "io"
    "fmt"
    "bufio"
    "strings"
    "flag"
    "sort"
    "compress/gzip"
    "path/filepath"

    "kingsford/arithc"
    "kingsford/bitio"
)

//===================================================================
// Kmer types

type Kmer   uint64

type KmerInfo struct {
    next    [4]uint32
}

type KmerHash   map[Kmer]*KmerInfo


//===================================================================
// Globals
//===================================================================

var (
    encodeFlags     *flag.FlagSet
    outFile         string
    refFile         string
    readFile        string
    globalK         int
    defaultInterval [4]uint32 = [...]uint32{2, 2, 2, 2}
    readStartInterval [4]uint32 = [...]uint32{2,2,2,2}

    defaultUsed     int
    contextExists   int
    smoothed        int
    flipped         int
)

// global variables for smoothing
var (
    lastSmooth      [4]uint64 = [...]uint64{1,2,3,4}
    charCount       uint64
    smoothFile      *os.File
    tmpByteSlice    []byte = make([]byte, 2)
)

const (
    pseudoCount         uint64 = 1
    observationWeight   uint64 = 10
    seenThreshold       uint32 = 2 // before this threshold, increment 1 and treat as unseen
    observationInc      uint32 = 2 // once above seenThreshold, increment by this on each observation

    smoothOption        bool = false
    flipReadsOption     bool = true
)

//===================================================================
// Int <-> String Kmer representations
//===================================================================

func acgt(a byte) byte {
    switch a {
        case 'A': return 0
        case 'N': return 0
        case 'C': return 1
        case 'G': return 2
        case 'T': return 3
    }
    panic("Bad character!")

    return 0  // invalid chars -> 'A'
}

func baseFromBits(a byte) byte {
    return "ACGT"[a]
}

func stringToKmer(kmer string) Kmer {
    var x uint64
    for _, c := range kmer {
        x = (x<<2) | uint64(acgt(byte(c)))
    }
    return Kmer(x)
}

func isACGT(c rune) bool {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T'
}

func kmerToString(kmer Kmer, k int) string {
    s := make([]byte, k)
    for i := 0; i < k; i++ {
        s[k - i - 1] = baseFromBits(byte(kmer & 0x3))
        kmer >>= 2
    }
    return string(s)
}

/* RC computes the reverse complement of the single nucleotide */
func RC(c byte) byte {
    switch c {
        case 'A': return 'T'
        case 'C': return 'G'
        case 'G': return 'C'
        case 'T': return 'A'
    }
    return 'T'
}

/* reverseComplement returns the reverse complement of the string */
func reverseComplement(r string) string {
    s := make([]byte, len(r))
    for i := 0; i < len(r); i++ {
        s[len(r) - i - 1] = RC(r[i])
    }
    return string(s)
}

//===================================================================

func newKmerHash() KmerHash {
    return make(KmerHash)
}

/* Read a multifasta file into a slice of strings */
func readFastaFile(fastaFile string) []string {
    // open the .gz fasta file that is the references
    log.Println("Reading Reference File...")
    inFasta, err := os.Open(fastaFile)
    DIE_ON_ERR(err, "Couldn't open fasta file %s", fastaFile)
    defer inFasta.Close()

    // wrap the gzip reader around it
    in, err := gzip.NewReader(inFasta)
    DIE_ON_ERR(err, "Couldn't open gzipped file %s", fastaFile)
    defer in.Close()

    out := make([]string,0)
    cur := make([]string,0)

    scanner := bufio.NewScanner(in)
    for scanner.Scan() {
        line := strings.TrimSpace(strings.ToUpper(scanner.Text()))
        if len(line) > 0 && line[0] == '>' {
            if len(cur) > 0 {
                out = append(out, strings.Join(cur, ""))
                cur = make([]string, 0)
            } 
        } else if len(line) > 0 {
            cur = append(cur, line)
        }
    }
    if err := scanner.Err(); err != nil {
        log.Fatal(err)
    }
    return out
}

/* Count the kmers and the distribution of the next character in the
given fasta file */
func countKmers(k int, fastaFile string) KmerHash {
    seqs := readFastaFile(fastaFile)
    hash := newKmerHash()

    log.Printf("Counting %v-mer transitions in reference file...\n", k)
    for _, s := range seqs {
        for i := 0; i < len(s) - k; i++ {
            km := stringToKmer(s[i : i + k])
            _, ok := hash[km]
            if !ok {
                hash[km] = &KmerInfo{}
            }
            hash[km].next[acgt(s[i + k])] += observationInc
        }
    }
    return hash
}

/* make all transition counts at max */
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

func min64(a uint64, b uint64) uint64 {
    if a < b { 
        return a 
    } 
    return b
}

func maxChar(dist [4]uint32) byte {
    curMax := uint32(0)
    curR := 0
    for i, m := range dist {
        if m > curMax {
            curMax = m
            curR = i
        }
    }
    return baseFromBits(byte(curR))
}

func writeVarLenInt(out *os.File, b uint64) {
    if b < 128 {
        tmpByteSlice[0] = byte(b)
        tmpByteSlice[1] = 0
        out.Write(tmpByteSlice[0:1])
    } else if b < uint64(2)<<31 {
        tmpByteSlice[0] =  byte(0x80 | (b >> 8))
        tmpByteSlice[1] = byte(0x00FF & b)
        out.Write(tmpByteSlice)
    }
}

func smoothError(hash KmerHash, context string, next byte) byte {
    charCount++
    if info, ok := hash[stringToKmer(context)]; ok {
        letterIdx := int(acgt(next))
        if info.next[letterIdx] < seenThreshold {
            smoothed++
            //fmt.Fprintf(smoothFile, "%d ", charCount - lastSmooth[letterIdx]-1)
            writeVarLenInt(smoothFile, charCount - lastSmooth[letterIdx]-1)
            lastSmooth[letterIdx] = charCount
            return maxChar(info.next)
        }
    }
    return next
}

/*
// require seeing an item twice before counting it
x transcript sequences are assigned 2
* when we use a character in a context, increment by 1
x w = observationWeight * min64(uint64(dist[i]-1), 1) 
x Smooth anytime count < 2
*/

func contextWeight(charIdx int, dist[4]uint32) (w uint64) {
    if dist[charIdx] >= seenThreshold {
        w = observationWeight * uint64(dist[charIdx]) / uint64(observationInc)
        //if w < pseudoCount { w = pseudoCount }
        return
    }
    w = pseudoCount
    return
}

func defaultWeight(charIdx int, dist [4]uint32) uint64 {
    return uint64(dist[charIdx])
}

type WeightXformFcn func(int,[4]uint32)uint64

/* return the partial sum for the given character */
func intervalFor(nextChar byte, dist [4]uint32, weightOf WeightXformFcn) (a uint64, b uint64, total uint64) { 
    letterIdx := int(acgt(nextChar))
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

/* Compute the interval for the given context */
func nextInterval(hash KmerHash, context string, next byte) (a uint64, b uint64, total uint64) {
    kidx := acgt(next)
    contextMer := stringToKmer(context)
    info, ok := hash[contextMer]
    // if the context exists, use that distribution
    if ok {
        contextExists++
        a, b, total = intervalFor(next, info.next, contextWeight)
        if info.next[kidx] >= seenThreshold {  // increment double if in the transcriptome
            info.next[kidx] += observationInc
        } else {
            info.next[kidx]++
        }
    } else {
        // if the context doesnt exist, use a simple default interval
        defaultUsed++
        if ok { contextExists-- }
        a, b, total = intervalFor(next, defaultInterval, defaultWeight)
        defaultInterval[kidx]++

        // add this to the context now
        hash[contextMer] = &KmerInfo{}
        hash[contextMer].next[kidx]++
    }
    return
}

/* count the number of observaions of kmers in the read */
func countMatchingObservations(hash KmerHash, r string) (n uint32) {
    context := r[ : globalK]
    for i := globalK; i <= len(r) - globalK; i++ {
        if H, ok := hash[stringToKmer(context)]; ok {
            n += H.next[acgt(r[i])]
        }
        context = context[1:] + string(r[i])
    }

    return
}

/* count the number of kmers present in the context */
func countMatchingContexts(hash KmerHash, r string) (n int) {
    for i := 0; i <= len(r) - globalK; i++ {
        if _, ok := hash[stringToKmer(r[i : i + globalK])]; ok {
            n++
        }
    }
    return
}

func readAndFlipReads(readFile string, hash KmerHash, flipReadsOption bool) []string {
    // open the read file
    log.Println("Reading and flipping reads...")
    in, err := os.Open(readFile)
    DIE_ON_ERR(err, "Couldn't open read file %s", readFile)
    defer in.Close()

    // put the reads into a global array, flipped if needed
    reads := make([]string, 0)
    scanner := bufio.NewScanner(in)
    for scanner.Scan() {
        // remove spaces and convert on-ACGT to 'A'
        r := strings.Replace(strings.TrimSpace(strings.ToUpper(scanner.Text())), "N", "A", -1)
        if flipReadsOption {
            n1 := countMatchingContexts(hash, r)
            rcr := reverseComplement(r)
            n2 := countMatchingContexts(hash, rcr)
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

/* encode a single read: uses 1 scheme for initial part, and 1 scheme for the rest */
func encodeSingleReadWithBucket(r string, hash KmerHash, coder *arithc.Encoder) {
    // encode rest using the reference probs
    context := r[ : globalK]
    for i := globalK; i < len(r); i++ {
        char := r[i]
        a, b, total := nextInterval(hash, context, char)
        err := coder.Encode(a, b, total)
        DIE_ON_ERR(err, "Error encoding read: %s", r)
        context = context[1:] + string(char)
    }
}

// return the buckets and their counts
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

func writeCounts(f io.Writer, counts []int) {
    log.Printf("Writing counts...")
    for _, c := range counts {
        fmt.Fprintf(f, "%d ", c)
    }
    log.Printf("Done; write %d counts.", len(counts))
}

func encodeWithBuckets(readFile, outBaseName string, hash KmerHash, coder *arithc.Encoder) int {

    // read the reads and flip as needed
    reads := readAndFlipReads(readFile, hash, flipReadsOption)

    // create the buckets and counts
    buckets, counts := listBuckets(reads)

    // write the bittree for the bucket out to a file
    outBT, err := os.Create(outBaseName + ".bittree")
    DIE_ON_ERR(err, "Couldn't create bucket file: %s", outBaseName + ".bittree")
    defer outBT.Close()

    // compress the file with gzip as we are writing it
    outBZ, err := gzip.NewWriterLevel(outBT, gzip.BestCompression)
    DIE_ON_ERR(err, "Couldn't create gzipper for bucket file")
    defer outBZ.Close()

    // create a writer that lets us write bits
    writer := bitio.NewWriter(outBZ)
    defer writer.Close()

    /*** The main work to encode the bucket names ***/
    encodeKmersToFile(buckets, writer)

    // write out the counts
    countF, err := os.Create(outBaseName + ".counts")
    DIE_ON_ERR(err, "Couldn't create counts file: %s", outBaseName + ".counts")
    defer countF.Close()

    // compress it as we are writing it
    countZ, err := gzip.NewWriterLevel(countF, gzip.BestCompression)
    DIE_ON_ERR(err, "Couldn't create gzipper for count file")
    defer countZ.Close()

    /*** The main work to encode the bucket counts ***/
    writeCounts(countZ, counts)
    
    /*** The main work to encode the read tails ***/
    log.Printf("Encoding reads...")
    for _, r := range reads {
        encodeSingleReadWithBucket(r, hash, coder)
    }
    log.Printf("done.")
    return len(reads)
}


//===============================================================================
// DECODING
//===============================================================================

func readBucketCounts(countsFN string) []int {
    log.Printf("Reading bucket counts from %v", countsFN)

    // open the count file
    c1, err := os.Open(countsFN)
    DIE_ON_ERR(err, "Couldn't open count file: %s", countsFN)
    defer c1.Close()

    // the count file is compressed with gzip; uncompress it as we read it
    c, err := gzip.NewReader(c1)
    DIE_ON_ERR(err, "Couldn't create gzip reader: %v")
    defer c.Close()

    counts := make([]int, 0)
    var n int
    err = nil
    for x := 1; err == nil && x > 0; {
        x, err = fmt.Fscanf(c, "%d", &n)
        if x > 0 && err == nil {
            counts = append(counts, n)
        }
    }
    log.Printf("done; read %d counts", len(counts))
    return counts
}

func dart(dist [4]uint32, target uint32, weightOf WeightXformFcn) (uint64, uint64, uint64) {
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

// look up an interval that contains this total
func lookup(hash KmerHash, context string, t uint64) (c uint64, d uint64, s uint64) {
    ctx := stringToKmer(context)
    info, ok := hash[ctx]
    if ok {
        c, d, s = dart(info.next, uint32(t), contextWeight)
    }
    if !ok || c == d {
        return dart(defaultInterval, uint32(t), defaultWeight)
    }
    return
}

func sumDist(d [4]uint32, weightOf WeightXformFcn) (total uint64) {
    for i := range d {
        total += uint64(weightOf(i, d))
    }
    return
}

func contextTotal(hash KmerHash, context string) uint64 {
    info, ok := hash[stringToKmer(context)]
    if ok {
        return sumDist(info.next, contextWeight)
    } else {
        return sumDist(defaultInterval, defaultWeight)
    }
}

func decodeReads(kmers []string, counts []int, hash KmerHash, readLen int, out io.Writer, decoder *arithc.Decoder) {
    log.Printf("Decoding reads...")

    buf := bufio.NewWriter(out)

    var context string
    lu := func(t uint64) (uint64, uint64, uint64) {
        return lookup(hash, context, t)
    }

    n := 0
    curBucket := 0
    bucketCount := 0
    for curBucket < len(kmers) {
        // write the bucket
        buf.Write([]byte(kmers[curBucket]))

        // write the reads
        context = kmers[curBucket]
        for i := 0; i < readLen - len(kmers[0]); i++ {
            // decode next symbol
            symb, err := decoder.Decode(contextTotal(hash, context), lu)
            DIE_ON_ERR(err, "Fatal error decoding!")

            // write it out
            next := baseFromBits(byte(symb))
            buf.WriteByte(next)

            // update hash counts (throws away the computed interval; just
            // called for side effects.)
            nextInterval(hash, context, next) 
            
            // update the new context
            context = context[1:] + string(next)
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

func DIE_ON_ERR(err error, msg string, args... interface{}) {
    if err != nil {
        log.Printf("Error: " + msg, args...)
        log.Fatalln(err)
    }
}

func init() {
    encodeFlags = flag.NewFlagSet("encode", flag.ContinueOnError)
    encodeFlags.StringVar(&refFile, "ref", "", "reference fasta filename")
    encodeFlags.StringVar(&outFile, "out", "", "output filename")
    encodeFlags.StringVar(&readFile, "reads", "", "reads filename")
    encodeFlags.IntVar(&globalK, "k", 16, "length of k")
}

const (
    ENCODE int = 1
    DECODE int = 2
)

func removeExtension(filename string) string {
    extension := filepath.Ext(filename)
    return strings.TrimRight(filename, extension)
}

func main() {
    log.SetPrefix("kpath: ")
    log.Println("Starting kpath version 5-27-14")

    // parse the command line
    if len(os.Args) < 2 {
        encodeFlags.PrintDefaults()
        os.Exit(1)
    }
    var mode int
    if os.Args[1][0] == 'e' {
        mode = ENCODE
    } else {
        mode = DECODE
    }
    encodeFlags.Parse(os.Args[2:])
    if globalK <= 0 {
        log.Fatalf("K must be specified as a small positive integer with -k")
    }

    // count the kmers in the reference
    hash := countKmers(globalK, refFile)
    log.Printf("There are %v unique %v-mers in the reference\n", 
        len(hash), globalK)
    capTransitionCounts(hash, 2)

    if mode == ENCODE {
        /* encode -k -ref -reads=FOO.seq -out=OUT
            will encode into OUT.{enc,bittree,counts} */
        log.Printf("Reading from %s", readFile)
        log.Printf("Writing to %s, %s, %s", outFile + ".enc", outFile + ".bittree", outFile + ".counts")

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
        kmers := decodeKmersFromFile(headsFN, globalK)
        sort.Strings(kmers)

        // read the bucket counts
        counts := readBucketCounts(countsFN)

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

        READLEN := 35
        decodeReads(kmers, counts, hash, READLEN, outF, decoder)
    }
    log.Printf("Default interval used %v times and context used %v times", 
        defaultUsed, contextExists)
}


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
