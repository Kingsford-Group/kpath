package main

import (
    "log"
    "os"
    "bufio"
    "strings"
    "flag"

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
    pseudoCount         uint64 = 0
    observationWeight   uint64 = 10
    seenThreshold       uint32 = 2 // before this threshold, increment 1 and treat as unseen
    observationInc      uint32 = 2 // once above seenThreshold, increment by this on each observation

    smoothOption        bool = false
)

//===================================================================
// Int <-> String Kmer representations
//===================================================================

func acgt(a byte) byte {
    switch a {
        case 'A': return 0
        case 'C': return 1
        case 'G': return 2
        case 'T': return 3
    }
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

func kmerToString(kmer Kmer, k int) string {
    s := make([]byte, k)
    for i := 0; i < k; i++ {
        s[k - i - 1] = byte(kmer & 0x3)
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
    log.Println("Reading Reference File...")
    in, err := os.Open(fastaFile)
    if err != nil {
        log.Fatalf("Couldn't open fasta file %s: %v", fastaFile, err)
    }
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


/* Return true if transition has non-zero prob */
func transitionHasNonZeroProb(hash KmerHash, context string, next byte) bool {
    h, ok := hash[stringToKmer(context)]
    return ok && h.next[acgt(next)] > 0
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
        w = observationWeight * uint64(dist[charIdx]) / uint64(observationInc) + pseudoCount
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
    } 
    
    // if context doesn't exist, or gave a 0 prob, use default
    if !ok || a == b {
        // if the context doesnt exist, use a simple default interval
        defaultUsed++
        if ok { contextExists-- }
        a, b, total = intervalFor(next, defaultInterval, defaultWeight)
        defaultInterval[kidx]++
    }

    // if context didn't exist, create it
    if !ok {
        // add this to the context now
        hash[contextMer] = &KmerInfo{}
        hash[contextMer].next[kidx]++
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


/* encode a single read: uses 1 scheme for initial part, and 1 scheme for the rest */
func encodeSingleRead(r string, hash KmerHash, coder *arithc.Encoder) error {
    // guess whether we should reverse complement the read
    n1 := countMatchingContexts(hash, r)
    rcr := reverseComplement(r)
    n2 := countMatchingContexts(hash, rcr)
    if n2 > n1 {
        r = rcr
        flipped++
    }

    var i int
    // for early characters in the read, use the default interval
    for i = 0; i < globalK; i++ {
        a, b, total := intervalFor(r[i], readStartInterval, defaultWeight)
        readStartInterval[acgt(r[i])]++
        err := coder.Encode(a, b, total)
        if err != nil { return err }
    }

    // encode rest using the reference probs
    context := r[ : globalK]
    for ; i < len(r); i++ {
        char := r[i]
        if smoothOption { 
            char = smoothError(hash, context, char)
        }
        a, b, total := nextInterval(hash, context, char)
        err := coder.Encode(a, b, total)
        if err != nil { return err }
        context = context[1:] + string(char)
    }
    return nil
}


/* Read each line as a read and encode it using encodeSingleRead */
func encodeReads(readFile string, hash KmerHash, coder *arithc.Encoder) (n int) {
    log.Println("Encoding reads...")
    in, err := os.Open(readFile)
    if err != nil {
        log.Fatalf("Couldn't open reads file %s: %v", readFile, err)
    }
    defer in.Close()

    scanner := bufio.NewScanner(in)
    for scanner.Scan() {
        line := strings.TrimSpace(strings.ToUpper(scanner.Text()))
        err := encodeSingleRead(line, hash, coder)
        if err != nil {
            log.Fatalf("Error encoding read %v\n%v", line, err)
        }
        n++
    }
    if err := scanner.Err(); err != nil {
        log.Fatalf("Couldn't read reads file: %v", err)
    }
    return n
}


//===================================================================
// Command line and main driver
//===================================================================

func init() {
    encodeFlags = flag.NewFlagSet("encode", flag.ContinueOnError)
    encodeFlags.StringVar(&outFile, "out", "", "output filename")
    encodeFlags.StringVar(&refFile, "ref", "", "reference fasta filename")
    encodeFlags.StringVar(&readFile, "reads", "", "reads filename")
    encodeFlags.IntVar(&globalK, "k", 0, "length of k")
}

func main() {
    log.Println("Starting kpath version 4-20-14")

    if len(os.Args) < 2 {
        encodeFlags.PrintDefaults()
        os.Exit(1)
    }
    encodeFlags.Parse(os.Args[1:])

    // count the kmers in the reference
    hash := countKmers(globalK, refFile)
    log.Printf("There are %v unique %v-mers in the refernece\n", 
        len(hash), globalK)
    capTransitionCounts(hash, 2)

    // create the output file
    outF, err := os.Create(outFile)
    if err != nil {
        log.Fatalf("Couldn't create output file %s: %v", outFile, err)
    }
    defer outF.Close()
    writer := bitio.NewWriter(outF)
    defer writer.Close()

    if smoothOption {
        smoothFile, err = os.Create("smoothed.txt")
        if err != nil {
            log.Fatalf("Couldn't create smoothed file: %v", err)
        }
        defer smoothFile.Close()
    }

    // create encoder
    encoder := arithc.NewEncoder(writer)
    defer encoder.Finish()

    // encode reads
    n := encodeReads(readFile, hash, encoder)
    log.Printf("Encoded %v reads; default interval used %v times (context used %v times)\n", 
        n, defaultUsed, contextExists)
    log.Printf("Smoothed (changed) %v characters\n", smoothed)
    log.Printf("Reads Flipped: %v\n", flipped)
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
