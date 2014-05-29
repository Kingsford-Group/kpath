package main

import (
	"bufio"
	"log"
	"os"
	"strings"
	"path/filepath"

	"kingsford/arithc"
)

// global variables for smoothing
var (
	lastSmooth   [len(ALPHA)]uint64 = [...]uint64{1, 2, 3, 4}
	charCount    uint64
	smoothFile   *os.File
	tmpByteSlice []byte = make([]byte, 2)
)

/* encode a single read: uses 1 scheme for initial part, and 1 scheme for the rest */
func encodeSingleRead(r string, hash KmerHash, coder *arithc.Encoder) error {
	// guess whether we should reverse complement the read
	n1 := countMatchingObservations(hash, r)
	rcr := reverseComplement(r)
	n2 := countMatchingObservations(hash, rcr)
	if n2 > n1 {
		r = rcr
		flipped++
	}

	var i int
	// for early characters in the read, use the readStart interval
	for i = 0; i < globalK; i++ {
		a, b, total := intervalFor(r[i], readStartInterval, defaultWeight)
		readStartInterval[acgt(r[i])]++
		err := coder.Encode(a, b, total)
		if err != nil {
			return err
		}
	}

	// encode rest using the reference probs
	context := r[:globalK]
    contextMer := stringToKmer(context)
	for ; i < len(r); i++ {
		char := acgt(r[i])
		if smoothOption {
			char = smoothError(hash, contextMer, char)
		}
		a, b, total := nextInterval(hash, contextMer, char)
		err := coder.Encode(a, b, total)
		if err != nil {
			return err
		}
		//context = context[1:] + string(char)
        contextMer = shiftKmer(contextMer, char)
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

/* Return true if transition has non-zero prob */
func transitionHasNonZeroProb(hash KmerHash, context string, next byte) bool {
	h, ok := hash[stringToKmer(context)]
	return ok && h.next[acgt(next)] > 0
}

func removeExtension(filename string) string {
	extension := filepath.Ext(filename)
	return strings.TrimRight(filename, extension)
}

func min64(a uint64, b uint64) uint64 {
	if a < b {
		return a
	}
	return b
}

func maxChar(dist [len(ALPHA)]uint32) byte {
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
		tmpByteSlice[0] = byte(0x80 | (b >> 8))
		tmpByteSlice[1] = byte(0x00FF & b)
		out.Write(tmpByteSlice)
	}
}

func smoothError(hash KmerHash, contextMer Kmer, next byte) byte {
	charCount++
	if info, ok := hash[contextMer]; ok {
		letterIdx := int(acgt(next))
		if info.next[letterIdx] < seenThreshold {
			smoothed++
			//fmt.Fprintf(smoothFile, "%d ", charCount - lastSmooth[letterIdx]-1)
			writeVarLenInt(smoothFile, charCount-lastSmooth[letterIdx]-1)
			lastSmooth[letterIdx] = charCount
			return maxChar(info.next)
		}
	}
	return next
}

// countMatchingContexts() counts the number of kmers present in the hash.
func countMatchingContexts(hash KmerHash, r string) (n int) {
    /*
	context := r[:globalK]
	for i := globalK; i <= len(r)-globalK; i++ {
		if _, ok := hash[stringToKmer(context)]; ok {
			n++
		}
		context = context[1:] + string(r[i])
	}
	return */

    contextMer := stringToKmer(r[:globalK])
	for i := 0; i <= len(r)-globalK; i++ {
        if _, ok := hash[contextMer]; ok {
            n++
        }
        if i+globalK < len(r) {
            contextMer = shiftKmer(contextMer, r[i+globalK])
        }
        /*
		if _, ok := hash[stringToKmer(r[i:i+globalK])]; ok {
			n++
		} */
	} 
	return
}


