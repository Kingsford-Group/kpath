package main

import (
    "log"
    "os"
    "bufio"
    "strings"

    "kingsford/arithc"
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

/* Return true if transition has non-zero prob */
func transitionHasNonZeroProb(hash KmerHash, context string, next byte) bool {
    h, ok := hash[stringToKmer(context)]
    return ok && h.next[acgt(next)] > 0
}
