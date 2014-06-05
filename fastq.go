package main

import (
    "os"
    "bufio"
    "strings"
    "fmt"
)

// Represents a fastQ record
type FastQ struct {
    Seq []byte
    Quals []byte
    NLocations []byte
    IsFlipped bool
}

// NewFastQ creates a new, empty fastq record
func NewFastQ(seq []byte, quals []byte) *FastQ {
    f := FastQ{ 
       Seq: make([]byte, len(seq)),
       Quals: make([]byte, len(quals)),
       NLocations: make([]byte, 0),
       IsFlipped: false,
    }
    copy(f.Seq, seq)
    copy(f.Quals, quals)
    f.RemoveNs()
    return &f
}

// RemoveNs replaces any 'N's in the sequence with 'A' and records the position
// of the Ns in NLocations.
func (q *FastQ) RemoveNs() {
    for i, c := range q.Seq {
        if c == 'N' {
            q.Seq[i] = 'A'
            q.NLocations = append(q.NLocations, byte(i))
        }
    }
}

// ReverseComplement() will reverse complement a FastQ record, including
// reversing its quality values and updating it's Nlocations.
func (q *FastQ) ReverseComplement() {
    // reverse complement the sequence
    q.Seq = []byte(reverseComplement(string(q.Seq)))

    // reverse the quality array
    for i, j := 0, len(q.Quals)-1; i < j; i, j = i+1, j-1 {
        q.Quals[i], q.Quals[j] = q.Quals[j], q.Quals[i]
    }
    
    // reverse complement the locations
    for i, v := range q.NLocations {
        q.NLocations[i] = byte(len(q.Seq)) - v - 1
    }

    // record that we flipped
    q.IsFlipped = true
}


// PrintFastQ prints out the fastq record (used only for debugging).
func PrintFastQ(q *FastQ) {
    fmt.Println(string(q.Seq))
    fmt.Println(string(q.Quals))
    fmt.Printf("%v\n", q.NLocations)
}


// ReadFastQ reads fastq records from the file and pushes them out along the
// given channel. It will remove Ns from the sequence and replace them with As.
func ReadFastQ(filename string, out chan<- *FastQ) {
    // open the file
    fmt.Println(filename)
	in, err := os.Open(filename)
	DIE_ON_ERR(err, "Couldn't open read file %s", filename)
	defer in.Close()

    const (BETWEEN int = iota; INSEQ; INQUALS)
    state := BETWEEN

    seq := make([]byte, 0)
    quals := make([]byte, 0)
    var emptyQuals = make([]byte, 0)

    scanner := bufio.NewScanner(in)
    for scanner.Scan() {
        // read a line, remove white space
        r := strings.TrimSpace(strings.ToUpper(scanner.Text()))
        if len(r) == 0 { continue }

        // depending on state, manage record
        switch {

        case state == BETWEEN && r[0] == '@':
            seq = seq[0:0]
            quals = quals[0:0]
            state = INSEQ

        case state == INSEQ && r[0] == '+':
            state = INQUALS

        case state == INSEQ:
            seq = append(seq, []byte(r)...)

        case state == INQUALS:
            quals = append(quals, []byte(r)...)

            if len(quals) >= len(seq) {
                state = BETWEEN
                if writeQualOption {
                    out <- NewFastQ(seq, quals)
                } else {
                    out <- NewFastQ(seq, emptyQuals)
                }
            }
        }
    }
	DIE_ON_ERR(scanner.Err(), "Couldn't read reads file to completion")
    close(out)
}

