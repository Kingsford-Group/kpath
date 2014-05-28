package main

import (
	"testing"
)

func TestKmerConversion(t *testing.T) {
	var mers = []string{"AAAAAAAAAAACAAAC", "ACAGACGTAGACGTA", "ACAG", "TTATAT"}
	for _, m := range mers {
		x := stringToKmer(m)
		y := kmerToString(x, len(m))
		if y != m {
			t.Fatalf("%s != %s for kmer %v\n", y, m, x)
		}
	}
}
