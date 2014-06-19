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

func TestKmerShift(t *testing.T) {
	globalK = 5
	setShiftKmerMask()
	m1 := kmerToString(shiftKmer(stringToKmer("TTCGT"), acgt(byte('G'))), globalK)
	if m1 != "TCGTG" {
		t.Fatalf("%s != %s (1)", m1, "TCGTG")
	}
	m2 := kmerToString(shiftKmer(stringToKmer("TTCGT"), acgt(byte('C'))), globalK)
	if m2 != "TCGTC" {
		t.Fatalf("%s != %s (2)", m1, "TCGTC")
	}
}
