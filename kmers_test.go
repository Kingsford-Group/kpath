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
