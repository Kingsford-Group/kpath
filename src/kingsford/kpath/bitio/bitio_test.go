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

package bitio

import (
	"bufio"
	"os"
	"testing"
)

func TestWriteFile(t *testing.T) {
	// output an output file
	f, err := os.Create("testout.bit")
	if err != nil {
		t.Fatalf("Couldn't open file: %v", err)
	}
	defer f.Close()

	// wrap it in a bit writer
	bw := NewWriter(f)
	defer bw.Close()

	// write out some bits
	data := []byte{0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1}
	for _, b := range data {
		err := bw.WriteBit(b)
		if err != nil {
			t.Fatalf("Couldn't write bit: %v", err)
		}
	}
}

func TestWriteRead(t *testing.T) {
	f, err := os.Create("testout.bit")
	if err != nil {
		t.Fatalf("Couldn't create file: %v", err)
	}

	// wrap it in a bit writer
	bw := NewWriter(f)

	// write out some bits
	data := []byte{0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1}
	for _, b := range data {
		err := bw.WriteBit(b)
		if err != nil {
			t.Fatalf("Couldn't write bit: %v", err)
		}
	}

	bw.Close()
	f.Close()

	f, err = os.Open("testout.bit")
	if err != nil {
		t.Fatalf("Could open the file we wrote %v", err)
	}
	defer f.Close()

	bi := bufio.NewReader(f)
	br := NewReader(bi)
	defer br.Close()

	for i := 0; i < len(data) && err == nil; i++ {
		b, err := br.ReadBit()
		if err != nil {
			t.Fatalf("Couldn't read! %v", err)
		}
		if data[i] != b {
			t.Fatalf("Bits differ!")
		}
	}
}
