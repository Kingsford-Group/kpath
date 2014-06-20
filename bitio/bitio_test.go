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
