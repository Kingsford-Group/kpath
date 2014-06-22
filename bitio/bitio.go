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
	"io"
)

type Writer struct {
	buf     *bufio.Writer // where we write; should this just be io.Writer?
	curByte byte          // our local (small) buffer
	nbits   byte          // number of bits written to in curByte so far
}

// NewWriter() creates a new bit writer backed by the given io.Writer.
func NewWriter(w io.Writer) *Writer {
	return &Writer{bufio.NewWriter(w), 0, 0}
}

// WriteBit() writes out a single bit to the stream.
func (bw *Writer) WriteBit(b byte) error {
	bw.curByte |= b << (7 - bw.nbits)
	bw.nbits++

	// if we've filled up this byte
	var err error
	if bw.nbits == 8 {
		err = bw.buf.WriteByte(bw.curByte)
		if err != nil {
			panic(err)
		}
		bw.curByte = 0
		bw.nbits = 0
	}
	return err
}

/*
// Implement the Writer interface
// Note: this treats every BYTE as a bit
func (bw *Writer) Write(bb []byte) error {
    for _,b := range bb {
        err := bw.WriteBit(b)
        if err != nil {
            return err
        }
    }
    return nil
}
*/

// Close() close the bit writer and flushes the last byte.
func (bw *Writer) Close() error {
	if bw.nbits > 0 {
		err := bw.buf.WriteByte(bw.curByte)
		if err != nil {
			return err
		}
	}
	return bw.buf.Flush()
}

// Represents a reader that can read at a bit level
type Reader struct {
	reader  *bufio.Reader
	curByte byte
	curBit  int
}

// NewReader() creates a new bit reader that wraps the underlying bufio reader.
func NewReader(r *bufio.Reader) *Reader {
	return &Reader{r, 0, -1}
}

// ReadBit() reads a single bit from the file.  Bits are represented as bytes
// that have value 0 or 1 returns non-nil error on EOF.
func (br *Reader) ReadBit() (byte, error) {
	if br.curBit < 0 {
		by, err := br.reader.ReadByte()
		br.curByte = by
		if err != nil {
			return 0, err
		}
		br.curBit = 7
	}

	b := (br.curByte >> uint(br.curBit)) & 1
	br.curBit--
	return b, nil
}

// Close() closes the file --- a nop for now since we don't need to close anything.
func (br *Reader) Close() error {
	return nil
}
