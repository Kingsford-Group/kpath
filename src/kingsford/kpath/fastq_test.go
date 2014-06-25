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
	"fmt"
	"testing"
)

func TestReadFastQ(t *testing.T) {
	records := make(chan *FastQ)

	fn := "test.fq"
	fmt.Println(fn)
	go ReadFastQ(fn, records)

	for fq := range records {
		PrintFastQ(fq)
	}
}
