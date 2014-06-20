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
