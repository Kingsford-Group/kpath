package arithc

import (
	"kingsford/kpath/bitio"
)

// A function that returns the current interval for the given symbol
// (c,d,total) meaning that the probabilty interval is
// [c/total, (c+d)/total)
type IntervalFunc func(int16) (uint64, uint64, uint64)

// Arithmetic coder for encoding
// Create encoder via NewEncoder
type Encoder struct {
	writer           *bitio.Writer
	width            uint64
	lo               uint64
	bits_outstanding uint64
}

const (
	moffetB         uint8  = 64
	halfInterval    uint64 = 1 << (moffetB - 1)
	quarterInterval uint64 = 1 << (moffetB - 2)
)

// NewEncoder() sreates a new arithmetic coder that will output to the given bit writer
func NewEncoder(bw *bitio.Writer) *Encoder {
	return &Encoder{bw, halfInterval, 0, 0}
}

// outputBitPlusFollow() outputs the bits we know for sure at this point
func (ac *Encoder) outputBitPlusFollow(b byte) error {
	ac.writer.WriteBit(b)
	if b != 0 {
		b = 0
	} else {
		b = 1
	}
	for ac.bits_outstanding > 0 {
		ac.writer.WriteBit(b)
		ac.bits_outstanding--
	}
	return nil
}

// renormalize() outputs the known bits and readjust the range
func (ac *Encoder) renormalize() error {
	for ac.width <= quarterInterval {
		// These cases must be in this order to handle the overflow situation
		if ac.lo >= halfInterval {
			ac.outputBitPlusFollow(1)
			ac.lo -= halfInterval
		} else if ac.lo+ac.width <= halfInterval {
			ac.outputBitPlusFollow(0)
		} else {
			ac.bits_outstanding++
			ac.lo -= quarterInterval
		}
		ac.lo <<= 1
		ac.width <<= 1
	}
	return nil
}

// Encode() is the main entry point for encoding; it encodes the symbol to the
// stream. (c,d) should be the symbols range in a discret distribution of
// support [0,total].
func (ac *Encoder) Encode(c, d, total uint64) error {
	//if d <= c { panic("0-length range in Encode!") }

	// update the range (lo, width)
	r := ac.width / total
	ac.lo += r * c
	if d < total {
		ac.width = r * (d - c)
	} else {
		ac.width -= c * r
	}
	return ac.renormalize()
}

// Finish() flushes any bits to the stream. Encoding is done after this is
// called. It output a number in the final interval
func (ac *Encoder) Finish() error {
	for i := uint8(1); i <= moffetB; i++ {
		err := ac.outputBitPlusFollow(byte((ac.lo >> (moffetB - i)) & 1))
		if err != nil {
			return err
		}
	}
	return nil
}

//============================================================

// Find a symbol s with interval (c,d) such that c < freq(s) < d
type LookupFunc func(uint64) (uint64, uint64, uint64)

type Decoder struct {
	reader *bitio.Reader
	inbuf  uint64
	width  uint64
}

// NewDecoder() creates a new decoder to read from a bit stream.
func NewDecoder(r *bitio.Reader) (*Decoder, error) {
	var d uint64
	for i := uint8(1); i <= moffetB; i++ {
		b, err := r.ReadBit()
		if err != nil {
			return nil, err
		}
		d = (d << 1) + uint64(b)
	}
	return &Decoder{r, d, halfInterval}, nil
}

// min64() computes the minimum of 2 uint64s.
func min64(a uint64, b uint64) uint64 {
	if a <= b {
		return a
	} else {
		return b
	}
}

// Decode() decodes a single symbol from the stream. Total is the support
// of the current distribution, and lookup is a a function that will return
// a range in the distribution that a given value falls into.
func (ad *Decoder) Decode(total uint64, lookup LookupFunc) (uint64, error) {
	r := ad.width / total
	v := min64(total-1, ad.inbuf/r)
	c, d, symb := lookup(v)

	// renormalize
	ad.inbuf -= r * c

	if d < total {
		ad.width = r * (d - c)
	} else {
		ad.width -= r * c
	}

	for ad.width <= quarterInterval {
		ad.width <<= 1
		b, err := ad.reader.ReadBit()
		if err != nil || (b != 0 && b != 1) {
			return 0, err
		}
		ad.inbuf = (ad.inbuf << 1) + uint64(b)
	}
	return symb, nil
}
