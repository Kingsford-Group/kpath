#!/bin/sh

PATHSIN="../randgraphs/out.paths"
GRAPHIN="random_zero.edg"

time ./pathenc encode -graph=$GRAPHIN -paths=$PATHSIN -out="test.acp"
time ./pathenc decode -graph=$GRAPHIN -paths="test.acp" -out="recover.txt"
cmp $PATHSIN recover.txt


# Current Times
# Encode: 20s
# Decode: 1m16s
