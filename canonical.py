#!/usr/bin/env python

import sys

def RC(x):
    rcmap = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    A = []
    for x in reversed(x):
        A.append(rcmap[x])
    return "".join(A) 

with open(sys.argv[1]) as inp:
    for line in inp:
        line = line.strip()
        rcr = RC(line)
        print line if line < rcr else rcr
