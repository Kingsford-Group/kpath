#!/bin/sh

VERSION=0.6.1
DIR=kpath-$VERSION
SRC=src/kingsford/kpath

rm -rf $DIR $DIR.tar.gz

mkdir -p $DIR/$SRC
cp *.go $DIR/$SRC

mkdir -p $DIR/$SRC/arithc
cp arithc/*.go $DIR/$SRC/arithc

mkdir -p $DIR/$SRC/bitio
cp bitio/*.go $DIR/$SRC/bitio

cp kpath-$VERSION-macosx $DIR
cp kpath-$VERSION-linux $DIR

cp README.txt $DIR
cp LICENSE.txt $DIR

tar czf $DIR.tar.gz $DIR

