#!/bin/sh

VERSION=0.6.1
DIR=kpath-$VERSION

rm -rf $DIR $DIR.tar.gz

mkdir $DIR
mkdir $DIR/arithc
mkdir $DIR/bitio

cp *.go $DIR
cp arithc/*.go $DIR/arithc
cp bitio/*.go $DIR/bitio
cp kpath-$VERSION-macosx $DIR
cp kpath-$VERSION-linux $DIR

cp README.txt $DIR
cp LICENSE.txt $DIR

tar czf $DIR.tar.gz $DIR

