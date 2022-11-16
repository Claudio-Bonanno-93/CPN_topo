#!/bin/bash

if [ "$#" -gt 2 ] || [ "$#" -eq 0 ]; then
	echo "Wrong number of arguments! Correct usage: $0 CPN_NCOLORS (TARGET_EXECUTABLE) If no target is given, all targets are compiled"
	exit 1
fi

if [ "$#" -eq 2 ]; then
	CPN_NCOLORS=$1
	TARGET_EXEC=$2
	echo "Compiling $2 with CPN_NCOLORS=$1"
fi

if [ "$#" -eq 1 ]; then
	CPN_NCOLORS=$1
	TARGET_EXEC=''
	echo "Compiling all src files with CPN_NCOLORS=$1"
fi

# change value of N in include/macro.h
aux=$( grep '#define N' include/macro.h )
sed "s/${aux}/#define N ${CPN_NCOLORS}/g" include/macro.h > temp
mv temp include/macro.h

# compile code
./configure CFLAGS='-O3' #CC=icc CFLAGS='-O3 -axCORE-AVX512 -mtune=skylake -ip -ipo' LIBS="-ldl -lz -lc"
make ${TARGET_EXEC} #-j 18
