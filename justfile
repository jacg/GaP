# -*-Makefile-*-

# (Re)compile and run
run:
	#!/usr/bin/env sh
	cmake -S . -B build &&
	cmake --build build -j &&
	build/GaP

# Remove all traces of the local copy of the example
clean:
	rm build -rf
