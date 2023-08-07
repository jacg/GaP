# -*-Makefile-*-

# (Re)compile and run
run ARG1='':
	#!/usr/bin/env sh
	cmake -S . -B build &&
	cmake --build build -j &&
	build/GaP {{ARG1}}

view:
	nix run --impure github:guibou/nixGL -- just
	
# Remove all traces of the local copy of the example
clean:
	rm build -rf

