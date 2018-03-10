#include <stdio.h>
#include "opt.h"

Options opt;

bool parse_options(int argc, char **argv)
{
	opt.cxx = (char*)"c++";

	for(int i=1; i<argc; i++) {
		if(opt.fname) {
			fprintf(stderr, "unexpected argument: %s\n", argv[i]);
			return false;
		}
		opt.fname = argv[i];
	}

	if(!opt.fname) {
		fprintf(stderr, "you need to pass the filename of a meshgen source file\n");
		return false;
	}
	return true;
}
