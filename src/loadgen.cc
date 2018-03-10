#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dlfcn.h>
#include "loadgen.h"
#include "opt.h"

#define SO_FNAME	"/tmp/meshgen_plugin.so"

#define CXXFLAGS	"-pedantic -Wall -g -Isrc "
#define LDFLAGS		"-shared "


MeshGen *load_meshgen(const char *src_fname)
{
	static char cmdline[2048];

	sprintf(cmdline, "%s " CXXFLAGS LDFLAGS "-o " SO_FNAME " %s", opt.cxx, src_fname);
	puts(cmdline);
	if(system(cmdline) != 0) {
		return 0;
	}

	MeshGen *gen = new MeshGen;
	if(!(gen->so = dlopen(SO_FNAME, RTLD_LAZY | RTLD_GLOBAL))) {
		fprintf(stderr, "failed to load shared object: %s: %s\n", SO_FNAME, dlerror());
		delete gen;
		return 0;
	}
	if(!(gen->generate = (MeshGenFunc)dlsym(gen->so, "generate"))) {
		fprintf(stderr, "failed to retrieve entry point from %s\n", SO_FNAME);
		dlclose(gen->so);
		delete gen;
		return 0;
	}

	return gen;
}

void free_meshgen(MeshGen *mg)
{
	if(mg) {
		dlclose(mg->so);
		delete mg;
	}
}
