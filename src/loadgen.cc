#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dlfcn.h>
#include "loadgen.h"
#include "opt.h"

#define CXXFLAGS	"-pedantic -Wall -O3 -fPIC"

MeshGen *load_meshgen(const char *src_fname)
{
	static char cmdline[2048];
	static char so_fname[64];

	sprintf(so_fname, "/tmp/meshgen%d.so", getpid());

	sprintf(cmdline, "%s -shared -o %s " CXXFLAGS " -I" PREFIX "/include/meshgen %s",
			opt.cxx, so_fname, src_fname);
	puts(cmdline);
	if(system(cmdline) != 0) {
		return 0;
	}

	MeshGen *gen = new MeshGen;
	if(!(gen->so = dlopen(so_fname, RTLD_LAZY | RTLD_GLOBAL))) {
		fprintf(stderr, "failed to load shared object: %s: %s\n", so_fname, dlerror());
		delete gen;
		return 0;
	}
	if(!(gen->generate = (MeshGenFunc)dlsym(gen->so, "generate"))) {
		fprintf(stderr, "failed to retrieve entry point from %s\n", so_fname);
		dlclose(gen->so);
		delete gen;
		return 0;
	}

	gen->sofile = new char[strlen(so_fname) + 1];
	strcpy(gen->sofile, so_fname);

	return gen;
}

void free_meshgen(MeshGen *mg)
{
	if(mg) {
		dlclose(mg->so);
		remove(mg->sofile);
		delete [] mg->sofile;
		delete mg;
	}
}
