#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "loadgen.h"
#include "opt.h"

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <dlfcn.h>
#endif

#define CXXFLAGS	"-pedantic -Wall -Wno-unused-function -O3 -fPIC -include object.h -include mesh.h -include meshgen.h"

MeshGen *load_meshgen(const char *src_fname)
{
	static char cmdline[2048];
	static char so_fname[64];

#ifdef WIN32
	strcpy(so_fname, "meshgen.dll");
	sprintf(cmdline, "%s -shared -o %s " CXXFLAGS " -I" PREFIX "/include/meshgen %s",
			opt.cxx, so_fname, src_fname);
#else
	sprintf(so_fname, "/tmp/meshgen%d.so", getpid());
	sprintf(cmdline, "%s -shared -o %s " CXXFLAGS " -I" PREFIX "/include/meshgen %s",
			opt.cxx, so_fname, src_fname);
#endif
	puts(cmdline);
	if(system(cmdline) != 0) {
		return 0;
	}

	MeshGen *gen = new MeshGen;
#ifdef WIN32
	if(!(gen->so = LoadLibrary(so_fname))) {
		fprintf(stderr, "failed to load dll: %s\n", so_fname);
		delete gen;
		return 0;
	}
	if(!(gen->generate = (MeshGenFunc)GetProcAddress((HMODULE)gen->so, "generate"))) {
		fprintf(stderr, "failed to retrieve entry point from %s\n", so_fname);
		CloseHandle(gen->so);
		delete gen;
		return 0;
	}
#else	/* UNIX */
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
#endif

	gen->sofile = new char[strlen(so_fname) + 1];
	strcpy(gen->sofile, so_fname);

	return gen;
}

void free_meshgen(MeshGen *mg)
{
	if(mg) {
#ifdef WIN32
		CloseHandle(mg->so);
#else
		dlclose(mg->so);
#endif
		remove(mg->sofile);
		delete [] mg->sofile;
		delete mg;
	}
}
