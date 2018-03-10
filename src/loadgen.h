#ifndef LOADGEN_H_
#define LOADGEN_H_

class Mesh;

typedef bool (*MeshGenFunc)(Mesh*, void*);

struct MeshGen {
	void *so;
	MeshGenFunc generate;
};

MeshGen *load_meshgen(const char *src_fname);
void free_meshgen(MeshGen *mg);

#endif	/* LOADGEN_H_ */
