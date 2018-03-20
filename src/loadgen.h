#ifndef LOADGEN_H_
#define LOADGEN_H_

class Mesh;
class Object;

typedef Object *(*MeshGenFunc)(void*);

struct MeshGen {
	void *so;
	MeshGenFunc generate;
	char *sofile;
};

MeshGen *load_meshgen(const char *src_fname);
void free_meshgen(MeshGen *mg);

#endif	/* LOADGEN_H_ */
