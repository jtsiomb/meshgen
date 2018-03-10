#include "mesh.h"
#include "meshgen.h"

extern "C" bool generate(Mesh *mesh, void *cls)
{
	gen_torus(mesh, 1.0, 0.25, 24, 12);

	Mesh *sph = new Mesh;
	gen_sphere(sph, 0.5, 16, 8);

	mesh->append(*sph);
	delete sph;

	return true;
}
