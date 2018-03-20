#include "object.h"
#include "mesh.h"
#include "meshgen.h"

#define TEXSZ	256

// The generate function must be declared as extern "C" for meshtool to
// be able to find it.
extern "C" Object *generate(void *cls)
{
	Mat4 xform;

	/* --- torus --- */
	Object *torus = new Object;
	torus->mesh = new Mesh;
	gen_torus(torus->mesh, 1.0, 0.25, 24, 12);
	// also scale the horizontal texcoords by 2
	xform.scale(2, 1, 1);
	torus->mesh->texcoord_apply_xform(xform);

	unsigned char *pixels = new unsigned char[TEXSZ * TEXSZ * 3];
	torus->texture.width = torus->texture.height = TEXSZ;
	torus->texture.fmt = PFMT_RGB;
	torus->texture.pixels = pixels;

	for(int i=0; i<TEXSZ; i++) {
		for(int j=0; j<TEXSZ; j++) {
			*pixels++ = i ^ j;
			*pixels++ = (i ^ j) << 1;
			*pixels++ = (i ^ j) << 2;
		}
	}


	/* --- sphere --- */
	Object *sph = new Object;
	sph->mesh = new Mesh;
	gen_sphere(sph->mesh, 0.5, 16, 8);
	sph->xform.translate(0, 1, 0);

	torus->next = sph;	// link sph after torus

	return torus;	// return pointer to the first node of the list
}
