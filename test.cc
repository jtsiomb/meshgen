#include "object.h"
#include "mesh.h"
#include "meshgen.h"

void gen_texture(unsigned char *pixels, int xsz, int ysz);

extern "C" Object *generate(void *cls)
{
	Mat4 xform;

	/* --- torus --- */
	Object *torus = new Object;
	torus->mesh = new Mesh;
	gen_torus(torus->mesh, 1.0, 0.25, 24, 12);
	xform.scale(2, 1, 1);
	torus->mesh->texcoord_apply_xform(xform);

	torus->texture.width = torus->texture.height = 256;
	torus->texture.fmt = PFMT_RGB;
	torus->texture.pixels = new unsigned char[torus->texture.width * torus->texture.height * 3];
	gen_texture((unsigned char*)torus->texture.pixels, torus->texture.width, torus->texture.height);


	/* --- sphere --- */
	Object *sph = new Object;
	sph->mesh = new Mesh;
	gen_sphere(sph->mesh, 0.5, 16, 8);

	sph->xform.translate(0, 1, 0);

	torus->next = sph;

	return torus;
}


void gen_texture(unsigned char *pixels, int xsz, int ysz)
{
	for(int i=0; i<ysz; i++) {
		for(int j=0; j<xsz; j++) {
			*pixels++ = i ^ j;
			*pixels++ = (i ^ j) << 1;
			*pixels++ = (i ^ j) << 2;
		}
	}
}
