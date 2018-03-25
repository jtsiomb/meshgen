#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include "object.h"
#include "mesh.h"

#ifdef WIN32
#include <malloc.h>
#include <winsock2.h>
#else
#include <alloca.h>
#include <arpa/inet.h>
#endif

Image::Image()
{
	pixels = 0;
	width = height = 0;
	fmt = PFMT_RGB;
}

Image::~Image()
{
}


Object::Object()
	: color(1, 1, 1), specular(0, 0, 0)
{
	mesh = 0;
	next = 0;
	shininess = 64.0f;
}

Object::~Object()
{
	delete mesh;
}

bool dump_objlist(Object *objlist, const char *fname)
{
	FILE *fp;
	int objidx, voffs = 0;
	Object *obj;

	char *mtlname = (char*)alloca(strlen(fname) + 4);
	strcpy(mtlname, fname);
	char *suffix = strrchr(mtlname, '.');
	if(suffix) *suffix = 0;
	strcat(suffix, ".mtl");

	if(!(fp = fopen(mtlname, "wb"))) {
		fprintf(stderr, "failed to open %s for writing: %s\n", mtlname, strerror(errno));
		return false;
	}
	fprintf(fp, "# meshtool - material for %s\n", fname);

	objidx = 0;
	obj = objlist;
	while(obj) {
		fprintf(fp, "newmtl mtl%03d\n", objidx);
		fprintf(fp, "Ka 0 0 0\n");
		fprintf(fp, "Kd %g %g %g\n", obj->color.x, obj->color.y, obj->color.z);
		fprintf(fp, "Ks %g %g %g\n", obj->specular.x, obj->specular.y, obj->specular.z);
		fprintf(fp, "Ns %g\n", 100.0f * obj->shininess / 128.0f);
		fprintf(fp, "d 1\n");
		fprintf(fp, "illum 2\n");

		if(obj->texture.pixels) {
			char texname[32];
			sprintf(texname, "tex%03d.ppm", objidx);
			if(dump_image(&obj->texture, texname)) {
				fprintf(fp, "map_Kd %s\n", texname);
			}
		}

		fputc('\n', fp);
		obj = obj->next;
		++objidx;
	}
	fclose(fp);


	if(!(fp = fopen(fname, "wb"))) {
		fprintf(stderr, "failed to open %s for writing: %s\n", fname, strerror(errno));
		return false;
	}
	fprintf(fp, "mtllib %s\n", mtlname);

	objidx = 0;
	obj = objlist;
	while(obj) {
		if(obj != objlist) {
			fprintf(fp, "\n# --------------------------------\n\n");
		}
		fprintf(fp, "o mesh%03d\n", objidx);
		fprintf(fp, "usemtl mtl%03d\n", objidx);

		for(int i=0; i<3; i++) {
			fprintf(fp, "# xform%d %g %g %g %g\n", i, obj->xform[0][i], obj->xform[1][i],
					obj->xform[2][i], obj->xform[3][i]);
		}

		Mesh m = *obj->mesh;
		m.apply_xform(obj->xform);
		printf("vertex offset: %d\n", voffs);
		m.dump_obj(fp, voffs);

		voffs += m.get_attrib_count(MESH_ATTR_VERTEX);

		obj = obj->next;
		++objidx;
	}
	fclose(fp);

	return true;
}

bool dump_image(Image *img, const char *fname)
{
	FILE *fp;
	unsigned char *ptr = (unsigned char*)img->pixels;
	float *fptr = (float*)img->pixels;
	uint16_t pix16[3];

	if(!(fp = fopen(fname, "wb"))) {
		fprintf(stderr, "failed to open %s for writing: %s\n", fname, strerror(errno));
		return false;
	}

	int max_val = img->fmt == PFMT_RGB_FLOAT || img->fmt == PFMT_RGBA_FLOAT ? 65535 : 255;

	fprintf(fp, "P6\n%d %d\n%d\n", img->width, img->height, max_val);
	for(int i=0; i<img->width * img->height; i++) {
		switch(img->fmt) {
		case PFMT_RGB:
		case PFMT_RGBA:
			fputc(ptr[0], fp);
			fputc(ptr[1], fp);
			fputc(ptr[2], fp);
			ptr += img->fmt == PFMT_RGB ? 3 : 4;
			break;

		case PFMT_RGB_FLOAT:
		case PFMT_RGBA_FLOAT:
			for(int j=0; j<3; j++) {
				pix16[j] = htons((int)(*fptr++ * 65535.0f));
			}
			fwrite(pix16, 2, 3, fp);
			if(img->fmt == PFMT_RGBA_FLOAT) ++fptr;
			break;
		}
	}

	fclose(fp);
	return true;
}
