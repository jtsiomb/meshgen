#ifndef OBJECT_H_
#define OBJECT_H_

#include <gmath/gmath.h>

class Mesh;

enum PixelFormat {
	PFMT_RGB, PFMT_RGBA,
	PFMT_RGB_FLOAT, PFMT_RGBA_FLOAT
};

class Image {
public:
	void *pixels;
	int width, height;
	enum PixelFormat fmt;

	Image();
	~Image();
};

class Object {
public:
	Object *next;

	Mat4 xform;
	Mesh *mesh;

	Vec3 color, specular;
	float shininess;

	Image texture;

	Object();
	~Object();
};

bool dump_objlist(Object *objlist, const char *fname);
bool dump_image(Image *img, const char *fname);

#endif	// OBJECT_H_
