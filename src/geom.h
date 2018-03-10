#ifndef GEOMOBJ_H_
#define GEOMOBJ_H_

#include "gmath/gmath.h"

class GeomObject;
class SceneNode;

struct HitPoint {
	float dist;				//< parametric distance along the ray
	Vec3 pos;			//< position of intersection (orig + dir * dist)
	Vec3 normal;			//< normal at the point of intersection
	const void *obj;		//< pointer to the intersected object
	const SceneNode *node;
	Ray ray;
};

class GeomObject {
public:
	virtual ~GeomObject();

	virtual void set_union(const GeomObject *obj1, const GeomObject *obj2) = 0;
	virtual void set_intersection(const GeomObject *obj1, const GeomObject *obj2) = 0;

	virtual bool intersect(const Ray &ray, HitPoint *hit = 0) const = 0;
};

class Sphere : public GeomObject {
public:
	Vec3 center;
	float radius;

	Sphere();
	Sphere(const Vec3 &center, float radius);

	void set_union(const GeomObject *obj1, const GeomObject *obj2);
	void set_intersection(const GeomObject *obj1, const GeomObject *obj2);

	bool intersect(const Ray &ray, HitPoint *hit = 0) const;
};

class AABox : public GeomObject {
public:
	Vec3 min, max;

	AABox();
	AABox(const Vec3 &min, const Vec3 &max);

	void set_union(const GeomObject *obj1, const GeomObject *obj2);
	void set_intersection(const GeomObject *obj1, const GeomObject *obj2);

	bool intersect(const Ray &ray, HitPoint *hit = 0) const;
};

class Plane : public GeomObject {
public:
	Vec3 pt, normal;

	Plane();
	Plane(const Vec3 &pt, const Vec3 &normal);
	Plane(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3);
	Plane(const Vec3 &normal, float dist);

	void set_union(const GeomObject *obj1, const GeomObject *obj2);
	void set_intersection(const GeomObject *obj1, const GeomObject *obj2);

	bool intersect(const Ray &ray, HitPoint *hit = 0) const;
};

float sphere_distance(const Vec3 &cent, float rad, const Vec3 &pt);
float capsule_distance(const Vec3 &a, float ra, const Vec3 &b, float rb, const Vec3 &pt);

#endif	// GEOMOBJ_H_
