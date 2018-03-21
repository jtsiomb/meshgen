#ifndef GEOMOBJ_H_
#define GEOMOBJ_H_

/* TODO:
 * - implement distance functions
 */

#include <gmath/gmath.h>

enum GeomObjectType {
	GOBJ_UNKNOWN,
	GOBJ_SPHERE,
	GOBJ_AABOX,
	GOBJ_BOX,
	GOBJ_PLANE,
	GOBJ_DISC
};

enum {
	AABOX_PLANE_PX,
	AABOX_PLANE_NX,
	AABOX_PLANE_PY,
	AABOX_PLANE_NY,
	AABOX_PLANE_PZ,
	AABOX_PLANE_NZ
};

class GeomObject;
class Plane;

struct HitPoint {
	float dist;			// parametric distance along the ray
	Vec3 pos;			// position of intersection (orig + dir * dist)
	Vec3 normal;		// normal at the point of intersection
	Ray ray, local_ray;
	const GeomObject *obj;	// pointer to the intersected geom-object
	void *data;			// place to hang extra data
};

class GeomObject {
public:
	GeomObjectType type;

	GeomObject();
	virtual ~GeomObject();

	virtual bool valid() const;
	virtual void invalidate();

	virtual bool intersect(const Ray &ray, HitPoint *hit = 0) const = 0;
	virtual bool contains(const Vec3 &pt) const = 0;

	virtual float distance(const Vec3 &v) const;
	virtual float signed_distance(const Vec3 &v) const;

	virtual float distance_sq(const Vec3 &v) const = 0;
	virtual float signed_distance_sq(const Vec3 &v) const = 0;
};

class Sphere : public GeomObject {
public:
	Vec3 center;
	float radius;

	Sphere();
	Sphere(const Vec3 &center, float radius);

	virtual bool valid() const;
	virtual void invalidate();

	virtual bool intersect(const Ray &ray, HitPoint *hit = 0) const;
	virtual bool contains(const Vec3 &pt) const;

	virtual float distance_sq(const Vec3 &v) const;
	virtual float signed_distance_sq(const Vec3 &v) const;
};

class AABox : public GeomObject {
public:
	Vec3 min, max;

	AABox();
	AABox(const Vec3 &min, const Vec3 &max);

	virtual bool valid() const;
	virtual void invalidate();

	virtual Vec3 get_corner(int idx) const;
	virtual Plane get_plane(int pidx) const;

	virtual bool intersect(const Ray &ray, HitPoint *hit = 0) const;
	virtual bool contains(const Vec3 &pt) const;

	virtual float distance_sq(const Vec3 &v) const;
	virtual float signed_distance_sq(const Vec3 &v) const;
};

class Box : public AABox {
public:
	Mat4 xform;

	Box();
	Box(const AABox &aabox, const Mat4 &xform);
	Box(const Vec3 &min, const Vec3 &max);
	Box(const Vec3 &min, const Vec3 &max, const Mat4 &xform);
	Box(const Vec3 &pos, const Vec3 &vi, const Vec3 &vj, const Vec3 &vk);
	Box(const Vec3 *varr, int vcount);

	virtual void invalidate();

	virtual Vec3 get_corner(int idx) const;

	virtual bool intersect(const Ray &ray, HitPoint *hit = 0) const;
	virtual bool contains(const Vec3 &pt) const;

	virtual float distance_sq(const Vec3 &v) const;
	virtual float signed_distance_sq(const Vec3 &v) const;
};


class Plane : public GeomObject {
public:
	Vec3 pt, normal;

	Plane();
	Plane(const Vec3 &pt, const Vec3 &normal);
	Plane(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3);
	Plane(const Vec3 &normal, float dist);

	virtual bool intersect(const Ray &ray, HitPoint *hit = 0) const;
	virtual bool contains(const Vec3 &pt) const;

	virtual float distance(const Vec3 &v) const;
	virtual float signed_distance(const Vec3 &v) const;

	virtual float distance_sq(const Vec3 &v) const;
	virtual float signed_distance_sq(const Vec3 &v) const;
};

class Disc : public Plane {
public:
	float radius;

	Disc();
	Disc(const Vec3 &pt, const Vec3 &normal, float rad);
	Disc(const Vec3 &normal, float dist, float rad);

	virtual bool valid() const;
	virtual void invalidate();

	virtual bool intersect(const Ray &ray, HitPoint *hit = 0) const;
	//! true if the projection of pt to the plane is contained within the disc radius
	virtual bool contains(const Vec3 &pt) const;

	virtual float distance_sq(const Vec3 &v) const;
	virtual float signed_distance_sq(const Vec3 &v) const;
};

//! project point to plane
Vec3 proj_point_plane(const Vec3 &pt, const Plane &plane);

//! calculate the bounding sphere of any object
bool calc_bounding_sphere(Sphere *sph, const GeomObject *obj);
//! calculate the bounding sphere of two objects
bool calc_bounding_sphere(Sphere *sph, const GeomObject *a, const GeomObject *b);
//! calculate the bounding sphere of multiple objects
bool calc_bounding_sphere(Sphere *sph, const GeomObject **objv, int num);
//! calculate the bounding sphere of multiple points, optionally transformed by `xform`
bool calc_bounding_sphere(Sphere *sph, const Vec3 *v, int num, const Mat4 &xform = Mat4::identity);

//! calculate the bounding axis-aligned box of any object
bool calc_bounding_aabox(AABox *box, const GeomObject *obj);
//! calculate the bounding axis-aligned box of two objects
bool calc_bounding_aabox(AABox *box, const GeomObject *a, const GeomObject *b);
//! calculate the bounding axis-aligned box of multiple objects
bool calc_bounding_aabox(AABox *box, const GeomObject **objv, int num);
//! calculate the bounding axis-aligned box of multiple points, optionally transformed by `xform`
bool calc_bounding_aabox(AABox *box, const Vec3 *v, int num, const Mat4 &xform = Mat4::identity);

//! calculate the bounding box of any object
bool calc_bounding_box(Box *box, const GeomObject *obj);

//! calculate the intersection plane of two spheres. false if there is no plane or infinite planes.
bool intersect_sphere_sphere(Plane *result, const Sphere &a, const Sphere &b);
//! calculate the intersection line of two planes. returns false if planes are exactly parallel.
bool intersect_plane_plane(Ray *result, const Plane &a, const Plane &b);
/*! calculate the intesection circle of a plane and a sphere. returns false if they don't intersect.
 * \{
 */
bool intersect_sphere_plane(Sphere *result, const Sphere &s, const Plane &p);
bool intersect_plane_sphere(Sphere *result, const Plane &p, const Sphere &s);
//! \}

//! calculate the intersection of two axis-aligned bounding boxes
bool intersect_aabox_aabox(AABox *res, const AABox &a, const AABox &b);

//! determine if a sphere and an axis-aligned box collide
bool collision_sphere_aabox(const Sphere &s, const AABox &b);

#endif	// GEOMOBJ_H_
