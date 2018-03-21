#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include "geom.h"

#define SPHERE(ptr)	((Sphere*)ptr)
#define AABOX(ptr)	((AABox*)ptr)
#define BOX(ptr)	((Box*)ptr)
#define PLANE(ptr)	((Plane*)ptr)

GeomObject::GeomObject()
{
	type = GOBJ_UNKNOWN;
}

GeomObject::~GeomObject()
{
}

bool GeomObject::valid() const
{
	return true;
}

void GeomObject::invalidate()
{
}

float GeomObject::distance(const Vec3 &v) const
{
	return sqrt(distance_sq(v));
}

float GeomObject::signed_distance(const Vec3 &v) const
{
	return sqrt(signed_distance_sq(v));
}

Sphere::Sphere()
{
	type = GOBJ_SPHERE;
	radius = 1.0;
}

Sphere::Sphere(const Vec3 &cent, float radius)
	: center(cent)
{
	type = GOBJ_SPHERE;
	this->radius = radius;
}

bool Sphere::valid() const
{
	return radius >= 0.0f;
}

void Sphere::invalidate()
{
	center = Vec3(0, 0, 0);
	radius = -1;
}

bool Sphere::intersect(const Ray &ray, HitPoint *hit) const
{
	float a = dot(ray.dir, ray.dir);
	float b = 2.0 * ray.dir.x * (ray.origin.x - center.x) +
		2.0 * ray.dir.y * (ray.origin.y - center.y) +
		2.0 * ray.dir.z * (ray.origin.z - center.z);
	float c = dot(ray.origin, ray.origin) + dot(center, center) -
		2.0 * dot(ray.origin, center) - radius * radius;

	float discr = b * b - 4.0 * a * c;
	if(discr < 1e-4) {
		return false;
	}

	float sqrt_discr = sqrt(discr);
	float t0 = (-b + sqrt_discr) / (2.0 * a);
	float t1 = (-b - sqrt_discr) / (2.0 * a);

	if(t0 < 1e-4)
		t0 = t1;
	if(t1 < 1e-4)
		t1 = t0;

	float t = t0 < t1 ? t0 : t1;
	if(t < 1e-4) {
		return false;
	}

	// fill the HitPoint structure
	if(hit) {
		hit->obj = this;
		hit->dist = t;
		hit->pos = ray.origin + ray.dir * t;
		hit->normal = (hit->pos - center) / radius;
	}
	return true;
}

bool Sphere::contains(const Vec3 &pt) const
{
	return length_sq(pt - center) <= radius * radius;
}

float Sphere::distance_sq(const Vec3 &v) const
{
	return std::max(length_sq(v - center) - radius * radius, 0.0f);
}

float Sphere::signed_distance_sq(const Vec3 &v) const
{
	return length_sq(v - center) - radius * radius;
}

AABox::AABox()
{
	type = GOBJ_AABOX;
}

AABox::AABox(const Vec3 &vmin, const Vec3 &vmax)
	: min(vmin), max(vmax)
{
	type = GOBJ_AABOX;
}

bool AABox::valid() const
{
	return min.x <= max.x && min.y <= max.y && min.z <= max.z;
}

void AABox::invalidate()
{
	min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
}

Vec3 AABox::get_corner(int idx) const
{
	Vec3 v[] = {min, max};
	static const int xidx[] = {0, 1, 1, 0, 0, 1, 1, 0};
	static const int yidx[] = {0, 0, 0, 0, 1, 1, 1, 1};
	static const int zidx[] = {0, 0, 1, 1, 0, 0, 1, 1};
	return Vec3(v[xidx[idx]].x, v[yidx[idx]].y, v[zidx[idx]].z);
}

Plane AABox::get_plane(int pidx) const
{
	Vec3 c = (max - min) * 0.5f;
	switch(pidx) {
	case AABOX_PLANE_PX:
		return Plane(Vec3(max.x, c.y, c.z), Vec3(1, 0, 0));
	case AABOX_PLANE_NX:
		return Plane(Vec3(min.x, c.y, c.z), Vec3(-1, 0, 0));
	case AABOX_PLANE_PY:
		return Plane(Vec3(c.x, max.x, c.z), Vec3(0, 1, 0));
	case AABOX_PLANE_NY:
		return Plane(Vec3(c.x, min.x, c.z), Vec3(0, -1, 0));
	case AABOX_PLANE_PZ:
		return Plane(Vec3(c.x, c.y, max.z), Vec3(0, 0, 1));
	case AABOX_PLANE_NZ:
		return Plane(Vec3(c.x, c.y, min.z), Vec3(0, 0, -1));
	default:
		break;
	}
	abort();
	return Plane();
}

bool AABox::intersect(const Ray &ray, HitPoint *hit) const
{
	Vec3 param[2] = {min, max};
	Vec3 inv_dir(1.0 / ray.dir.x, 1.0 / ray.dir.y, 1.0 / ray.dir.z);
	int sign[3] = {inv_dir.x < 0, inv_dir.y < 0, inv_dir.z < 0};

	float tmin = (param[sign[0]].x - ray.origin.x) * inv_dir.x;
	float tmax = (param[1 - sign[0]].x - ray.origin.x) * inv_dir.x;
	float tymin = (param[sign[1]].y - ray.origin.y) * inv_dir.y;
	float tymax = (param[1 - sign[1]].y - ray.origin.y) * inv_dir.y;

	if(tmin > tymax || tymin > tmax) {
		return false;
	}
	if(tymin > tmin) {
		tmin = tymin;
	}
	if(tymax < tmax) {
		tmax = tymax;
	}

	float tzmin = (param[sign[2]].z - ray.origin.z) * inv_dir.z;
	float tzmax = (param[1 - sign[2]].z - ray.origin.z) * inv_dir.z;

	if(tmin > tzmax || tzmin > tmax) {
		return false;
	}
	if(tzmin > tmin) {
		tmin = tzmin;
	}
	if(tzmax < tmax) {
		tmax = tzmax;
	}

	float t = tmin < 1e-4 ? tmax : tmin;
	if(t >= 1e-4) {

		if(hit) {
			hit->obj = this;
			hit->dist = t;
			hit->pos = ray.origin + ray.dir * t;

			float min_dist = FLT_MAX;
			Vec3 offs = min + (max - min) / 2.0;
			Vec3 local_hit = hit->pos - offs;

			static const Vec3 axis[] = {
				Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)
			};
			//int tcidx[][2] = {{2, 1}, {0, 2}, {0, 1}};

			for(int i=0; i<3; i++) {
				float dist = fabs((max[i] - offs[i]) - fabs(local_hit[i]));
				if(dist < min_dist) {
					min_dist = dist;
					hit->normal = axis[i] * (local_hit[i] < 0.0 ? 1.0 : -1.0);
					//hit->texcoord = Vec2(hit->pos[tcidx[i][0]], hit->pos[tcidx[i][1]]);
				}
			}
		}
		return true;
	}
	return false;
}

bool AABox::contains(const Vec3 &v) const
{
	return v.x >= min.x && v.y >= min.y && v.z >= min.z &&
		v.x <= max.x && v.y <= max.y && v.z <= max.z;
}

#define SQ(x)	((x) * (x))
float AABox::distance_sq(const Vec3 &v) const
{
	float dsq = 0.0f;

	for(int i=0; i<3; i++) {
		if(v[i] < min[i]) dsq += SQ(min[i] - v[i]);
		if(v[i] > max[i]) dsq += SQ(v[i] - max[i]);
	}
	return dsq;
}

float AABox::signed_distance_sq(const Vec3 &v) const
{
	if(!contains(v)) {
		return distance_sq(v);
	}

	float dsq = 0.0f;
	for(int i=0; i<3; i++) {
		float dmin = v[i] - min[i];
		float dmax = max[i] - v[i];
		dsq -= dmin < dmax ? SQ(dmin) : SQ(dmax);
	}
	return dsq;
}

Box::Box()
{
	type = GOBJ_BOX;
}

Box::Box(const AABox &aabox, const Mat4 &xform)
	: xform(xform)
{
	type = GOBJ_BOX;
	min = aabox.min;
	max = aabox.max;
}

void Box::invalidate()
{
	AABox::invalidate();
	xform = Mat4::identity;
}

Box::Box(const Vec3 &min, const Vec3 &max)
	: AABox(min, max)
{
	type = GOBJ_BOX;
}

Box::Box(const Vec3 &min, const Vec3 &max, const Mat4 &xform)
	: AABox(min, max), xform(xform)
{
	type = GOBJ_BOX;
}

// XXX all this shit is completely untested
Box::Box(const Vec3 &pos, const Vec3 &vi, const Vec3 &vj, const Vec3 &vk)
{
	type = GOBJ_BOX;
	float ilen = length(vi);
	float jlen = length(vj);
	float klen = length(vk);

	min = Vec3(-ilen, -jlen, -klen);
	max = Vec3(ilen, jlen, klen);

	float si = ilen == 0.0 ? 1.0 : 1.0 / ilen;
	float sj = jlen == 0.0 ? 1.0 : 1.0 / jlen;
	float sk = klen == 0.0 ? 1.0 : 1.0 / klen;

	xform = Mat4(vi * si, vj * sj, vk * sk);
	xform.translate(pos);
}

Box::Box(const Vec3 *varr, int vcount)
{
	type = GOBJ_BOX;
	calc_bounding_aabox(this, varr, vcount);
}

Vec3 Box::get_corner(int idx) const
{
	return xform * AABox::get_corner(idx);
}

bool Box::intersect(const Ray &ray, HitPoint *hit) const
{
	Mat4 inv_xform = inverse(xform);
	Mat4 dir_inv_xform = inv_xform.upper3x3();
	Mat4 dir_xform = transpose(dir_inv_xform);
	Ray local_ray = Ray(inv_xform * ray.origin, dir_inv_xform * ray.dir);

	bool res = AABox::intersect(local_ray, hit);
	if(!res || !hit) return res;

	hit->pos = xform * hit->pos;
	hit->normal = dir_xform * hit->normal;
	hit->local_ray = local_ray;
	hit->ray = ray;
	return true;
}

bool Box::contains(const Vec3 &pt) const
{
	// XXX is it faster to extract 6 planes and do dot products? sounds marginal
	return AABox::contains(inverse(xform) * pt);
}

float Box::distance_sq(const Vec3 &v) const
{
	return 0.0f;	// TODO
}

float Box::signed_distance_sq(const Vec3 &v) const
{
	return 0.0f;	// TODO
}

Plane::Plane()
	: normal(0.0, 1.0, 0.0)
{
	type = GOBJ_PLANE;
}

Plane::Plane(const Vec3 &p, const Vec3 &norm)
	: pt(p)
{
	type = GOBJ_PLANE;
	normal = normalize(norm);
}

Plane::Plane(const Vec3 &p1, const Vec3 &p2, const Vec3 &p3)
	: pt(p1)
{
	type = GOBJ_PLANE;
	normal = normalize(cross(p2 - p1, p3 - p1));
}

Plane::Plane(const Vec3 &normal, float dist)
{
	type = GOBJ_PLANE;
	this->normal = normalize(normal);
	pt = this->normal * dist;
}

bool Plane::intersect(const Ray &ray, HitPoint *hit) const
{
	float ndotdir = dot(normal, ray.dir);
	if(fabs(ndotdir) < 1e-4) {
		return false;
	}

	if(hit) {
		Vec3 ptdir = pt - ray.origin;
		float t = dot(normal, ptdir) / ndotdir;

		hit->dist = t;
		hit->pos = ray.origin + ray.dir * t;
		hit->normal = normal;
		hit->obj = this;
	}
	return true;
}

bool Plane::contains(const Vec3 &v) const
{
	return dot(v, normal) <= 0.0;
}

float Plane::distance(const Vec3 &v) const
{
	return std::max(dot(v - pt, normal), 0.0f);
}

float Plane::signed_distance(const Vec3 &v) const
{
	return dot(v - pt, normal);
}

float Plane::distance_sq(const Vec3 &v) const
{
	float d = distance(v);
	return d * d;
}

float Plane::signed_distance_sq(const Vec3 &v) const
{
	float d = distance(v);
	return dot(v, normal) >= 0.0f ? d * d : -(d * d);
}


Disc::Disc()
{
	type = GOBJ_DISC;
	radius = 1.0;
}

Disc::Disc(const Vec3 &pt, const Vec3 &normal, float rad)
	: Plane(pt, normal)
{
	type = GOBJ_DISC;
	radius = rad;
}

Disc::Disc(const Vec3 &normal, float dist, float rad)
	: Plane(normal, dist)
{
	type = GOBJ_DISC;
	radius = rad;
}

bool Disc::valid() const
{
	return radius >= 0.0f;
}

void Disc::invalidate()
{
	radius = -1;
}

bool Disc::intersect(const Ray &ray, HitPoint *hit) const
{
	HitPoint phit;
	if(Plane::intersect(ray, &phit)) {
		if(length_sq(phit.pos - pt) <= radius * radius) {
			*hit = phit;
			return true;
		}
	}
	return false;
}

bool Disc::contains(const Vec3 &pt) const
{
	Vec3 pj = proj_point_plane(pt, *this);
	return length_sq(pj - this->pt) <= radius * radius;
}

float Disc::distance_sq(const Vec3 &v) const
{
	return 0.0;	// TODO
}

float Disc::signed_distance_sq(const Vec3 &v) const
{
	return 0.0;	// TODO
}


Vec3 proj_point_plane(const Vec3 &pt, const Plane &plane)
{
	float dist = plane.signed_distance(pt);
	return pt - plane.normal * dist;
}

// ---- bounding sphere calculations ----

bool calc_bounding_sphere(Sphere *sph, const GeomObject *obj)
{
	if(!obj->valid()) {
		sph->invalidate();
		return true;
	}

	switch(obj->type) {
	case GOBJ_SPHERE:
		*sph = *(Sphere*)obj;
		break;

	case GOBJ_AABOX:
		sph->center = (AABOX(obj)->min + AABOX(obj)->max) * 0.5;
		sph->radius = length(AABOX(obj)->max - AABOX(obj)->min) * 0.5;
		break;

	case GOBJ_BOX:
		sph->center = (BOX(obj)->min + BOX(obj)->max) * 0.5 + BOX(obj)->xform.get_translation();
		sph->radius = length(BOX(obj)->max - BOX(obj)->min) * 0.5;
		break;

	case GOBJ_PLANE:
	default:
		return false;
	}
	return true;
}

bool calc_bounding_sphere(Sphere *sph, const GeomObject *a, const GeomObject *b)
{
	Sphere bsa, bsb;

	if(!calc_bounding_sphere(&bsa, a) || !calc_bounding_sphere(&bsb, b)) {
		return false;
	}

	float dist = length(bsa.center - bsb.center);
	float surf_dist = dist - (bsa.radius + bsb.radius);
	float d1 = bsa.radius + surf_dist / 2.0;
	float d2 = bsb.radius + surf_dist / 2.0;
	float t = d1 / (d1 + d2);

	if(t < 0.0) t = 0.0;
	if(t > 1.0) t = 1.0;

	sph->center = bsa.center * t + bsb.center * (1.0 - t);
	sph->radius = std::max(dist * t + bsb.radius, dist * (1.0f - t) + bsa.radius);
	return true;
}

bool calc_bounding_sphere(Sphere *sph, const GeomObject **objv, int num)
{
	if(num <= 0) return false;

	if(!calc_bounding_sphere(sph, objv[0])) {
		return false;
	}

	for(int i=1; i<num; i++) {
		if(!calc_bounding_sphere(sph, sph, objv[i])) {
			return false;
		}
	}
	return true;
}

bool calc_bounding_sphere(Sphere *sph, const Vec3 *v, int num, const Mat4 &xform)
{
	if(num <= 0) return false;

	sph->center = Vec3(0.0, 0.0, 0.0);
	for(int i=0; i<num; i++) {
		sph->center += xform * v[i];
	}
	sph->center /= (float)num;

	float rad_sq = 0.0f;
	for(int i=0; i<num; i++) {
		Vec3 dir = xform * v[i] - sph->center;
		rad_sq = std::max(rad_sq, dot(dir, dir));
	}
	sph->radius = sqrt(rad_sq);
	return true;
}

bool calc_bounding_aabox(AABox *box, const GeomObject *obj)
{
	if(!obj->valid()) {
		box->invalidate();
		return true;
	}

	switch(obj->type) {
	case GOBJ_AABOX:
		*box = *(AABox*)obj;
		break;

	case GOBJ_BOX:
		{
			Vec3 v[8];
			for(int i=0; i<8; i++) {
				v[i] = BOX(obj)->get_corner(i);
			}
			calc_bounding_aabox(box, v, 8);
		}
		break;

	case GOBJ_SPHERE:
		{
			float r = SPHERE(obj)->radius;
			box->min = SPHERE(obj)->center - Vec3(r, r, r);
			box->max = SPHERE(obj)->center + Vec3(r, r, r);
		}
		break;

	case GOBJ_PLANE:
	default:
		return false;
	}
	return true;
}

bool calc_bounding_aabox(AABox *box, const GeomObject *a, const GeomObject *b)
{
	AABox bba, bbb;

	if(!calc_bounding_aabox(&bba, a) || !calc_bounding_aabox(&bbb, b)) {
		return false;
	}

	for(int i=0; i<3; i++) {
		box->min[i] = std::min(bba.min[i], bbb.min[i]);
		box->max[i] = std::max(bba.max[i], bbb.max[i]);
	}
	return true;
}

bool calc_bounding_aabox(AABox *box, const GeomObject **objv, int num)
{
	if(num <= 0) return false;

	if(!calc_bounding_aabox(box, objv[0])) {
		return false;
	}

	for(int i=1; i<num; i++) {
		if(!calc_bounding_aabox(box, box, objv[i])) {
			return false;
		}
	}
	return true;
}

bool calc_bounding_aabox(AABox *box, const Vec3 *v, int num, const Mat4 &xform)
{
	if(num <= 0) return false;

	box->min = box->max = xform * v[0];
	for(int i=1; i<num; i++) {
		Vec3 p = xform * v[i];

		for(int j=0; j<3; j++) {
			box->min[j] = std::min(box->min[j], p[j]);
			box->max[j] = std::max(box->max[j], p[j]);
		}
	}
	return true;
}

bool calc_bounding_box(Box *box, const GeomObject *obj)
{
	if(!obj->valid()) {
		box->invalidate();
		return true;
	}

	switch(obj->type) {
	case GOBJ_BOX:
		*box = *(Box*)obj;
		break;

	case GOBJ_AABOX:
		box->min = BOX(obj)->min;
		box->max = BOX(obj)->max;
		box->xform = Mat4::identity;
		break;

	case GOBJ_SPHERE:
		{
			float r = SPHERE(obj)->radius;
			box->min = SPHERE(obj)->center - Vec3(r, r, r);
			box->max = SPHERE(obj)->center + Vec3(r, r, r);
			box->xform = Mat4::identity;
		}
		break;

	case GOBJ_PLANE:
	default:
		return false;
	}
	return true;
}

bool intersect_sphere_sphere(Disc *result, const Sphere &a, const Sphere &b)
{
	Vec3 dir = b.center - a.center;

	float dist_sq = length_sq(dir);
	if(dist_sq <= 1e-8) return false;

	float rsum = a.radius + b.radius;
	float rdif = fabs(a.radius - b.radius);
	if(dist_sq > rsum * rsum || dist_sq < rdif * rdif) {
		return false;
	}

	float dist = sqrt(dist_sq);
	float t = (dist_sq + a.radius * a.radius - b.radius * b.radius) / (2.0 * sqrt(dist_sq));

	result->pt = a.center + dir * t;
	result->normal = dir / dist;
	result->radius = sin(acos(t)) * a.radius;
	return true;
}

bool intersect_plane_plane(Ray *result, const Plane &a, const Plane &b)
{
	return false;	// TODO
}

bool intersect_sphere_plane(Sphere *result, const Sphere &s, const Plane &p)
{
	return false;	// TODO
}

bool intersect_plane_sphere(Sphere *result, const Plane &p, const Sphere &s)
{
	return false;	// TODO
}

bool intersect_aabox_aabox(AABox *res, const AABox &a, const AABox &b)
{
	for(int i=0; i<3; i++) {
		res->min[i] = std::max(a.min[i], b.min[i]);
		res->max[i] = std::min(a.max[i], b.max[i]);

		if(res->max[i] < res->min[i]) {
			res->max[i] = res->min[i];
		}
	}
	return res->min.x != res->max.x && res->min.y != res->max.y && res->min.z != res->max.z;
}

bool collision_sphere_aabox(const Sphere &s, const AABox &b)
{
	return b.distance_sq(s.center) <= s.radius * s.radius;
}
