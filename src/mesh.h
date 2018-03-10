#ifndef MESH_H_
#define MESH_H_

#include <stdio.h>
#include <string>
#include <vector>
#include "gmath/gmath.h"
#include "geom.h"

enum {
	MESH_ATTR_VERTEX,
	MESH_ATTR_NORMAL,
	MESH_ATTR_TANGENT,
	MESH_ATTR_TEXCOORD,
	MESH_ATTR_COLOR,
	MESH_ATTR_BONEWEIGHTS,
	MESH_ATTR_BONEIDX,

	NUM_MESH_ATTR
};

// intersection mode flags
enum {
	ISECT_DEFAULT	= 0,	// default (whole mesh, all intersections)
	ISECT_FRONT		= 1,	// front-faces only
	ISECT_FACE		= 2,	// return intersected face pointer instead of mesh
	ISECT_VERTICES	= 4		// return (?) TODO
};

//class XFormNode;


class Triangle {
public:
	Vec3 v[3];
	Vec3 normal;
	bool normal_valid;
	int id;

	Triangle();
	Triangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);
	Triangle(int n, const Vec3 *varr, const unsigned int *idxarr = 0);

	/// calculate normal (quite expensive)
	void calc_normal();
	const Vec3 &get_normal() const;

	void transform(const Mat4 &xform);

	void draw() const;
	void draw_wire() const;

	/// calculate barycentric coordinates of a point
	Vec3 calc_barycentric(const Vec3 &pos) const;

	bool intersect(const Ray &ray, HitPoint *hit = 0) const;
};


class Mesh {
private:
	std::string name;
	unsigned int nverts, nfaces;

	// current value for each attribute for the immedate mode
	// interface.
	Vec4 cur_val[NUM_MESH_ATTR];

	unsigned int buffer_objects[NUM_MESH_ATTR + 1];

	// vertex attribute data and buffer objects
	struct VertexAttrib {
		int nelem;					// number of elements per attribute range: [1, 4]
		std::vector<float> data;
		unsigned int vbo;
		mutable bool vbo_valid;		// if this is false, the vbo needs updating from the data
		mutable bool data_valid;	// if this is false, the data needs to be pulled from the vbo
		//int sdr_loc;
	} vattr[NUM_MESH_ATTR];

	static int global_sdr_loc[NUM_MESH_ATTR];

	//std::vector<XFormNode*> bones;	// bones affecting this mesh

	// index data and buffer object
	std::vector<unsigned int> idata;
	unsigned int ibo;
	mutable bool ibo_valid;
	mutable bool idata_valid;

	// index buffer object for wireframe rendering (constructed on demand)
	unsigned int wire_ibo;
	mutable bool wire_ibo_valid;

	// axis-aligned bounding box
	mutable AABox aabb;
	mutable bool aabb_valid;

	// bounding sphere
	mutable Sphere bsph;
	mutable bool bsph_valid;

	// keeps the last intersected face
	mutable Triangle hitface;
	// keeps the last intersected vertex position
	mutable Vec3 hitvert;

	void calc_aabb();
	void calc_bsph();

	static unsigned int intersect_mode;
	static float vertex_sel_dist;

	static float vis_vecsize;

	/// update the VBOs after data has changed (invalid vbo/ibo)
	void update_buffers();
	/// construct/update the wireframe index buffer (called from draw_wire).
	void update_wire_ibo();

	mutable int cur_sdr;
	bool pre_draw() const;
	void post_draw() const;


public:
	static bool use_custom_sdr_attr;

	Mesh();
	~Mesh();

	Mesh(const Mesh &rhs);
	Mesh &operator =(const Mesh &rhs);
	bool clone(const Mesh &m);

	void set_name(const char *name);
	const char *get_name() const;

	bool has_attrib(int attr) const;
	bool is_indexed() const;

	// clears everything about this mesh, and returns to the newly constructed state
	void clear();

	// access the vertex attribute data
	// if vdata == 0, space is just allocated
	float *set_attrib_data(int attrib, int nelem, unsigned int num, const float *vdata = 0); // invalidates vbo
	float *get_attrib_data(int attrib);	// invalidates vbo
	const float *get_attrib_data(int attrib) const;

	// simple access to any particular attribute
	void set_attrib(int attrib, int idx, const Vec4 &v); // invalidates vbo
	Vec4 get_attrib(int attrib, int idx) const;

	int get_attrib_count(int attrib) const;

	// ... same for index data
	unsigned int *set_index_data(int num, const unsigned int *indices = 0); // invalidates ibo
	unsigned int *get_index_data();	// invalidates ibo
	const unsigned int *get_index_data() const;

	int get_index_count() const;

	void append(const Mesh &mesh);

	// immediate-mode style mesh construction interface
	void vertex(float x, float y, float z);
	void normal(float nx, float ny, float nz);
	void tangent(float tx, float ty, float tz);
	void texcoord(float u, float v, float w);
	void boneweights(float w1, float w2, float w3, float w4);
	void boneidx(int idx1, int idx2, int idx3, int idx4);

	int get_poly_count() const;

	/* apply a transformation to the vertices and its inverse-transpose
	 * to the normals and tangents.
	 */
	void apply_xform(const Mat4 &xform);
	void apply_xform(const Mat4 &xform, const Mat4 &dir_xform);

	void flip();	// both faces and normals
	void flip_faces();
	void flip_normals();

	void explode();	// undo all vertex sharing

	void calc_face_normals(); // this is only guaranteed to work on an exploded mesh

	// adds a bone and returns its index
	/*int add_bone(XFormNode *bone);
	const XFormNode *get_bone(int idx) const;
	int get_bones_count() const;*/

	// access the shader attribute locations
	static void set_attrib_location(int attr, int loc);
	static int get_attrib_location(int attr);
	static void clear_attrib_locations();

	static void set_vis_vecsize(float sz);
	static float get_vis_vecsize();

	void draw() const;
	void draw_wire() const;
	void draw_vertices() const;
	void draw_normals() const;
	void draw_tangents() const;

	/** get the bounding box in local space. The result will be cached, and subsequent
	 * calls will return the same box. The cache gets invalidated by any functions that can affect
	 * the vertex data (non-const variant of get_attrib_data(MESH_ATTR_VERTEX, ...) included).
	 * @{ */
	void get_aabbox(Vec3 *vmin, Vec3 *vmax) const;
	const AABox &get_aabbox() const;
	/// @}

	/** get the bounding sphere in local space. The result will be cached, and subsequent
	 * calls will return the same box. The cache gets invalidated by any functions that can affect
	 * the vertex data (non-const variant of get_attrib_data(MESH_ATTR_VERTEX, ...) included).
	 * @{ */
	float get_bsphere(Vec3 *center, float *rad) const;
	const Sphere &get_bsphere() const;

	static void set_intersect_mode(unsigned int mode);
	static unsigned int get_intersect_mode();
	static void set_vertex_select_distance(float dist);
	static float get_vertex_select_distance();

	/** Find the intersection between the mesh and a ray.
	 * XXX Brute force at the moment, not intended to be used for anything other than picking in tools.
	 *     If you intend to use it in a speed-critical part of the code, you'll *have* to optimize it!
	 */
	bool intersect(const Ray &ray, HitPoint *hit = 0) const;

	// texture coordinate manipulation
	void texcoord_apply_xform(const Mat4 &xform);
	void texcoord_gen_plane(const Vec3 &norm, const Vec3 &tang);
	void texcoord_gen_box();
	void texcoord_gen_cylinder();

	bool dump(const char *fname) const;
	bool dump(FILE *fp) const;
	bool dump_obj(const char *fname) const;
	bool dump_obj(FILE *fp) const;
};

#endif	// MESH_H_
