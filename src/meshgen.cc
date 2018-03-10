#include <stdio.h>
#include "meshgen.h"
#include "mesh.h"

// -------- sphere --------

#define SURAD(u)	((u) * 2.0 * M_PI)
#define SVRAD(v)	((v) * M_PI)

static Vec3 sphvec(float theta, float phi)
{
	return Vec3(sin(theta) * sin(phi),
			cos(phi),
			cos(theta) * sin(phi));
}

void gen_sphere(Mesh *mesh, float rad, int usub, int vsub, float urange, float vrange)
{
	if(urange == 0.0 || vrange == 0.0) return;

	if(usub < 4) usub = 4;
	if(vsub < 2) vsub = 2;

	int uverts = usub + 1;
	int vverts = vsub + 1;

	int num_verts = uverts * vverts;
	int num_quads = usub * vsub;
	int num_tri = num_quads * 2;

	mesh->clear();
	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);
	unsigned int *idxarr = mesh->set_index_data(num_tri * 3, 0);

	float du = urange / (float)(uverts - 1);
	float dv = vrange / (float)(vverts - 1);

	float u = 0.0;
	for(int i=0; i<uverts; i++) {
		float theta = u * 2.0 * M_PI;

		float v = 0.0;
		for(int j=0; j<vverts; j++) {
			float phi = v * M_PI;

			Vec3 pos = sphvec(theta, phi);

			*varr++ = pos * rad;
			*narr++ = pos;
			*tarr++ = normalize(sphvec(theta + 0.1f, (float)M_PI / 2.0f) - sphvec(theta - 0.1f, (float)M_PI / 2.0f));
			*uvarr++ = Vec2(u / urange, v / vrange);

			if(i < usub && j < vsub) {
				int idx = i * vverts + j;
				*idxarr++ = idx;
				*idxarr++ = idx + 1;
				*idxarr++ = idx + vverts + 1;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts + 1;
				*idxarr++ = idx + vverts;
			}

			v += dv;
		}
		u += du;
	}
}

// ------ geosphere ------
#define PHI		1.618034

static Vec3 icosa_pt[] = {
	Vec3(PHI, 1, 0),
	Vec3(-PHI, 1, 0),
	Vec3(PHI, -1, 0),
	Vec3(-PHI, -1, 0),
	Vec3(1, 0, PHI),
	Vec3(1, 0, -PHI),
	Vec3(-1, 0, PHI),
	Vec3(-1, 0, -PHI),
	Vec3(0, PHI, 1),
	Vec3(0, -PHI, 1),
	Vec3(0, PHI, -1),
	Vec3(0, -PHI, -1)
};
enum { P11, P12, P13, P14, P21, P22, P23, P24, P31, P32, P33, P34 };
static int icosa_idx[] = {
	P11, P31, P21,
	P11, P22, P33,
	P13, P21, P32,
	P13, P34, P22,
	P12, P23, P31,
	P12, P33, P24,
	P14, P32, P23,
	P14, P24, P34,

	P11, P33, P31,
	P12, P31, P33,
	P13, P32, P34,
	P14, P34, P32,

	P21, P13, P11,
	P22, P11, P13,
	P23, P12, P14,
	P24, P14, P12,

	P31, P23, P21,
	P32, P21, P23,
	P33, P22, P24,
	P34, P24, P22
};

static void geosphere(std::vector<Vec3> *verts, const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, int iter)
{
	if(!iter) {
		verts->push_back(v1);
		verts->push_back(v2);
		verts->push_back(v3);
		return;
	}

	Vec3 v12 = normalize(v1 + v2);
	Vec3 v23 = normalize(v2 + v3);
	Vec3 v31 = normalize(v3 + v1);

	geosphere(verts, v1, v12, v31, iter - 1);
	geosphere(verts, v2, v23, v12, iter - 1);
	geosphere(verts, v3, v31, v23, iter - 1);
	geosphere(verts, v12, v23, v31, iter - 1);
}

void gen_geosphere(Mesh *mesh, float rad, int subdiv, bool hemi)
{
	int num_tri = (sizeof icosa_idx / sizeof *icosa_idx) / 3;

	std::vector<Vec3> verts;
	for(int i=0; i<num_tri; i++) {
		Vec3 v[3];

		for(int j=0; j<3; j++) {
			int vidx = icosa_idx[i * 3 + j];
			v[j] = normalize(icosa_pt[vidx]);
		}

		if(hemi && (v[0].y < 0.0 || v[1].y < 0.0 || v[2].y < 0.0)) {
			continue;
		}

		geosphere(&verts, v[0], v[1], v[2], subdiv);
	}

	int num_verts = (int)verts.size();

	mesh->clear();
	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);

	for(int i=0; i<num_verts; i++) {
		*varr++ = verts[i] * rad;
		*narr++ = verts[i];

		float theta = atan2(verts[i].z, verts[i].x);
		float phi = acos(verts[i].y);

		*tarr++ = normalize(sphvec(theta + 0.1f, (float)M_PI / 2.0f) - sphvec(theta - 0.1f, (float)M_PI / 2.0f));

		float u = 0.5 * theta / M_PI + 0.5;
		float v = phi / M_PI;
		*uvarr++ = Vec2(u, v);
	}
}

// -------- torus -----------
static Vec3 torusvec(float theta, float phi, float mr, float rr)
{
	theta = -theta;

	float rx = -cos(phi) * rr + mr;
	float ry = sin(phi) * rr;
	float rz = 0.0;

	float x = rx * sin(theta) + rz * cos(theta);
	float y = ry;
	float z = -rx * cos(theta) + rz * sin(theta);

	return Vec3(x, y, z);
}

void gen_torus(Mesh *mesh, float mainrad, float ringrad, int usub, int vsub, float urange, float vrange)
{
	if(usub < 4) usub = 4;
	if(vsub < 2) vsub = 2;

	int uverts = usub + 1;
	int vverts = vsub + 1;

	int num_verts = uverts * vverts;
	int num_quads = usub * vsub;
	int num_tri = num_quads * 2;

	mesh->clear();
	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);
	unsigned int *idxarr = mesh->set_index_data(num_tri * 3, 0);

	float du = urange / (float)(uverts - 1);
	float dv = vrange / (float)(vverts - 1);

	float u = 0.0;
	for(int i=0; i<uverts; i++) {
		float theta = u * 2.0 * M_PI;

		float v = 0.0;
		for(int j=0; j<vverts; j++) {
			float phi = v * 2.0 * M_PI;

			Vec3 pos = torusvec(theta, phi, mainrad, ringrad);
			Vec3 cent = torusvec(theta, phi, mainrad, 0.0);

			*varr++ = pos;
			*narr++ = (pos - cent) / ringrad;

			Vec3 pprev = torusvec(theta - 0.1f, phi, mainrad, ringrad);
			Vec3 pnext = torusvec(theta + 0.1f, phi, mainrad, ringrad);

			*tarr++ = normalize(pnext - pprev);
			*uvarr++ = Vec2(u * urange, v * vrange);

			if(i < usub && j < vsub) {
				int idx = i * vverts + j;
				*idxarr++ = idx;
				*idxarr++ = idx + 1;
				*idxarr++ = idx + vverts + 1;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts + 1;
				*idxarr++ = idx + vverts;
			}

			v += dv;
		}
		u += du;
	}
}


// -------- cylinder --------

static Vec3 cylvec(float theta, float height)
{
	return Vec3(sin(theta), height, cos(theta));
}

void gen_cylinder(Mesh *mesh, float rad, float height, int usub, int vsub, int capsub, float urange, float vrange)
{
	if(usub < 4) usub = 4;
	if(vsub < 1) vsub = 1;

	int uverts = usub + 1;
	int vverts = vsub + 1;

	int num_body_verts = uverts * vverts;
	int num_body_quads = usub * vsub;
	int num_body_tri = num_body_quads * 2;

	int capvverts = capsub ? capsub + 1 : 0;
	int num_cap_verts = uverts * capvverts;
	int num_cap_quads = usub * capsub;
	int num_cap_tri = num_cap_quads * 2;

	int num_verts = num_body_verts + num_cap_verts * 2;
	int num_tri = num_body_tri + num_cap_tri * 2;

	mesh->clear();
	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);
	unsigned int *idxarr = mesh->set_index_data(num_tri * 3, 0);

	float du = urange / (float)(uverts - 1);
	float dv = vrange / (float)(vverts - 1);

	float u = 0.0;
	for(int i=0; i<uverts; i++) {
		float theta = SURAD(u);

		float v = 0.0;
		for(int j=0; j<vverts; j++) {
			float y = (v - 0.5) * height;
			Vec3 pos = cylvec(theta, y);

			*varr++ = Vec3(pos.x * rad, pos.y, pos.z * rad);
			*narr++ = Vec3(pos.x, 0.0, pos.z);
			*tarr++ = normalize(cylvec(theta + 0.1, 0.0) - cylvec(theta - 0.1, 0.0));
			*uvarr++ = Vec2(u * urange, v * vrange);

			if(i < usub && j < vsub) {
				int idx = i * vverts + j;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts + 1;
				*idxarr++ = idx + 1;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts;
				*idxarr++ = idx + vverts + 1;
			}

			v += dv;
		}
		u += du;
	}


	// now the cap!
	if(!capsub) {
		return;
	}

	dv = 1.0 / (float)(capvverts - 1);

	u = 0.0;
	for(int i=0; i<uverts; i++) {
		float theta = SURAD(u);

		float v = 0.0;
		for(int j=0; j<capvverts; j++) {
			float r = v * rad;

			Vec3 pos = cylvec(theta, height / 2.0) * r;
			pos.y = height / 2.0;
			Vec3 tang = normalize(cylvec(theta + 0.1, 0.0) - cylvec(theta - 0.1, 0.0));

			*varr++ = pos;
			*narr++ = Vec3(0, 1, 0);
			*tarr++ = tang;
			*uvarr++ = Vec2(u * urange, v);

			pos.y = -height / 2.0;
			*varr++ = pos;
			*narr++ = Vec3(0, -1, 0);
			*tarr++ = -tang;
			*uvarr++ = Vec2(u * urange, v);

			if(i < usub && j < capsub) {
				unsigned int idx = num_body_verts + (i * capvverts + j) * 2;

				unsigned int vidx[4] = {
					idx,
					idx + capvverts * 2,
					idx + (capvverts + 1) * 2,
					idx + 2
				};

				*idxarr++ = vidx[0];
				*idxarr++ = vidx[2];
				*idxarr++ = vidx[1];
				*idxarr++ = vidx[0];
				*idxarr++ = vidx[3];
				*idxarr++ = vidx[2];

				*idxarr++ = vidx[0] + 1;
				*idxarr++ = vidx[1] + 1;
				*idxarr++ = vidx[2] + 1;
				*idxarr++ = vidx[0] + 1;
				*idxarr++ = vidx[2] + 1;
				*idxarr++ = vidx[3] + 1;
			}

			v += dv;
		}
		u += du;
	}
}

// -------- cone --------

static Vec3 conevec(float theta, float y, float height)
{
	float scale = 1.0 - y / height;
	return Vec3(sin(theta) * scale, y, cos(theta) * scale);
}

void gen_cone(Mesh *mesh, float rad, float height, int usub, int vsub, int capsub, float urange, float vrange)
{
	if(usub < 4) usub = 4;
	if(vsub < 1) vsub = 1;

	int uverts = usub + 1;
	int vverts = vsub + 1;

	int num_body_verts = uverts * vverts;
	int num_body_quads = usub * vsub;
	int num_body_tri = num_body_quads * 2;

	int capvverts = capsub ? capsub + 1 : 0;
	int num_cap_verts = uverts * capvverts;
	int num_cap_quads = usub * capsub;
	int num_cap_tri = num_cap_quads * 2;

	int num_verts = num_body_verts + num_cap_verts;
	int num_tri = num_body_tri + num_cap_tri;

	mesh->clear();
	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);
	unsigned int *idxarr = mesh->set_index_data(num_tri * 3, 0);

	float du = urange / (float)(uverts - 1);
	float dv = vrange / (float)(vverts - 1);

	float u = 0.0;
	for(int i=0; i<uverts; i++) {
		float theta = SURAD(u);

		float v = 0.0;
		for(int j=0; j<vverts; j++) {
			float y = v * height;
			Vec3 pos = conevec(theta, y, height);

			Vec3 tang = normalize(conevec(theta + 0.1, 0.0, height) - conevec(theta - 0.1, 0.0, height));
			Vec3 bitang = normalize(conevec(theta, y + 0.1, height) - pos);

			*varr++ = Vec3(pos.x * rad, pos.y, pos.z * rad);
			*narr++ = cross(tang, bitang);
			*tarr++ = tang;
			*uvarr++ = Vec2(u * urange, v * vrange);

			if(i < usub && j < vsub) {
				int idx = i * vverts + j;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts + 1;
				*idxarr++ = idx + 1;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts;
				*idxarr++ = idx + vverts + 1;
			}

			v += dv;
		}
		u += du;
	}


	// now the bottom cap!
	if(!capsub) {
		return;
	}

	dv = 1.0 / (float)(capvverts - 1);

	u = 0.0;
	for(int i=0; i<uverts; i++) {
		float theta = SURAD(u);

		float v = 0.0;
		for(int j=0; j<capvverts; j++) {
			float r = v * rad;

			Vec3 pos = conevec(theta, 0.0, height) * r;
			Vec3 tang = normalize(cylvec(theta + 0.1, 0.0) - cylvec(theta - 0.1, 0.0));

			*varr++ = pos;
			*narr++ = Vec3(0, -1, 0);
			*tarr++ = tang;
			*uvarr++ = Vec2(u * urange, v);

			if(i < usub && j < capsub) {
				unsigned int idx = num_body_verts + i * capvverts + j;

				unsigned int vidx[4] = {
					idx,
					idx + capvverts,
					idx + (capvverts + 1),
					idx + 1
				};

				*idxarr++ = vidx[0];
				*idxarr++ = vidx[1];
				*idxarr++ = vidx[2];
				*idxarr++ = vidx[0];
				*idxarr++ = vidx[2];
				*idxarr++ = vidx[3];
			}

			v += dv;
		}
		u += du;
	}
}


// -------- plane --------

void gen_plane(Mesh *mesh, float width, float height, int usub, int vsub)
{
	gen_heightmap(mesh, width, height, usub, vsub, 0);
}


// ----- heightmap ------

void gen_heightmap(Mesh *mesh, float width, float height, int usub, int vsub, float (*hf)(float, float, void*), void *hfdata)
{
	if(usub < 1) usub = 1;
	if(vsub < 1) vsub = 1;

	mesh->clear();

	int uverts = usub + 1;
	int vverts = vsub + 1;
	int num_verts = uverts * vverts;

	int num_quads = usub * vsub;
	int num_tri = num_quads * 2;

	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);
	unsigned int *idxarr = mesh->set_index_data(num_tri * 3, 0);

	float du = 1.0 / (float)usub;
	float dv = 1.0 / (float)vsub;

	float u = 0.0;
	for(int i=0; i<uverts; i++) {
		float v = 0.0;
		for(int j=0; j<vverts; j++) {
			float x = (u - 0.5) * width;
			float y = (v - 0.5) * height;
			float z = hf ? hf(u, v, hfdata) : 0.0;

			Vec3 normal = Vec3(0, 0, 1);
			if(hf) {
				float u1z = hf(u + du, v, hfdata);
				float v1z = hf(u, v + dv, hfdata);

				Vec3 tang = Vec3(du * width, 0, u1z - z);
				Vec3 bitan = Vec3(0, dv * height, v1z - z);
				normal = normalize(cross(tang, bitan));
			}

			*varr++ = Vec3(x, y, z);
			*narr++ = normal;
			*tarr++ = Vec3(1, 0, 0);
			*uvarr++ = Vec2(u, v);

			if(i < usub && j < vsub) {
				int idx = i * vverts + j;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts + 1;
				*idxarr++ = idx + 1;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts;
				*idxarr++ = idx + vverts + 1;
			}

			v += dv;
		}
		u += du;
	}
}

// ----- box ------
void gen_box(Mesh *mesh, float xsz, float ysz, float zsz, int usub, int vsub)
{
	static const float face_angles[][2] = {
		{0, 0},
		{M_PI / 2.0, 0},
		{M_PI, 0},
		{3.0 * M_PI / 2.0, 0},
		{0, M_PI / 2.0},
		{0, -M_PI / 2.0}
	};

	if(usub < 1) usub = 1;
	if(vsub < 1) vsub = 1;

	mesh->clear();

	for(int i=0; i<6; i++) {
		Mat4 xform, dir_xform;
		Mesh m;

		gen_plane(&m, 1, 1, usub, vsub);
		xform.translate(Vec3(0, 0, 0.5));
		xform.rotate(Vec3(face_angles[i][1], face_angles[i][0], 0));
		dir_xform = xform;
		m.apply_xform(xform, dir_xform);

		mesh->append(m);
	}

	Mat4 scale;
	scale.scaling(xsz, ysz, zsz);
	mesh->apply_xform(scale, Mat4::identity);
}

/*
void gen_box(Mesh *mesh, float xsz, float ysz, float zsz)
{
	mesh->clear();

	const int num_faces = 6;
	int num_verts = num_faces * 4;
	int num_tri = num_faces * 2;

	float x = xsz / 2.0;
	float y = ysz / 2.0;
	float z = zsz / 2.0;

	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);
	unsigned int *idxarr = mesh->set_index_data(num_tri * 3, 0);

	static const Vec2 uv[] = { Vec2(0, 0), Vec2(1, 0), Vec2(1, 1), Vec2(0, 1) };

	// front
	for(int i=0; i<4; i++) {
		*narr++ = Vec3(0, 0, 1);
		*tarr++ = Vec3(1, 0, 0);
		*uvarr++ = uv[i];
	}
	*varr++ = Vec3(-x, -y, z);
	*varr++ = Vec3(x, -y, z);
	*varr++ = Vec3(x, y, z);
	*varr++ = Vec3(-x, y, z);
	// right
	for(int i=0; i<4; i++) {
		*narr++ = Vec3(1, 0, 0);
		*tarr++ = Vec3(0, 0, -1);
		*uvarr++ = uv[i];
	}
	*varr++ = Vec3(x, -y, z);
	*varr++ = Vec3(x, -y, -z);
	*varr++ = Vec3(x, y, -z);
	*varr++ = Vec3(x, y, z);
	// back
	for(int i=0; i<4; i++) {
		*narr++ = Vec3(0, 0, -1);
		*tarr++ = Vec3(-1, 0, 0);
		*uvarr++ = uv[i];
	}
	*varr++ = Vec3(x, -y, -z);
	*varr++ = Vec3(-x, -y, -z);
	*varr++ = Vec3(-x, y, -z);
	*varr++ = Vec3(x, y, -z);
	// left
	for(int i=0; i<4; i++) {
		*narr++ = Vec3(-1, 0, 0);
		*tarr++ = Vec3(0, 0, 1);
		*uvarr++ = uv[i];
	}
	*varr++ = Vec3(-x, -y, -z);
	*varr++ = Vec3(-x, -y, z);
	*varr++ = Vec3(-x, y, z);
	*varr++ = Vec3(-x, y, -z);
	// top
	for(int i=0; i<4; i++) {
		*narr++ = Vec3(0, 1, 0);
		*tarr++ = Vec3(1, 0, 0);
		*uvarr++ = uv[i];
	}
	*varr++ = Vec3(-x, y, z);
	*varr++ = Vec3(x, y, z);
	*varr++ = Vec3(x, y, -z);
	*varr++ = Vec3(-x, y, -z);
	// bottom
	for(int i=0; i<4; i++) {
		*narr++ = Vec3(0, -1, 0);
		*tarr++ = Vec3(1, 0, 0);
		*uvarr++ = uv[i];
	}
	*varr++ = Vec3(-x, -y, -z);
	*varr++ = Vec3(x, -y, -z);
	*varr++ = Vec3(x, -y, z);
	*varr++ = Vec3(-x, -y, z);

	// index array
	static const int faceidx[] = {0, 1, 2, 0, 2, 3};
	for(int i=0; i<num_faces; i++) {
		for(int j=0; j<6; j++) {
			*idxarr++ = faceidx[j] + i * 4;
		}
	}
}
*/

static inline Vec3 rev_vert(float u, float v, Vec2 (*rf)(float, float, void*), void *cls)
{
	Vec2 pos = rf(u, v, cls);

	float angle = u * 2.0 * M_PI;
	float x = pos.x * cos(angle);
	float y = pos.y;
	float z = pos.x * sin(angle);

	return Vec3(x, y, z);
}

// ------ surface of revolution -------
void gen_revol(Mesh *mesh, int usub, int vsub, Vec2 (*rfunc)(float, float, void*), void *cls)
{
	gen_revol(mesh, usub, vsub, rfunc, 0, cls);
}

void gen_revol(Mesh *mesh, int usub, int vsub, Vec2 (*rfunc)(float, float, void*),
		Vec2 (*nfunc)(float, float, void*), void *cls)
{
	if(!rfunc) return;
	if(usub < 3) usub = 3;
	if(vsub < 1) vsub = 1;

	mesh->clear();

	int uverts = usub + 1;
	int vverts = vsub + 1;
	int num_verts = uverts * vverts;

	int num_quads = usub * vsub;
	int num_tri = num_quads * 2;

	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);
	unsigned int *idxarr = mesh->set_index_data(num_tri * 3, 0);

	float du = 1.0 / (float)(uverts - 1);
	float dv = 1.0 / (float)(vverts - 1);

	float u = 0.0;
	for(int i=0; i<uverts; i++) {
		float v = 0.0;
		for(int j=0; j<vverts; j++) {
			Vec3 pos = rev_vert(u, v, rfunc, cls);

			Vec3 nextu = rev_vert(fmod(u + du, 1.0), v, rfunc, cls);
			Vec3 tang = nextu - pos;
			if(length_sq(tang) < 1e-6) {
				float new_v = v > 0.5 ? v - dv * 0.25 : v + dv * 0.25;
				nextu = rev_vert(fmod(u + du, 1.0), new_v, rfunc, cls);
				tang = nextu - pos;
			}

			Vec3 normal;
			if(nfunc) {
				normal = rev_vert(u, v, nfunc, cls);
			} else {
				Vec3 nextv = rev_vert(u, v + dv, rfunc, cls);
				Vec3 bitan = nextv - pos;
				if(length_sq(bitan) < 1e-6) {
					nextv = rev_vert(u, v - dv, rfunc, cls);
					bitan = pos - nextv;
				}

				normal = cross(tang, bitan);
			}

			*varr++ = pos;
			*narr++ = normalize(normal);
			*tarr++ = normalize(tang);
			*uvarr++ = Vec2(u, v);

			if(i < usub && j < vsub) {
				int idx = i * vverts + j;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts + 1;
				*idxarr++ = idx + 1;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts;
				*idxarr++ = idx + vverts + 1;
			}

			v += dv;
		}
		u += du;
	}
}


static inline Vec3 sweep_vert(float u, float v, float height, Vec2 (*sf)(float, float, void*), void *cls)
{
	Vec2 pos = sf(u, v, cls);

	float x = pos.x;
	float y = v * height;
	float z = pos.y;

	return Vec3(x, y, z);
}

// ---- sweep shape along a path ----
void gen_sweep(Mesh *mesh, float height, int usub, int vsub, Vec2 (*sfunc)(float, float, void*), void *cls)
{
	if(!sfunc) return;
	if(usub < 3) usub = 3;
	if(vsub < 1) vsub = 1;

	mesh->clear();

	int uverts = usub + 1;
	int vverts = vsub + 1;
	int num_verts = uverts * vverts;

	int num_quads = usub * vsub;
	int num_tri = num_quads * 2;

	Vec3 *varr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_VERTEX, 3, num_verts, 0);
	Vec3 *narr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_NORMAL, 3, num_verts, 0);
	Vec3 *tarr = (Vec3*)mesh->set_attrib_data(MESH_ATTR_TANGENT, 3, num_verts, 0);
	Vec2 *uvarr = (Vec2*)mesh->set_attrib_data(MESH_ATTR_TEXCOORD, 2, num_verts, 0);
	unsigned int *idxarr = mesh->set_index_data(num_tri * 3, 0);

	float du = 1.0 / (float)(uverts - 1);
	float dv = 1.0 / (float)(vverts - 1);

	float u = 0.0;
	for(int i=0; i<uverts; i++) {
		float v = 0.0;
		for(int j=0; j<vverts; j++) {
			Vec3 pos = sweep_vert(u, v, height, sfunc, cls);

			Vec3 nextu = sweep_vert(fmod(u + du, 1.0), v, height, sfunc, cls);
			Vec3 tang = nextu - pos;
			if(length_sq(tang) < 1e-6) {
				float new_v = v > 0.5 ? v - dv * 0.25 : v + dv * 0.25;
				nextu = sweep_vert(fmod(u + du, 1.0), new_v, height, sfunc, cls);
				tang = nextu - pos;
			}

			Vec3 normal;
			Vec3 nextv = sweep_vert(u, v + dv, height, sfunc, cls);
			Vec3 bitan = nextv - pos;
			if(length_sq(bitan) < 1e-6) {
				nextv = sweep_vert(u, v - dv, height, sfunc, cls);
				bitan = pos - nextv;
			}

			normal = cross(tang, bitan);

			*varr++ = pos;
			*narr++ = normalize(normal);
			*tarr++ = normalize(tang);
			*uvarr++ = Vec2(u, v);

			if(i < usub && j < vsub) {
				int idx = i * vverts + j;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts + 1;
				*idxarr++ = idx + 1;

				*idxarr++ = idx;
				*idxarr++ = idx + vverts;
				*idxarr++ = idx + vverts + 1;
			}

			v += dv;
		}
		u += du;
	}
}
