#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include "opengl.h"
#include "mesh.h"
//#include "logger.h"
//#include "xform_node.h"

#define USE_OLDGL

bool Mesh::use_custom_sdr_attr = true;
int Mesh::global_sdr_loc[NUM_MESH_ATTR] = { 0, 1, 2, 3, 4, 5, 6, 7 };
unsigned int Mesh::intersect_mode = ISECT_DEFAULT;
float Mesh::vertex_sel_dist = 0.01;
float Mesh::vis_vecsize = 1.0;

Mesh::Mesh()
{
	clear();

	glGenBuffers(NUM_MESH_ATTR + 1, buffer_objects);

	for(int i=0; i<NUM_MESH_ATTR; i++) {
		vattr[i].vbo = buffer_objects[i];
	}
	ibo = buffer_objects[NUM_MESH_ATTR];
	wire_ibo = 0;
}

Mesh::~Mesh()
{
	glDeleteBuffers(NUM_MESH_ATTR + 1, buffer_objects);

	if(wire_ibo) {
		glDeleteBuffers(1, &wire_ibo);
	}
}

Mesh::Mesh(const Mesh &rhs)
{
	clear();

	glGenBuffers(NUM_MESH_ATTR + 1, buffer_objects);

	for(int i=0; i<NUM_MESH_ATTR; i++) {
		vattr[i].vbo = buffer_objects[i];
	}
	ibo = buffer_objects[NUM_MESH_ATTR];
	wire_ibo = 0;

	clone(rhs);
}

Mesh &Mesh::operator =(const Mesh &rhs)
{
	if(&rhs != this) {
		clone(rhs);
	}
	return *this;
}

bool Mesh::clone(const Mesh &m)
{
	clear();

	for(int i=0; i<NUM_MESH_ATTR; i++) {
		if(m.has_attrib(i)) {
			m.get_attrib_data(i);	// force validation of the actual data on the source mesh

			vattr[i].nelem = m.vattr[i].nelem;
			vattr[i].data = m.vattr[i].data;	// copy the actual data
			vattr[i].data_valid = true;
		}
	}

	if(m.is_indexed()) {
		m.get_index_data();		// again, force validation

		// copy the index data
		idata = m.idata;
		idata_valid = true;
	}

	name = m.name;
	nverts = m.nverts;
	nfaces = m.nfaces;

	//bones = m.bones;

	memcpy(cur_val, m.cur_val, sizeof cur_val);

	aabb = m.aabb;
	aabb_valid = m.aabb_valid;
	bsph = m.bsph;
	bsph_valid = m.bsph_valid;

	hitface = m.hitface;
	hitvert = m.hitvert;

	intersect_mode = m.intersect_mode;
	vertex_sel_dist = m.vertex_sel_dist;
	vis_vecsize = m.vis_vecsize;

	return true;
}

void Mesh::set_name(const char *name)
{
	this->name = name;
}

const char *Mesh::get_name() const
{
	return name.c_str();
}

bool Mesh::has_attrib(int attr) const
{
	if(attr < 0 || attr >= NUM_MESH_ATTR) {
		return false;
	}

	// if neither of these is valid, then nobody has set this attribute
	return vattr[attr].vbo_valid || vattr[attr].data_valid;
}

bool Mesh::is_indexed() const
{
	return ibo_valid || idata_valid;
}

void Mesh::clear()
{
	//bones.clear();

	for(int i=0; i<NUM_MESH_ATTR; i++) {
		vattr[i].nelem = 0;
		vattr[i].vbo_valid = false;
		vattr[i].data_valid = false;
		//vattr[i].sdr_loc = -1;
		vattr[i].data.clear();
	}
	ibo_valid = idata_valid = false;
	idata.clear();

	wire_ibo_valid = false;

	nverts = nfaces = 0;

	bsph_valid = false;
	aabb_valid = false;
}

float *Mesh::set_attrib_data(int attrib, int nelem, unsigned int num, const float *data)
{
	if(attrib < 0 || attrib >= NUM_MESH_ATTR) {
		fprintf(stderr, "%s: invalid attrib: %d\n", __FUNCTION__, attrib);
		return 0;
	}

	if(nverts && num != nverts) {
		fprintf(stderr, "%s: attribute count missmatch (%d instead of %d)\n", __FUNCTION__, num, nverts);
		return 0;
	}
	nverts = num;

	vattr[attrib].data.clear();
	vattr[attrib].nelem = nelem;
	vattr[attrib].data.resize(num * nelem);

	if(data) {
		memcpy(&vattr[attrib].data[0], data, num * nelem * sizeof *data);
	}

	vattr[attrib].data_valid = true;
	vattr[attrib].vbo_valid = false;
	return &vattr[attrib].data[0];
}

float *Mesh::get_attrib_data(int attrib)
{
	if(attrib < 0 || attrib >= NUM_MESH_ATTR) {
		fprintf(stderr, "%s: invalid attrib: %d\n", __FUNCTION__, attrib);
		return 0;
	}

	vattr[attrib].vbo_valid = false;
	return (float*)((const Mesh*)this)->get_attrib_data(attrib);
}

const float *Mesh::get_attrib_data(int attrib) const
{
	if(attrib < 0 || attrib >= NUM_MESH_ATTR) {
		fprintf(stderr, "%s: invalid attrib: %d\n", __FUNCTION__, attrib);
		return 0;
	}

	if(!vattr[attrib].data_valid) {
#if GL_ES_VERSION_2_0
		fprintf(stderr, "%s: can't read back attrib data on CrippledGL ES\n", __FUNCTION__);
		return 0;
#else
		if(!vattr[attrib].vbo_valid) {
			fprintf(stderr, "%s: unavailable attrib: %d\n", __FUNCTION__, attrib);
			return 0;
		}

		// local data copy is unavailable, grab the data from the vbo
		Mesh *m = (Mesh*)this;
		m->vattr[attrib].data.resize(nverts * vattr[attrib].nelem);

		glBindBuffer(GL_ARRAY_BUFFER, vattr[attrib].vbo);
		void *data = glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
		memcpy(&m->vattr[attrib].data[0], data, nverts * vattr[attrib].nelem * sizeof(float));
		glUnmapBuffer(GL_ARRAY_BUFFER);

		vattr[attrib].data_valid = true;
#endif
	}

	return &vattr[attrib].data[0];
}

void Mesh::set_attrib(int attrib, int idx, const Vec4 &v)
{
	float *data = get_attrib_data(attrib);
	if(data) {
		data += idx * vattr[attrib].nelem;
		for(int i=0; i<vattr[attrib].nelem; i++) {
			data[i] = v[i];
		}
	}
}

Vec4 Mesh::get_attrib(int attrib, int idx) const
{
	Vec4 v(0.0, 0.0, 0.0, 1.0);
	const float *data = get_attrib_data(attrib);
	if(data) {
		data += idx * vattr[attrib].nelem;
		for(int i=0; i<vattr[attrib].nelem; i++) {
			v[i] = data[i];
		}
	}
	return v;
}

int Mesh::get_attrib_count(int attrib) const
{
	return has_attrib(attrib) ? nverts : 0;
}


unsigned int *Mesh::set_index_data(int num, const unsigned int *indices)
{
	int nidx = nfaces * 3;
	if(nidx && num != nidx) {
		fprintf(stderr, "%s: index count mismatch (%d instead of %d)\n", __FUNCTION__, num, nidx);
		return 0;
	}
	nfaces = num / 3;

	idata.clear();
	idata.resize(num);

	if(indices) {
		memcpy(&idata[0], indices, num * sizeof *indices);
	}

	idata_valid = true;
	ibo_valid = false;

	return &idata[0];
}

unsigned int *Mesh::get_index_data()
{
	ibo_valid = false;
	return (unsigned int*)((const Mesh*)this)->get_index_data();
}

const unsigned int *Mesh::get_index_data() const
{
	if(!idata_valid) {
#if GL_ES_VERSION_2_0
		fprintf(stderr, "%s: can't read back index data in CrippledGL ES\n", __FUNCTION__);
		return 0;
#else
		if(!ibo_valid) {
			fprintf(stderr, "%s: indices unavailable\n", __FUNCTION__);
			return 0;
		}

		// local data copy is unavailable, gram the data from the ibo
		Mesh *m = (Mesh*)this;
		int nidx = nfaces * 3;
		m->idata.resize(nidx);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		void *data = glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_READ_ONLY);
		memcpy(&m->idata[0], data, nidx * sizeof(unsigned int));
		glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);

		idata_valid = true;
#endif
	}

	return &idata[0];
}

int Mesh::get_index_count() const
{
	return nfaces * 3;
}

void Mesh::append(const Mesh &mesh)
{
	unsigned int idxoffs = nverts;

	if(!nverts) {
		clone(mesh);
		return;
	}

	nverts += mesh.nverts;
	nfaces += mesh.nfaces;

	for(int i=0; i<NUM_MESH_ATTR; i++) {
		if(has_attrib(i) && mesh.has_attrib(i)) {
			// force validating the data arrays
			get_attrib_data(i);
			mesh.get_attrib_data(i);

			// append the mesh data
			vattr[i].data.insert(vattr[i].data.end(), mesh.vattr[i].data.begin(), mesh.vattr[i].data.end());
		}
	}

	if(ibo_valid || idata_valid) {
		// make index arrays valid
		get_index_data();
		mesh.get_index_data();

		size_t orig_sz = idata.size();

		idata.insert(idata.end(), mesh.idata.begin(), mesh.idata.end());

		// fixup all the new indices
		for(size_t i=orig_sz; i<idata.size(); i++) {
			idata[i] += idxoffs;
		}
	}

	// fuck everything
	wire_ibo_valid = false;
	aabb_valid = false;
	bsph_valid = false;
}

// assemble a complete vertex by adding all the useful attributes
void Mesh::vertex(float x, float y, float z)
{
	cur_val[MESH_ATTR_VERTEX] = Vec4(x, y, z, 1.0f);
	vattr[MESH_ATTR_VERTEX].data_valid = true;
	vattr[MESH_ATTR_VERTEX].nelem = 3;

	for(int i=0; i<NUM_MESH_ATTR; i++) {
		if(vattr[i].data_valid) {
			for(int j=0; j<vattr[MESH_ATTR_VERTEX].nelem; j++) {
				vattr[i].data.push_back(cur_val[i][j]);
			}
		}
		vattr[i].vbo_valid = false;
	}

	if(idata_valid) {
		idata.clear();
	}
	ibo_valid = idata_valid = false;
}

void Mesh::normal(float nx, float ny, float nz)
{
	cur_val[MESH_ATTR_NORMAL] = Vec4(nx, ny, nz, 1.0f);
	vattr[MESH_ATTR_NORMAL].data_valid = true;
	vattr[MESH_ATTR_NORMAL].nelem = 3;
}

void Mesh::tangent(float tx, float ty, float tz)
{
	cur_val[MESH_ATTR_TANGENT] = Vec4(tx, ty, tz, 1.0f);
	vattr[MESH_ATTR_TANGENT].data_valid = true;
	vattr[MESH_ATTR_TANGENT].nelem = 3;
}

void Mesh::texcoord(float u, float v, float w)
{
	cur_val[MESH_ATTR_TEXCOORD] = Vec4(u, v, w, 1.0f);
	vattr[MESH_ATTR_TEXCOORD].data_valid = true;
	vattr[MESH_ATTR_TEXCOORD].nelem = 3;
}

void Mesh::boneweights(float w1, float w2, float w3, float w4)
{
	cur_val[MESH_ATTR_BONEWEIGHTS] = Vec4(w1, w2, w3, w4);
	vattr[MESH_ATTR_BONEWEIGHTS].data_valid = true;
	vattr[MESH_ATTR_BONEWEIGHTS].nelem = 4;
}

void Mesh::boneidx(int idx1, int idx2, int idx3, int idx4)
{
	cur_val[MESH_ATTR_BONEIDX] = Vec4(idx1, idx2, idx3, idx4);
	vattr[MESH_ATTR_BONEIDX].data_valid = true;
	vattr[MESH_ATTR_BONEIDX].nelem = 4;
}

int Mesh::get_poly_count() const
{
	if(nfaces) {
		return nfaces;
	}
	if(nverts) {
		return nverts / 3;
	}
	return 0;
}

/// static function
void Mesh::set_attrib_location(int attr, int loc)
{
	if(attr < 0 || attr >= NUM_MESH_ATTR) {
		return;
	}
	Mesh::global_sdr_loc[attr] = loc;
}

/// static function
int Mesh::get_attrib_location(int attr)
{
	if(attr < 0 || attr >= NUM_MESH_ATTR) {
		return -1;
	}
	return Mesh::global_sdr_loc[attr];
}

/// static function
void Mesh::clear_attrib_locations()
{
	for(int i=0; i<NUM_MESH_ATTR; i++) {
		Mesh::global_sdr_loc[i] = -1;
	}
}

/// static function
void Mesh::set_vis_vecsize(float sz)
{
	Mesh::vis_vecsize = sz;
}

float Mesh::get_vis_vecsize()
{
	return Mesh::vis_vecsize;
}

void Mesh::apply_xform(const Mat4 &xform)
{
	Mat4 dir_xform = xform.upper3x3();
	apply_xform(xform, dir_xform);
}

void Mesh::apply_xform(const Mat4 &xform, const Mat4 &dir_xform)
{
	for(unsigned int i=0; i<nverts; i++) {
		Vec4 v = get_attrib(MESH_ATTR_VERTEX, i);
		set_attrib(MESH_ATTR_VERTEX, i, xform * v);

		if(has_attrib(MESH_ATTR_NORMAL)) {
			Vec3 n = Vec3(get_attrib(MESH_ATTR_NORMAL, i));
			set_attrib(MESH_ATTR_NORMAL, i, Vec4(dir_xform * n));
		}
		if(has_attrib(MESH_ATTR_TANGENT)) {
			Vec3 t = Vec3(get_attrib(MESH_ATTR_TANGENT, i));
			set_attrib(MESH_ATTR_TANGENT, i, Vec4(dir_xform * t));
		}
	}
}

void Mesh::flip()
{
	flip_faces();
	flip_normals();
}

void Mesh::flip_faces()
{
	if(is_indexed()) {
		unsigned int *indices = get_index_data();
		if(!indices) return;

		int idxnum = get_index_count();
		for(int i=0; i<idxnum; i+=3) {
			unsigned int tmp = indices[i + 2];
			indices[i + 2] = indices[i + 1];
			indices[i + 1] = tmp;
		}

	} else {
		Vec3 *verts = (Vec3*)get_attrib_data(MESH_ATTR_VERTEX);
		if(!verts) return;

		int vnum = get_attrib_count(MESH_ATTR_VERTEX);
		for(int i=0; i<vnum; i+=3) {
			Vec3 tmp = verts[i + 2];
			verts[i + 2] = verts[i + 1];
			verts[i + 1] = tmp;
		}
	}
}

void Mesh::flip_normals()
{
	Vec3 *normals = (Vec3*)get_attrib_data(MESH_ATTR_NORMAL);
	if(!normals) return;

	int num = get_attrib_count(MESH_ATTR_NORMAL);
	for(int i=0; i<num; i++) {
		normals[i] = -normals[i];
	}
}

void Mesh::explode()
{
	if(!is_indexed()) return;	// no vertex sharing is possible in unindexed meshes

	unsigned int *indices = get_index_data();
	assert(indices);

	int idxnum = get_index_count();
	int nnverts = idxnum;

	nverts = nnverts;

	for(int i=0; i<NUM_MESH_ATTR; i++) {
		if(!has_attrib(i)) continue;

		const float *srcbuf = get_attrib_data(i);

		float *tmpbuf = new float[nnverts * vattr[i].nelem];
		float *dstptr = tmpbuf;
		for(int j=0; j<idxnum; j++) {
			unsigned int idx = indices[j];
			const float *srcptr = srcbuf + idx * vattr[i].nelem;

			for(int k=0; k<vattr[i].nelem; k++) {
				*dstptr++ = *srcptr++;
			}
		}
		set_attrib_data(i, vattr[i].nelem, nnverts, tmpbuf);
		delete [] tmpbuf;
	}

	ibo_valid = false;
	idata_valid = false;
	idata.clear();

	nfaces = idxnum / 3;
}

void Mesh::calc_face_normals()
{
	const Vec3 *varr = (Vec3*)get_attrib_data(MESH_ATTR_VERTEX);
	Vec3 *narr = (Vec3*)get_attrib_data(MESH_ATTR_NORMAL);
	if(!varr) {
		return;
	}

	if(is_indexed()) {
		const unsigned int *idxarr = get_index_data();

		for(unsigned int i=0; i<nfaces; i++) {
			Triangle face(i, varr, idxarr);
			face.calc_normal();

			for(int j=0; j<3; j++) {
				unsigned int idx = *idxarr++;
				narr[idx] = face.normal;
			}
		}
	} else {
		// non-indexed
		int nfaces = nverts / 3;

		for(int i=0; i<nfaces; i++) {
			Triangle face(varr[0], varr[1], varr[2]);
			face.calc_normal();

			narr[0] = narr[1] = narr[2] = face.normal;
			varr += vattr[MESH_ATTR_VERTEX].nelem;
			narr += vattr[MESH_ATTR_NORMAL].nelem;
		}
	}
}

/*
int Mesh::add_bone(XFormNode *bone)
{
	int idx = bones.size();
	bones.push_back(bone);
	return idx;
}

const XFormNode *Mesh::get_bone(int idx) const
{
	if(idx < 0 || idx >= (int)bones.size()) {
		return 0;
	}
	return bones[idx];
}

int Mesh::get_bones_count() const
{
	return (int)bones.size();
}
*/

bool Mesh::pre_draw() const
{
	cur_sdr = 0;
	glGetIntegerv(GL_CURRENT_PROGRAM, &cur_sdr);

	((Mesh*)this)->update_buffers();

	if(!vattr[MESH_ATTR_VERTEX].vbo_valid) {
		fprintf(stderr, "%s: invalid vertex buffer\n", __FUNCTION__);
		return false;
	}

	if(cur_sdr && use_custom_sdr_attr) {
		// rendering with shaders
		if(global_sdr_loc[MESH_ATTR_VERTEX] == -1) {
			fprintf(stderr, "%s: shader attribute location for vertices unset\n", __FUNCTION__);
			return false;
		}

		for(int i=0; i<NUM_MESH_ATTR; i++) {
			int loc = global_sdr_loc[i];
			if(loc >= 0 && vattr[i].vbo_valid) {
				glBindBuffer(GL_ARRAY_BUFFER, vattr[i].vbo);
				glVertexAttribPointer(loc, vattr[i].nelem, GL_FLOAT, GL_FALSE, 0, 0);
				glEnableVertexAttribArray(loc);
			}
		}
	} else {
#ifndef GL_ES_VERSION_2_0
		// rendering with fixed-function (not available in GLES2)
		glBindBuffer(GL_ARRAY_BUFFER, vattr[MESH_ATTR_VERTEX].vbo);
		glVertexPointer(vattr[MESH_ATTR_VERTEX].nelem, GL_FLOAT, 0, 0);
		glEnableClientState(GL_VERTEX_ARRAY);

		if(vattr[MESH_ATTR_NORMAL].vbo_valid) {
			glBindBuffer(GL_ARRAY_BUFFER, vattr[MESH_ATTR_NORMAL].vbo);
			glNormalPointer(GL_FLOAT, 0, 0);
			glEnableClientState(GL_NORMAL_ARRAY);
		}
		if(vattr[MESH_ATTR_TEXCOORD].vbo_valid) {
			glBindBuffer(GL_ARRAY_BUFFER, vattr[MESH_ATTR_TEXCOORD].vbo);
			glTexCoordPointer(vattr[MESH_ATTR_TEXCOORD].nelem, GL_FLOAT, 0, 0);
			glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		}
		if(vattr[MESH_ATTR_COLOR].vbo_valid) {
			glBindBuffer(GL_ARRAY_BUFFER, vattr[MESH_ATTR_COLOR].vbo);
			glColorPointer(vattr[MESH_ATTR_COLOR].nelem, GL_FLOAT, 0, 0);
			glEnableClientState(GL_COLOR_ARRAY);
		}
		if(vattr[MESH_ATTR_TEXCOORD2].vbo_valid) {
			glClientActiveTexture(GL_TEXTURE1);
			glBindBuffer(GL_ARRAY_BUFFER, vattr[MESH_ATTR_TEXCOORD2].vbo);
			glTexCoordPointer(vattr[MESH_ATTR_TEXCOORD2].nelem, GL_FLOAT, 0, 0);
			glEnableClientState(GL_TEXTURE_COORD_ARRAY);
			glClientActiveTexture(GL_TEXTURE0);
		}
#endif
	}
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	return true;
}

void Mesh::draw() const
{
	if(!pre_draw()) return;

	if(ibo_valid) {
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glDrawElements(GL_TRIANGLES, nfaces * 3, GL_UNSIGNED_INT, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	} else {
		glDrawArrays(GL_TRIANGLES, 0, nverts);
	}

	post_draw();
}

void Mesh::post_draw() const
{
	if(cur_sdr && use_custom_sdr_attr) {
		// rendered with shaders
		for(int i=0; i<NUM_MESH_ATTR; i++) {
			int loc = global_sdr_loc[i];
			if(loc >= 0 && vattr[i].vbo_valid) {
				glDisableVertexAttribArray(loc);
			}
		}
	} else {
#ifndef GL_ES_VERSION_2_0
		// rendered with fixed-function
		glDisableClientState(GL_VERTEX_ARRAY);
		if(vattr[MESH_ATTR_NORMAL].vbo_valid) {
			glDisableClientState(GL_NORMAL_ARRAY);
		}
		if(vattr[MESH_ATTR_TEXCOORD].vbo_valid) {
			glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		}
		if(vattr[MESH_ATTR_COLOR].vbo_valid) {
			glDisableClientState(GL_COLOR_ARRAY);
		}
		if(vattr[MESH_ATTR_TEXCOORD2].vbo_valid) {
			glClientActiveTexture(GL_TEXTURE1);
			glDisableClientState(GL_TEXTURE_COORD_ARRAY);
			glClientActiveTexture(GL_TEXTURE0);
		}
#endif
	}
}

void Mesh::draw_wire() const
{
	if(!pre_draw()) return;

	((Mesh*)this)->update_wire_ibo();

	int num_faces = get_poly_count();
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, wire_ibo);
	glDrawElements(GL_LINES, num_faces * 6, GL_UNSIGNED_INT, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	post_draw();
}

void Mesh::draw_vertices() const
{
	if(!pre_draw()) return;

	glDrawArrays(GL_POINTS, 0, nverts);

	post_draw();
}

void Mesh::draw_normals() const
{
#ifdef USE_OLDGL
	int cur_sdr = 0;
	glGetIntegerv(GL_CURRENT_PROGRAM, &cur_sdr);

	Vec3 *varr = (Vec3*)get_attrib_data(MESH_ATTR_VERTEX);
	Vec3 *norm = (Vec3*)get_attrib_data(MESH_ATTR_NORMAL);
	if(!varr || !norm) {
		return;
	}

	glBegin(GL_LINES);
	if(cur_sdr && use_custom_sdr_attr) {
		int vert_loc = global_sdr_loc[MESH_ATTR_VERTEX];
		if(vert_loc < 0) {
			glEnd();
			return;
		}

		for(size_t i=0; i<nverts; i++) {
			glVertexAttrib3f(vert_loc, varr[i].x, varr[i].y, varr[i].z);
			Vec3 end = varr[i] + norm[i] * vis_vecsize;
			glVertexAttrib3f(vert_loc, end.x, end.y, end.z);
		}
	} else {
		for(size_t i=0; i<nverts; i++) {
			glVertex3f(varr[i].x, varr[i].y, varr[i].z);
			Vec3 end = varr[i] + norm[i] * vis_vecsize;
			glVertex3f(end.x, end.y, end.z);
		}
	}
	glEnd();
#endif	// USE_OLDGL
}

void Mesh::draw_tangents() const
{
#ifdef USE_OLDGL
	int cur_sdr = 0;
	glGetIntegerv(GL_CURRENT_PROGRAM, &cur_sdr);

	Vec3 *varr = (Vec3*)get_attrib_data(MESH_ATTR_VERTEX);
	Vec3 *tang = (Vec3*)get_attrib_data(MESH_ATTR_TANGENT);
	if(!varr || !tang) {
		return;
	}

	glBegin(GL_LINES);
	if(cur_sdr && use_custom_sdr_attr) {
		int vert_loc = global_sdr_loc[MESH_ATTR_VERTEX];
		if(vert_loc < 0) {
			glEnd();
			return;
		}

		for(size_t i=0; i<nverts; i++) {
			glVertexAttrib3f(vert_loc, varr[i].x, varr[i].y, varr[i].z);
			Vec3 end = varr[i] + tang[i] * vis_vecsize;
			glVertexAttrib3f(vert_loc, end.x, end.y, end.z);
		}
	} else {
		for(size_t i=0; i<nverts; i++) {
			glVertex3f(varr[i].x, varr[i].y, varr[i].z);
			Vec3 end = varr[i] + tang[i] * vis_vecsize;
			glVertex3f(end.x, end.y, end.z);
		}
	}
	glEnd();
#endif	// USE_OLDGL
}

void Mesh::get_aabbox(Vec3 *vmin, Vec3 *vmax) const
{
	if(!aabb_valid) {
		((Mesh*)this)->calc_aabb();
	}
	*vmin = aabb.min;
	*vmax = aabb.max;
}

const AABox &Mesh::get_aabbox() const
{
	if(!aabb_valid) {
		((Mesh*)this)->calc_aabb();
	}
	return aabb;
}

float Mesh::get_bsphere(Vec3 *center, float *rad) const
{
	if(!bsph_valid) {
		((Mesh*)this)->calc_bsph();
	}
	*center = bsph.center;
	*rad = bsph.radius;
	return bsph.radius;
}

const Sphere &Mesh::get_bsphere() const
{
	if(!bsph_valid) {
		((Mesh*)this)->calc_bsph();
	}
	return bsph;
}

/// static function
void Mesh::set_intersect_mode(unsigned int mode)
{
	Mesh::intersect_mode = mode;
}

/// static function
unsigned int Mesh::get_intersect_mode()
{
	return Mesh::intersect_mode;
}

/// static function
void Mesh::set_vertex_select_distance(float dist)
{
	Mesh::vertex_sel_dist = dist;
}

/// static function
float Mesh::get_vertex_select_distance()
{
	return Mesh::vertex_sel_dist;
}

bool Mesh::intersect(const Ray &ray, HitPoint *hit) const
{
	assert((Mesh::intersect_mode & (ISECT_VERTICES | ISECT_FACE)) != (ISECT_VERTICES | ISECT_FACE));

	const Vec3 *varr = (Vec3*)get_attrib_data(MESH_ATTR_VERTEX);
	const Vec3 *narr = (Vec3*)get_attrib_data(MESH_ATTR_NORMAL);
	if(!varr) {
		return false;
	}
	const unsigned int *idxarr = get_index_data();

	// first test with the bounding box
	AABox box;
	get_aabbox(&box.min, &box.max);
	if(!box.intersect(ray)) {
		return false;
	}

	HitPoint nearest_hit;
	nearest_hit.dist = FLT_MAX;
	nearest_hit.data = 0;

	if(Mesh::intersect_mode & ISECT_VERTICES) {
		// we asked for "intersections" with the vertices of the mesh
		long nearest_vidx = -1;
		float thres_sq = Mesh::vertex_sel_dist * Mesh::vertex_sel_dist;

		for(unsigned int i=0; i<nverts; i++) {

			if((Mesh::intersect_mode & ISECT_FRONT) && dot(narr[i], ray.dir) > 0) {
				continue;
			}

			// project the vertex onto the ray line
			float t = dot(varr[i] - ray.origin, ray.dir);
			Vec3 vproj = ray.origin + ray.dir * t;

			float dist_sq = length_sq(vproj - varr[i]);
			if(dist_sq < thres_sq) {
				if(!hit) {
					return true;
				}
				if(t < nearest_hit.dist) {
					nearest_hit.dist = t;
					nearest_vidx = i;
				}
			}
		}

		if(nearest_vidx != -1) {
			hitvert = varr[nearest_vidx];
			nearest_hit.data = &hitvert;
		}

	} else {
		// regular intersection test with polygons

		for(unsigned int i=0; i<nfaces; i++) {
			Triangle face(i, varr, idxarr);

			// ignore back-facing polygons if the mode flags include ISECT_FRONT
			if((Mesh::intersect_mode & ISECT_FRONT) && dot(face.get_normal(), ray.dir) > 0) {
				continue;
			}

			HitPoint fhit;
			if(face.intersect(ray, hit ? &fhit : 0)) {
				if(!hit) {
					return true;
				}
				if(fhit.dist < nearest_hit.dist) {
					nearest_hit = fhit;
					hitface = face;
				}
			}
		}
	}

	if(nearest_hit.data) {
		if(hit) {
			*hit = nearest_hit;

			// if we are interested in the mesh and not the faces set obj to this
			if(Mesh::intersect_mode & ISECT_FACE) {
				hit->data = &hitface;
			} else if(Mesh::intersect_mode & ISECT_VERTICES) {
				hit->data = &hitvert;
			} else {
				hit->data = (void*)this;
			}
		}
		return true;
	}
	return false;
}


// texture coordinate manipulation
void Mesh::texcoord_apply_xform(const Mat4 &xform)
{
	if(!has_attrib(MESH_ATTR_TEXCOORD)) {
		return;
	}

	for(unsigned int i=0; i<nverts; i++) {
		Vec4 tc = get_attrib(MESH_ATTR_TEXCOORD, i);
		set_attrib(MESH_ATTR_TEXCOORD, i, xform * tc);
	}
}

void Mesh::texcoord_gen_plane(const Vec3 &norm, const Vec3 &tang)
{
	if(!nverts) return;

	if(!has_attrib(MESH_ATTR_TEXCOORD)) {
		// allocate texture coordinate attribute array
		set_attrib_data(MESH_ATTR_TEXCOORD, 2, nverts);
	}

	Vec3 n = normalize(norm);
	Vec3 b = normalize(cross(n, tang));
	Vec3 t = cross(b, n);

	for(unsigned int i=0; i<nverts; i++) {
		Vec3 pos = Vec3(get_attrib(MESH_ATTR_VERTEX, i));

		// distance along the tangent direction
		float u = dot(pos, t);
		// distance along the bitangent direction
		float v = dot(pos, b);

		set_attrib(MESH_ATTR_TEXCOORD, i, Vec4(u, v, 0, 1));
	}
}

void Mesh::texcoord_gen_box()
{
	if(!nverts || !has_attrib(MESH_ATTR_NORMAL)) return;

	if(!has_attrib(MESH_ATTR_TEXCOORD)) {
		// allocate texture coordinate attribute array
		set_attrib_data(MESH_ATTR_TEXCOORD, 2, nverts);
	}

	for(unsigned int i=0; i<nverts; i++) {
		Vec3 pos = Vec3(get_attrib(MESH_ATTR_VERTEX, i)) * 0.5 + Vec3(0.5, 0.5, 0.5);
		Vec3 norm = Vec3(get_attrib(MESH_ATTR_NORMAL, i));

		float abs_nx = fabs(norm.x);
		float abs_ny = fabs(norm.y);
		float abs_nz = fabs(norm.z);
		int dom = abs_nx > abs_ny && abs_nx > abs_nz ? 0 : (abs_ny > abs_nz ? 1 : 2);

		float uv[2], *uvptr = uv;
		for(int j=0; j<3; j++) {
			if(j == dom) continue;	// skip dominant axis

			*uvptr++ = pos[j];
		}
		set_attrib(MESH_ATTR_TEXCOORD, i, Vec4(uv[0], uv[1], 0, 1));
	}
}

void Mesh::texcoord_gen_cylinder()
{
	if(!nverts) return;

	if(!has_attrib(MESH_ATTR_TEXCOORD)) {
		// allocate texture coordinate attribute array
		set_attrib_data(MESH_ATTR_TEXCOORD, 2, nverts);
	}

	for(unsigned int i=0; i<nverts; i++) {
		Vec3 pos = Vec3(get_attrib(MESH_ATTR_VERTEX, i));

		float rho = sqrt(pos.x * pos.x + pos.z * pos.z);
		float theta = rho == 0.0 ? 0.0 : atan2(pos.z / rho, pos.x / rho);

		float u = theta / (2.0 * M_PI) + 0.5;
		float v = pos.y;

		set_attrib(MESH_ATTR_TEXCOORD, i, Vec4(u, v, 0, 1));
	}
}


bool Mesh::dump(const char *fname) const
{
	FILE *fp = fopen(fname, "wb");
	if(fp) {
		bool res = dump(fp);
		fclose(fp);
		return res;
	}
	return false;
}

bool Mesh::dump(FILE *fp) const
{
	if(!has_attrib(MESH_ATTR_VERTEX)) {
		return false;
	}

	fprintf(fp, "VERTEX ATTRIBUTES\n");
	static const char *label[] = { "pos", "nor", "tan", "tex", "col", "bw", "bid", "tex2" };
	static const char *elemfmt[] = { 0, " %s(%g)", " %s(%g, %g)", " %s(%g, %g, %g)", " %s(%g, %g, %g, %g)", 0 };

	for(int i=0; i<(int)nverts; i++) {
		fprintf(fp, "%5u:", i);
		for(int j=0; j<NUM_MESH_ATTR; j++) {
			if(has_attrib(j)) {
				Vec4 v = get_attrib(j, i);
				int nelem = vattr[j].nelem;
				fprintf(fp, elemfmt[nelem], label[j], v.x, v.y, v.z, v.w);
			}
		}
		fputc('\n', fp);
	}

	if(is_indexed()) {
		const unsigned int *idx = get_index_data();
		int numidx = get_index_count();
		int numtri = numidx / 3;
		assert(numidx % 3 == 0);

		fprintf(fp, "FACES\n");

		for(int i=0; i<numtri; i++) {
			fprintf(fp, "%5d: %d %d %d\n", i, idx[0], idx[1], idx[2]);
			idx += 3;
		}
	}
	return true;
}

bool Mesh::dump_obj(const char *fname) const
{
	FILE *fp = fopen(fname, "wb");
	if(fp) {
		bool res = dump_obj(fp);
		fclose(fp);
		return res;
	}
	return false;
}

bool Mesh::dump_obj(FILE *fp, int voffs) const
{
	if(!has_attrib(MESH_ATTR_VERTEX)) {
		return false;
	}

	for(int i=0; i<(int)nverts; i++) {
		Vec4 v = get_attrib(MESH_ATTR_VERTEX, i);
		fprintf(fp, "v %f %f %f\n", v.x, v.y, v.z);
	}

	if(has_attrib(MESH_ATTR_NORMAL)) {
		for(int i=0; i<(int)nverts; i++) {
			Vec4 v = get_attrib(MESH_ATTR_NORMAL, i);
			fprintf(fp, "vn %f %f %f\n", v.x, v.y, v.z);
		}
	}

	if(has_attrib(MESH_ATTR_TEXCOORD)) {
		for(int i=0; i<(int)nverts; i++) {
			Vec4 v = get_attrib(MESH_ATTR_TEXCOORD, i);
			fprintf(fp, "vt %f %f\n", v.x, v.y);
		}
	}

	if(is_indexed()) {
		const unsigned int *idxptr = get_index_data();
		int numidx = get_index_count();
		int numtri = numidx / 3;
		assert(numidx % 3 == 0);

		for(int i=0; i<numtri; i++) {
			fputc('f', fp);
			for(int j=0; j<3; j++) {
				unsigned int idx = *idxptr++ + 1 + voffs;
				fprintf(fp, " %u/%u/%u", idx, idx, idx);
			}
			fputc('\n', fp);
		}
	} else {
		int numtri = nverts / 3;
		unsigned int idx = 1 + voffs;
		for(int i=0; i<numtri; i++) {
			fputc('f', fp);
			for(int j=0; j<3; j++) {
				fprintf(fp, " %u/%u/%u", idx, idx, idx);
				++idx;
			}
			fputc('\n', fp);
		}
	}
	return true;
}

// ------ private member functions ------

void Mesh::calc_aabb()
{
	// the cast is to force calling the const version which doesn't invalidate
	if(!((const Mesh*)this)->get_attrib_data(MESH_ATTR_VERTEX)) {
		return;
	}

	aabb.min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	aabb.max = -aabb.min;

	for(unsigned int i=0; i<nverts; i++) {
		Vec4 v = get_attrib(MESH_ATTR_VERTEX, i);
		for(int j=0; j<3; j++) {
			if(v[j] < aabb.min[j]) {
				aabb.min[j] = v[j];
			}
			if(v[j] > aabb.max[j]) {
				aabb.max[j] = v[j];
			}
		}
	}
	aabb_valid = true;
}

void Mesh::calc_bsph()
{
	// the cast is to force calling the const version which doesn't invalidate
	if(!((const Mesh*)this)->get_attrib_data(MESH_ATTR_VERTEX)) {
		return;
	}

	Vec3 v;
	bsph.center = Vec3(0, 0, 0);

	// first find the center
	for(unsigned int i=0; i<nverts; i++) {
		v = Vec3(get_attrib(MESH_ATTR_VERTEX, i));
		bsph.center += v;
	}
	bsph.center /= (float)nverts;

	bsph.radius = 0.0f;
	for(unsigned int i=0; i<nverts; i++) {
		v = Vec3(get_attrib(MESH_ATTR_VERTEX, i));
		float dist_sq = length_sq(v - bsph.center);
		if(dist_sq > bsph.radius) {
			bsph.radius = dist_sq;
		}
	}
	bsph.radius = sqrt(bsph.radius);

	bsph_valid = true;
}

void Mesh::update_buffers()
{
	for(int i=0; i<NUM_MESH_ATTR; i++) {
		if(has_attrib(i) && !vattr[i].vbo_valid) {
			glBindBuffer(GL_ARRAY_BUFFER, vattr[i].vbo);
			glBufferData(GL_ARRAY_BUFFER, nverts * vattr[i].nelem * sizeof(float), &vattr[i].data[0], GL_STATIC_DRAW);
			vattr[i].vbo_valid = true;
		}
	}
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	if(idata_valid && !ibo_valid) {
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, nfaces * 3 * sizeof(unsigned int), &idata[0], GL_STATIC_DRAW);
		ibo_valid = true;
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void Mesh::update_wire_ibo()
{
	update_buffers();

	if(wire_ibo_valid) {
		return;
	}

	if(!wire_ibo) {
		glGenBuffers(1, &wire_ibo);
	}

	int num_faces = get_poly_count();

	unsigned int *wire_idxarr = new unsigned int[num_faces * 6];
	unsigned int *dest = wire_idxarr;

	if(ibo_valid) {
		// we're dealing with an indexed mesh
		const unsigned int *idxarr = ((const Mesh*)this)->get_index_data();

		for(int i=0; i<num_faces; i++) {
			*dest++ = idxarr[0];
			*dest++ = idxarr[1];
			*dest++ = idxarr[1];
			*dest++ = idxarr[2];
			*dest++ = idxarr[2];
			*dest++ = idxarr[0];
			idxarr += 3;
		}
	} else {
		// not an indexed mesh ...
		for(int i=0; i<num_faces; i++) {
			int vidx = i * 3;
			*dest++ = vidx;
			*dest++ = vidx + 1;
			*dest++ = vidx + 1;
			*dest++ = vidx + 2;
			*dest++ = vidx + 2;
			*dest++ = vidx;
		}
	}

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, wire_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_faces * 6 * sizeof(unsigned int), wire_idxarr, GL_STATIC_DRAW);
	delete [] wire_idxarr;
	wire_ibo_valid = true;
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


// ------ class Triangle ------
Triangle::Triangle()
{
	normal_valid = false;
	id = -1;
}

Triangle::Triangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
{
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;
	normal_valid = false;
	id = -1;
}

Triangle::Triangle(int n, const Vec3 *varr, const unsigned int *idxarr)
{
	if(idxarr) {
		v[0] = varr[idxarr[n * 3]];
		v[1] = varr[idxarr[n * 3 + 1]];
		v[2] = varr[idxarr[n * 3 + 2]];
	} else {
		v[0] = varr[n * 3];
		v[1] = varr[n * 3 + 1];
		v[2] = varr[n * 3 + 2];
	}
	normal_valid = false;
	id = n;
}

void Triangle::calc_normal()
{
	normal = normalize(cross(v[1] - v[0], v[2] - v[0]));
	normal_valid = true;
}

const Vec3 &Triangle::get_normal() const
{
	if(!normal_valid) {
		((Triangle*)this)->calc_normal();
	}
	return normal;
}

void Triangle::transform(const Mat4 &xform)
{
	v[0] = xform * v[0];
	v[1] = xform * v[1];
	v[2] = xform * v[2];
	normal_valid = false;
}

void Triangle::draw() const
{
	Vec3 n[3];
	n[0] = n[1] = n[2] = get_normal();

	int vloc = Mesh::get_attrib_location(MESH_ATTR_VERTEX);
	int nloc = Mesh::get_attrib_location(MESH_ATTR_NORMAL);

	glEnableVertexAttribArray(vloc);
	glVertexAttribPointer(vloc, 3, GL_FLOAT, GL_FALSE, 0, &v[0].x);
	glVertexAttribPointer(nloc, 3, GL_FLOAT, GL_FALSE, 0, &n[0].x);

	glDrawArrays(GL_TRIANGLES, 0, 3);

	glDisableVertexAttribArray(vloc);
	glDisableVertexAttribArray(nloc);
}

void Triangle::draw_wire() const
{
	static const int idxarr[] = {0, 1, 1, 2, 2, 0};
	int vloc = Mesh::get_attrib_location(MESH_ATTR_VERTEX);

	glEnableVertexAttribArray(vloc);
	glVertexAttribPointer(vloc, 3, GL_FLOAT, GL_FALSE, 0, &v[0].x);

	glDrawElements(GL_LINES, 6, GL_UNSIGNED_INT, idxarr);

	glDisableVertexAttribArray(vloc);
}

Vec3 Triangle::calc_barycentric(const Vec3 &pos) const
{
	Vec3 norm = get_normal();

	float area_sq = fabs(dot(cross(v[1] - v[0], v[2] - v[0]), norm));
	if(area_sq < 1e-5) {
		return Vec3(0, 0, 0);
	}

	float asq0 = fabs(dot(cross(v[1] - pos, v[2] - pos), norm));
	float asq1 = fabs(dot(cross(v[2] - pos, v[0] - pos), norm));
	float asq2 = fabs(dot(cross(v[0] - pos, v[1] - pos), norm));

	return Vec3(asq0 / area_sq, asq1 / area_sq, asq2 / area_sq);
}

bool Triangle::intersect(const Ray &ray, HitPoint *hit) const
{
	Vec3 normal = get_normal();

	float ndotdir = dot(ray.dir, normal);
	if(fabs(ndotdir) < 1e-4) {
		return false;
	}

	Vec3 vertdir = v[0] - ray.origin;
	float t = dot(normal, vertdir) / ndotdir;

	Vec3 pos = ray.origin + ray.dir * t;
	Vec3 bary = calc_barycentric(pos);

	if(bary.x + bary.y + bary.z > 1.00001) {
		return false;
	}

	if(hit) {
		hit->dist = t;
		hit->pos = ray.origin + ray.dir * t;
		hit->normal = normal;
		hit->data = (void*)this;
	}
	return true;
}
