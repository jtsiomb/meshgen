#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include "opengl.h"
#include <GL/glut.h>
#include <gmath/gmath.h>
#include "mesh.h"
#include "opt.h"
#include "loadgen.h"

static bool init();
static void cleanup();
static void display();
static void draw_grid(float delta, int count);
static void reshape(int x, int y);
static void keydown(unsigned char key, int x, int y);
static void mouse(int bn, int st, int x, int y);
static void motion(int x, int y);
static bool reload();

static Mesh *mesh;
static bool opt_grid = true;

static float cam_theta, cam_phi = 25, cam_dist = 5;
static Vec3 cam_pos = Vec3(0, 1, 0);

static int prev_x, prev_y;
static bool bnstate[8];


int main(int argc, char **argv)
{
	glutInit(&argc, argv);

	if(!parse_options(argc, argv)) {
		return 1;
	}

	glutInitWindowSize(1280, 800);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE | GLUT_MULTISAMPLE);
	glutCreateWindow("Procedural mesh editor");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keydown);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);

	if(!init()) {
		return 1;
	}
	atexit(cleanup);

	glutMainLoop();
	return 0;
}


static bool init()
{
	glewInit();

	glClearColor(0.12, 0.12, 0.12, 1);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);

	glEnable(GL_MULTISAMPLE);

	reload();

	return true;
}

static void cleanup()
{
}

static void set_light(int idx, const Vec3 &pos, const Vec3 &color)
{
	unsigned int lt = GL_LIGHT0 + idx;
	float posv[] = { pos.x, pos.y, pos.z, 1 };
	float colv[] = { color.x, color.y, color.z, 1 };

	glEnable(lt);
	glLightfv(lt, GL_POSITION, posv);
	glLightfv(lt, GL_DIFFUSE, colv);
	glLightfv(lt, GL_SPECULAR, colv);
}


static void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0, 0, -cam_dist);
	glRotatef(cam_phi, 1, 0, 0);
	glRotatef(cam_theta, 0, 1, 0);
	glTranslatef(-cam_pos.x, -cam_pos.y, -cam_pos.z);

	static const Vec3 lpos[] = { Vec3(-50, 75, 100), Vec3(100, 0, 30), Vec3(-10, -10, 60) };
	set_light(0, lpos[0], Vec3(1.0, 0.8, 0.7) * 0.8);
	set_light(1, lpos[1], Vec3(0.6, 0.7, 1.0) * 0.6);
	set_light(2, lpos[2], Vec3(0.8, 1.0, 0.8) * 0.3);

	if(opt_grid) {
		draw_grid(0.5, 10);
	}

	if(mesh) {
		mesh->draw();
	}

	glutSwapBuffers();
	assert(glGetError() == GL_NO_ERROR);
}

static void draw_grid(float delta, int count)
{
	float max_dist = delta * count;

	glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT);

	glDisable(GL_LIGHTING);

	glLineWidth(2.0);

	glBegin(GL_LINES);
	glColor3f(1, 0, 0);
	glVertex3f(0, 0, -max_dist);
	glVertex3f(0, 0, max_dist);
	glColor3f(0, 1, 0);
	glVertex3f(-max_dist, 0, 0);
	glVertex3f(max_dist, 0, 0);
	glEnd();

	glLineWidth(1.0);

	glBegin(GL_LINES);
	glColor3f(0.4, 0.4, 0.4);

	float d = delta;
	for(int i=0; i<count; i++) {
		glVertex3f(-d, 0, -max_dist);
		glVertex3f(-d, 0, max_dist);
		glVertex3f(d, 0, -max_dist);
		glVertex3f(d, 0, max_dist);
		glVertex3f(-max_dist, 0, -d);
		glVertex3f(max_dist, 0, -d);
		glVertex3f(-max_dist, 0, d);
		glVertex3f(max_dist, 0, d);
		d += delta;
	}
	glEnd();

	glPopAttrib();
}

static void reshape(int x, int y)
{
	glViewport(0, 0, x, y);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(50, (float)x / (float)y, 0.5, 500.0);
}

static void keydown(unsigned char key, int x, int y)
{
	switch(key) {
	case 27:
		exit(0);

	case 'r':
	case 'R':
		reload();
		glutPostRedisplay();
		break;

	case 'w':
	case 'W':
		if(mesh) {
			if(mesh->dump_obj("mesh.obj")) {
				printf("Mesh saved: mesh.obj\n");
			}
		}
		break;

	case 'g':
		opt_grid = !opt_grid;
		glutPostRedisplay();
		break;
	}
}

static void mouse(int bn, int st, int x, int y)
{
	bnstate[bn - GLUT_LEFT] = st == GLUT_DOWN;
	prev_x = x;
	prev_y = y;
}

static void motion(int x, int y)
{
	int dx = x - prev_x;
	int dy = y - prev_y;
	prev_x = x;
	prev_y = y;

	if(!dx && !dy) return;

	if(bnstate[0]) {
		cam_theta += dx * 0.5;
		cam_phi += dy * 0.5;

		if(cam_phi < -90) cam_phi = -90;
		if(cam_phi > 90) cam_phi = 90;
		glutPostRedisplay();
	}
	if(bnstate[1]) {
		float panx = dx * 0.01;
		float pany = dy * 0.01;

		float theta = deg_to_rad(cam_theta);
		float phi = deg_to_rad(cam_phi);

		Mat4 cmat;
		cmat.pre_rotate_x(phi);
		cmat.pre_rotate_y(theta);

		Vec3 right = Vec3(panx, 0, 0) * cmat;
		Vec3 up = Vec3(0, pany, 0) * cmat;
		cam_pos += right - up;

		glutPostRedisplay();
	}
	if(bnstate[2]) {
		cam_dist += dy * 0.1;

		if(cam_dist < 0) cam_dist = 0;
		glutPostRedisplay();
	}
}


static bool reload()
{
	MeshGen *mgen = load_meshgen(opt.fname);
	if(!mgen) return false;

	Mesh *newmesh = new Mesh;
	if(!mgen->generate(newmesh, 0)) {
		fprintf(stderr, "mesh generation failed\n");
		delete newmesh;
		free_meshgen(mgen);
		return false;
	}
	free_meshgen(mgen);

	delete mesh;
	mesh = newmesh;
	return true;
}
