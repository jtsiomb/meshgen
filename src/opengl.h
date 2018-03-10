#ifndef OPENGL_H_
#define OPENGL_H_

#ifdef HAVE_CONFIG_H_
#include "config.h"
#endif

#if defined(IPHONE) || defined(__IPHONE__)
#include <OpenGLES/ES2/gl.h>

#define glClearDepth	glClearDepthf
#define GLDEF
#include "sanegl.h"

#elif defined(ANDROID)
#include <GLES2/gl2.h>
#include <GLES2/gl2ext.h>
#define GLDEF
#include "sanegl.h"

#else

#include <GL/glew.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif	/* __APPLE__ */

#endif	/* IPHONE */

#endif	/* OPENGL_H_ */
