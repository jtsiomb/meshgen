#include <stdio.h>
#include <GL/freeglut.h>
#include <windows.h>
#include "mainloop.h"

static bool redraw_pending;

void main_loop(struct resman *rman)
{
	for(;;) {
		if(!redraw_pending) {
			int num_handles;
			unsigned int idx;
			void **handles = resman_get_wait_handles(rman, &num_handles);

			idx = MsgWaitForMultipleObjects(num_handles, handles, 0, INFINITE, QS_ALLEVENTS);
			if(idx == WAIT_FAILED) {
				unsigned int err = GetLastError();
				fprintf(stderr, "failed to wait for events: %u\n", err);
			}

			if((int)idx < num_handles) {
				glutPostRedisplay();
			}
		}

		redraw_pending = false;
		glutMainLoopEvent();
	}
}

void post_redisplay()
{
	redraw_pending = true;
	glutPostRedisplay();
}
