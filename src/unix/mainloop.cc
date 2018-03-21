#include <errno.h>
#include <unistd.h>
#include <sys/select.h>
#include <GL/freeglut.h>
#include <GL/glx.h>
#include <X11/Xlib.h>
#include "mainloop.h"

static bool redraw_pending;

void main_loop(struct resman *rman)
{
	Display *dpy = glXGetCurrentDisplay();
	int xfd = ConnectionNumber(dpy);

	for(;;) {
		if(!redraw_pending) {
			int num_rman_fds;
			int *rman_fds = resman_get_wait_fds(rman, &num_rman_fds);

			fd_set rdset;
			int maxfd = xfd;

			FD_ZERO(&rdset);
			FD_SET(xfd, &rdset);

			for(int i=0; i<num_rman_fds; i++) {
				FD_SET(rman_fds[i], &rdset);
				if(rman_fds[i] > maxfd) maxfd = rman_fds[i];
			}

			int res;
			while((res = select(maxfd + 1, &rdset, 0, 0, 0)) == -1 && errno == EINTR);

			if((res > 0 && !FD_ISSET(xfd, &rdset)) || res > 1) {
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
