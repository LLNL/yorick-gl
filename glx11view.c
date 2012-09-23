/*
 * $Id: glx11view.c,v 1.1 2005-09-18 22:07:54 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#undef ONE_SIDED

#include <stdlib.h>
#include <stdio.h>
#include "glBasic.h"
#include "glcode.h"
#include "glfunc.h"
#include "glMouse.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

extern void MySwapBuffers(void);

#define SPHERE_LIST 1
#define DRAW_LIST 2

int have_sphere= 0;

extern glWinProp *glCurrWin3d;
extern glWinProp *glWin3dList[8];

/* //////////////////////////////////////////////////////////////// */

void yglPrepList(void)
{
  if(!glCurrWin3d) return;
  if(glCurrWin3d->use_list) {
    if( (glCurrWin3d->seq_num > 0) && (glCurrWin3d->seq_num > glCurrWin3d->list_num) ) {
      if(glCurrWin3d->have_gl_list) glDeleteLists(glCurrWin3d->the_gl_list, 1);
      glCurrWin3d->have_gl_list= 0;
      glNewList(glCurrWin3d->the_gl_list, GL_COMPILE);
    }
  }
}

void yglPrepDraw(glWinProp *theWin3d)
{
  glWinProp *res;

  if(!theWin3d) {
    char *tmp;
    tmp= getenv("DISPLAY");
    res= yglMakWin(tmp, 500, 500, "3D window 0");
    DEMAND(res, "failed to create 3D window")
    /* reserve an OpenGL display list for this window */
    glCurrWin3d->the_gl_list = glGenLists(1);
	theWin3d= glCurrWin3d;
    glWin3dList[0]= theWin3d;
  }
  DEMAND(theWin3d, "Unable to make OpenGL window")

  /* Make the rendering context current */
  yglMakeCurrent(theWin3d);

  /* Clear the window with current clearing color */
  glClearColor(theWin3d->back_red, theWin3d->back_green, 
               theWin3d->back_blue, theWin3d->back_alpha);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

#ifdef ONE_SIDED
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
#else
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
#endif
}

void yglFinFrame(void)
{
  /* all drawing complete. swap the back buffer to the screen */
  glFlush();
  /* Swap our scene to the front */
  MySwapBuffers();

  /* Allow other rendering contexts to co-exist */
  yglMakeCurrent(0);
}

void yglFinCache(void)
{
  /* close the display list at the end of "cached" mode drawing */
  if(!glCurrWin3d) return;
  if(glCurrWin3d->use_list) {
    if(!glCurrWin3d->have_gl_list) {
      glEndList();
      glCurrWin3d->have_gl_list= 1;
      glCurrWin3d->list_num= glCurrWin3d->seq_num;
    }
    glCallList(glCurrWin3d->the_gl_list);
  }
}

void yglFinDirect(void)
{
  /* nothing to do at the end of direct mode drawing */
}

void yglForceDraw(void)
{
  /* If this could be done in the "windows way" I would just
     post an invalidate. However, yorick does not drop into
     the event loop until it is completely idle, so instead 
     call a special OnDraw() clone function.
  */
  if(glCurrWin3d) yglDoDraw(glCurrWin3d);
}

void yglDoDraw(glWinProp *theWin3d)
{
  /* Make the rendering context current
     This call doesn't do anything useful when printing, because
     the only display list that exists is tied to the on screen window.
  */
  yglMakeCurrent(theWin3d);

  if(theWin3d->have_gl_list) {
    glCallList(theWin3d->the_gl_list);
    /* Flush drawing commands */
    glFlush();
  }
	
  /* Swap our scene to the front */
  MySwapBuffers();

  /* Allow other rendering contexts to co-exist */
  yglMakeCurrent(0);

  /* clear any pending redraw requests */
  /*  ValidateRect(NULL); */
}


void yglGetPixels(long nx, long ny, unsigned char *pix)
{
  GLint align_old;

  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  /* Make the rendering context current */
  yglMakeCurrent(glCurrWin3d);

  /* Read the screen buffer 
     WARNING: the pixel alignment must be set properly or the 
     call may write to bytes outside the array and give an
     improper phase for pixels in the array. */
  glReadBuffer(GL_FRONT);
  glGetIntegerv(GL_PACK_ALIGNMENT,&align_old);
  glPixelStorei(GL_PACK_ALIGNMENT,1);
  glReadPixels(0, 0, nx, ny, GL_RGB, GL_UNSIGNED_BYTE, pix);
  glPixelStorei(GL_PACK_ALIGNMENT,align_old);
  ygl_fpemask(1);
}

void yglPutPixels(long nx, long ny, unsigned char *pix)
{
  if(!glCurrWin3d) return;
  /* Make the rendering context current */
  ygl_fpemask(0);
  yglMakeCurrent(glCurrWin3d);

  /* save coord xform */
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0.0, (GLfloat) nx, 0.0, (GLfloat)ny, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();	/* save existing modelview matrix */
  glLoadIdentity();
  /* position at the lower left corner and draw the pixels. */
  glRasterPos2i(0, 0);
  glDrawPixels(nx, ny, GL_RGB, GL_UNSIGNED_BYTE, pix);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  ygl_fpemask(1);
}
