/*
 * $Id: glx11setup.c,v 1.2 2006-03-25 03:12:29 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "glcode.h"
#include "glfunc.h"
#include "glMouse.h"
#include "glBasic.h"
#include "binio.h"
#include "pstdlib.h"
#include <stdio.h>

extern void g_disconnect(p_scr *s);
extern int p_wincount(p_scr *s);

static void gl_on_expose(void *c, int *xy);
static void gl_on_destroy(void *c);
static void gl_on_resize(void *c,int w,int h);
static void gl_on_focus(void *c,int in);
static void gl_on_key(void *c,int k,int md);
static void gl_on_click(void *c,int b,int md,int x,int y, unsigned long ms);
static void gl_on_motion(void *c,int md,int x,int y);
static void gl_on_deselect(void *c);

extern void MySwapBuffers(void);
extern void yglInitlights(void);
#if 0
extern void myReshape(int w, int h);
#endif
extern void yglInitRC(void *pData);
extern void ygl_update_3d(void);

extern int the_gl_list;

static int ygl_mouse_pos_x= 0;
static int ygl_mouse_pos_y= 0;
static int ygl_mouse_is_down= 0;
static int scr_no_win= 0;

glWinProp *glCurrWin3d= 0;
glWinProp *glWin3dList[8];

/* ------------------------------------------------------------------------ */

extern void (*g_on_keyline)(char *msg);  /* gist/xfancy.c */
static char my_msg[96];
static int my_msglen = 0;

static g_callbacks gl_x11_on = {
  "GL top level", 0, gl_on_destroy, gl_on_resize, gl_on_focus,
  gl_on_key, gl_on_click, gl_on_motion, gl_on_deselect };

static g_callbacks gl_glx_on = {
  "GL window", gl_on_expose, 0, 0, 0, 0, 0, 0, 0 };


int yglWin3d(int num, int w, int h)
{
  int i;
  glWinProp *res;
  char titlstr[80];

  if(!glCurrWin3d) {
    /* make sure the list of window pointers is empty */
    for(i= 0; i <= 7; i++) glWin3dList[i]= 0;
  }
  if(num < 0 || num > 7) return 1;  /* bad window number */
  if(glWin3dList[num]) {
    /* change the current window */
    glCurrWin3d= glWin3dList[num];
  } else {
    /* create a new 3D window */
    sprintf(titlstr, "3D window %d", num);
    res= yglMakWin(0, w, h, titlstr);
    if(!res) {
      /* failed to create the 3D window */
      return 2;
	}
    /* reserve an OpenGL display list for this window */
    glCurrWin3d->the_gl_list = glGenLists(1);
    glWin3dList[num]= glCurrWin3d;
  }
  return 0;  /* indicate success */
}

int yglWinKill3d(int num)
{
  glWinProp *theWin3d;

  if(num < 0 || num > 7) return 1;  /* bad window number */
  if(glWin3dList[num]) {
    theWin3d= glWin3dList[num];
    shutdown3d(theWin3d);
    glWin3dList[num]= 0;
    /* update the current 3D window only if necessary */
    if(glCurrWin3d == theWin3d) resetcurrwin3d();
  } else {
    /* tried to delete a non-existent 3D window */
    return 2;
  }
  return 0;  /* indicate success */
}

int isWin3d(glWinProp *win)
{
  /* return 1 if this is the top level window for a 3D "canvas" */
  if(win && win->on == &gl_x11_on) return 1;
  /* return 2 if this is the glx window for a 3D "canvas" */
  if(win && win->on == &gl_glx_on) return 2;
  return 0;
}

int yglCurrWin3d(void)
{
  int i;

  for(i= 0; i <= 7; i++) {
    if(glWin3dList[i] == glCurrWin3d) {
      return i;
    }
  }
  return 0;
}

glWinProp *yglMakWin(char *displayName, int w, int h, char *title)
{
  p_glwin *gl_win = 0;
  p_win *top_win = 0;
  extern p_scr *g_connect(char *displayName);  /* gist/xbasic.h */
  extern int gist_input_hint;
  glWinProp *theWin3d;
  glWinProp *oldWin3d= glCurrWin3d;
  glInnerWinProp *theInner;
  p_scr *scr_gl;

  scr_gl = g_connect(displayName);
  if (!scr_gl) return 0;
  theWin3d= (glWinProp *) p_malloc(sizeof(glWinProp));
  if(!theWin3d) return 0;
  theInner= &(theWin3d->innerWin);
  theWin3d->inner= theInner;
  theWin3d->on= &gl_x11_on;
  theWin3d->s= scr_gl;
  theInner->topwin= theWin3d;
  theInner->on= &gl_glx_on;
  yglInitWin3d(theWin3d);
  glCurrWin3d= theWin3d;

  top_win = p_window(scr_gl, w, h, title, P_BG,
                    (gist_input_hint?0:P_NOKEY) | P_RGBMODEL, theWin3d);
  if (!top_win) {
    glCurrWin3d= oldWin3d;
    p_free(theWin3d);
    return 0;
  }
  gl_win = p_glcreate(top_win, w, h, 0, 0, theInner);
  if (!gl_win) {
    p_win *w = top_win;
    glCurrWin3d= oldWin3d;
    top_win = 0;
    p_destroy(w);
    p_free(theWin3d);
    return 0;
  }
  theWin3d->gl_win= gl_win;
  theWin3d->top_win= top_win;
  theWin3d->width= w;
  theWin3d->hite= h;

  /* Finish Set-up */
  p_glcurrent(gl_win);
  yglInitRC(0);

  return theWin3d;
}

int shutdown3d(glWinProp *theWin3d)
{
  int num;
  p_glwin *wg = theWin3d->gl_win;
  p_win *topwin = theWin3d->top_win;
  p_scr *s;

  if(isWin3d(theWin3d) != 1) {
    /* should never see anything but a top level window for the 3D canvas, so return */
    return -1;
  }
  /* set a flag indicating that there are (for the moment)
     no connections to this display */
/*  s= (theWin3d->top_win)->s; */
  s= theWin3d->s;
  if (s && !p_wincount(s)) scr_no_win= 1;

  num= winnum3d(theWin3d);
  if(num < 0) {
    /* unknown window, so do nothing */
    return -2;
  }
  theWin3d->top_win = 0;
  theWin3d->dirty= 0; /* stop any redraw attempts */
  if (wg) p_gldestroy(wg);
  theWin3d->gl_win = 0;
  if (topwin) p_destroy(topwin);
  /* In Windows, the 3D window might have been destroyed 
     by an asynchronous on_destroy() triggered by
     the p_gldestroy() call. Protect against a double free() */
  p_free(theWin3d);
  glWin3dList[num]= 0;
  return 0;  /* indicate success */
}

void resetcurrwin3d(void)
{
  int i;

  glCurrWin3d= 0; /* indicate no 3D window left */
  for(i= 7; i >= 0; i--) {
    if(glWin3dList[i]) {
      glCurrWin3d= glWin3dList[i];
      break;
    }
  }
}

int winnum3d(glWinProp *theWin3d)
{
  int i;

  for(i= 7; i >= 0; i--) {
    if(glWin3dList[i] == theWin3d) {
      return i;
    }
  }
  return -1;  /* indicate not found */
}

void MySwapBuffers(void)
{
  if (glCurrWin3d && glCurrWin3d->gl_win) p_glswap(glCurrWin3d->gl_win);
}

void yglMakeCurrent(glWinProp *theWin3d)
{
  if (theWin3d && theWin3d->gl_win) p_glcurrent(theWin3d->gl_win);
}

void yglInitWin3d(glWinProp *theWin)
{
  theWin->name= 0;
  theWin->gl_win= 0;
  theWin->top_win= 0;

  theWin->dirty= 0;
  theWin->back_red= 0.0f;
  theWin->back_green= 0.3f;
  theWin->back_blue= 0.8f;
  theWin->back_alpha= 1.0f;
  theWin->cage_red= 0.3f;
  theWin->cage_green= 0.3f;
  theWin->cage_blue= 0.6f;
  theWin->cage_alpha= 1.0f;
  theWin->grid_red= 0.5f;
  theWin->grid_green= 0.5f;
  theWin->grid_blue= 0.3f;
  theWin->grid_alpha= 1.0f;
  theWin->cage_xmin= 0.0;
  theWin->cage_xmax= 1.0;
  theWin->cage_ymin= 0.0;
  theWin->cage_ymax= 1.0;
  theWin->cage_zmin= 0.0;
  theWin->cage_zmax= 1.0;
  theWin->cage_style= 0;
  theWin->num_xgrid= 3;
  theWin->num_ygrid=3;
  theWin->num_zgrid= 3;

  /* Lighting components */
  theWin->ambientLight[0]= 0.2f;
  theWin->ambientLight[1]= 0.2f;
  theWin->ambientLight[2]= 0.2f;
  theWin->ambientLight[3]= 1.0f;
  theWin->diffuseLight[0]= 0.6f;
  theWin->diffuseLight[1]= 0.6f;
  theWin->diffuseLight[2]= 0.6f;
  theWin->diffuseLight[3]= 1.0f;
  theWin->specularLight[0]= 1.0f;
  theWin->specularLight[1]= 1.0f;
  theWin->specularLight[2]= 1.0f;
  theWin->specularLight[3]= 1.0f;
  /* if the last parameter is zero, the light is at infinity */
  theWin->positionLight[0]= 0.0f;
  theWin->positionLight[1]= 0.0f;
  theWin->positionLight[2]= 1.0f;
  theWin->positionLight[3]= 0.0f;
  theWin->light_model= GL_TRUE;
  theWin->mat_spec= 1.0f;
  theWin->shade_model= GL_SMOOTH;
  theWin->cull_mode= 0;
  theWin->poly_sides= GL_FRONT_AND_BACK;
  theWin->poly_mode= GL_FILL;
  theWin->mat_color= GL_AMBIENT_AND_DIFFUSE;

  /* parameters describing the color of the light */
  theWin->curr_ambient[0]= 0.2f;
  theWin->curr_ambient[1]= 0.2f;
  theWin->curr_ambient[2]= 0.2f;
  theWin->curr_ambient[3]= 1.0f;
  theWin->curr_diffuse[0]= 0.6f;
  theWin->curr_diffuse[1]= 0.6f;
  theWin->curr_diffuse[2]= 0.6f;
  theWin->curr_diffuse[3]= 1.0f;
  theWin->curr_specular[0]= 1.0f;
  theWin->curr_specular[1]= 1.0f;
  theWin->curr_specular[2]= 1.0f;
  theWin->curr_specular[3]= 1.0f;
  /* if the last parameter is zero, the light is at infinity */
  theWin->curr_position[0]= 0.0f;
  theWin->curr_position[1]= 0.0f;
  theWin->curr_position[2]= 1.0f;
  theWin->curr_position[3]= 0.0f;
  theWin->curr_light_model= GL_TRUE;

  /* current specular color */
  theWin->curr_mat_spec[0]= 1.0f;
  theWin->curr_mat_spec[1]= 1.0f;
  theWin->curr_mat_spec[2]= 1.0f;
  theWin->curr_mat_spec[3]= 1.0f;
  theWin->curr_shade_model= GL_SMOOTH;
  theWin->curr_cull_mode= 0;
  theWin->curr_poly_sides= GL_FRONT_AND_BACK;
  theWin->curr_poly_mode= GL_FILL;
  theWin->curr_mat_color= GL_AMBIENT_AND_DIFFUSE;

  theWin->eye[0]= 10.0;
  theWin->eye[1]= 0.0;
  theWin->eye[2]= 0.0;
  theWin->center[0]= 0.0;
  theWin->center[1]= 0.0;
  theWin->center[2]= 0.0;
  theWin->up[0]= 0.0;
  theWin->up[1]= 0.0;
  theWin->up[2]= 1.0;
  theWin->view[0]= 1.0;
  theWin->view[1]= 0.0;
  theWin->view[2]= 0.0;
  theWin->viewdist= 10.0;
  theWin->fov= 35.0;
  theWin->width=500;
  theWin->hite=500;
  theWin->curr_alpha=0.5;
  theWin->have_gl_list= 0;
  theWin->use_list= 1;
  theWin->use_array= 0;
  theWin->seq_num= 0;
  theWin->list_num= -1;
  theWin->cage_seq_num= 0;
  theWin->cage_state= -1;
  theWin->cursor_action= ACTION_ROTATE;
  theWin->BoxSeqNum= -1;
  theWin->tex3dChecked= 0;
  theWin->hasTex3d= 0;
  theWin->hasTex3dExt= 0;
  theWin->myglTexImage3D_ptr= 0;
  theWin->myGL_TEXTURE_3D= 0;
  theWin->myGL_PROXY_TEXTURE_3D= 0;
  theWin->myglBindTexture3D= 0;
  theWin->hascubetex= -1;
  theWin->hasTexExt= 0;
}



void ygl_update_3d(void)
{
  int i;

  /* check 3D windows for dirty bits and redraw as needed */
  if(scr_no_win) {
    /* on X11 systems, need to disconnect from the display
       when there are no longer any open windows on it */
    g_disconnect(0);
    scr_no_win= 0;
  }
  for(i= 0; i <= 7; i++) {
    if(glWin3dList[i]) {
      /* Need to draw the window if an event like a mouse action, 
         an expose or a redraw occurred. Also need to draw if yorgl's 
         display list has been altered. */
      if(glWin3dList[i]->dirty || (glWin3dList[i]->seq_num > glWin3dList[i]->list_num) ) {
        yglDraw3d(glWin3dList[i]);
      }
    }
  }
}

static void gl_on_expose(void *c, int *xy)
{
  glWinProp *theWin3d;
  glInnerWinProp *theInner= (glInnerWinProp *)c;

  if(isWin3d((glWinProp *)theInner) != 2) {
    /* only need to respond to expose events on the inner window */
    return;
  }
  theWin3d= theInner->topwin;
  theWin3d->dirty= 1;
}

static void gl_on_destroy(void *c)
{
  glWinProp *theWin3d= (glWinProp *)c;

  if(isWin3d(theWin3d) != 1) {
    /* should never see anything but a top level window for the 3D canvas, so return */
    return;
  }
  shutdown3d(theWin3d);
  /* update the current 3D window only if necessary */
  if(glCurrWin3d == theWin3d) resetcurrwin3d();
}

static void gl_on_resize(void *c,int w,int h)
{
  glWinProp *theWin3d= (glWinProp *)c;

  if(isWin3d(theWin3d) != 1) {
    /* should never see anything but a top level window for the 3D canvas, so return */
    return;
  }
  if (theWin3d->gl_win) {
    yglResize(theWin3d, w, h);
    theWin3d->dirty= 1;
  }
}

static void gl_on_focus(void *c,int in)
{
  /* nothing to do */
}

static void gl_on_deselect(void *c)
{
  /* nothing to do */
}

#undef DBG_MOUSE
#ifdef DBG_MOUSE
#include "play.h"
static mouse_count= 0;
static char msg_buf[120];
#endif

/* b=1 left b=2 middle b=3 right */
static void gl_on_click(void *c,int b,int md,int x,int y, unsigned long ms)
{
  glWinProp *theWin3d= (glWinProp *)c;
  glWinProp *oldWin3d= glCurrWin3d;

  if(isWin3d(theWin3d) != 1) {
    /* should never see anything but a top level window for the 3D canvas, so return */
    return;
  }
  if (!theWin3d->gl_win) return;
  oldWin3d= glCurrWin3d;
  glCurrWin3d= theWin3d;
  if (md & (1<<(b+2))) {
#ifdef DBG_MOUSE
    sprintf(msg_buf, "Mouse released. %d motion events while down\n",
            mouse_count);
    p_stdout(msg_buf);
#endif
    DoLButtonUp(ygl_mouse_is_down, x, y, theWin3d);
    ygl_mouse_is_down = 0;
  } else if (!ygl_mouse_is_down) {
#ifdef DBG_MOUSE
    mouse_count= 0;
    p_stdout("Mouse button pressed\n");
#endif
    if (b==1) {
      if (md & P_SHIFT) b = 3;
      else if (md & P_CONTROL) b = 2;
    }
    ygl_mouse_is_down = b;
    DoLButtonDown(b, x, y);
    if (!glCurrWin3d->always_show_obj) glCurrWin3d->object_on= 0;
  }
  glCurrWin3d= oldWin3d;
}

static void gl_on_motion(void *c,int md,int x,int y)
{
  glWinProp *theWin3d= (glWinProp *)c;
  glWinProp *oldWin3d= glCurrWin3d;

  if(isWin3d(theWin3d) != 1) {
    /* should never see anything but a top level window for the 3D canvas, so return */
    return;
  }
  if (!theWin3d->gl_win) return;
#ifdef DBG_MOUSE
  mouse_count++;
#endif
  ygl_mouse_pos_x = x;
  ygl_mouse_pos_y = y;
  /* IGNORE mouse motion unless the left button is down */
  if (ygl_mouse_is_down) {
    /* carry out any necessary rotation */
    oldWin3d= glCurrWin3d;
    glCurrWin3d= theWin3d;
    yglMouseMove(ygl_mouse_is_down, x, y, theWin3d);
    glCurrWin3d= oldWin3d;
  }
}

static void gl_on_key(void *c,int k,int md)
{
  glWinProp *theWin3d= (glWinProp *)c;

  if(isWin3d(theWin3d) != 1) {
    /* should never see anything but a top level window for the 3D canvas, so return */
    return;
  }
  if (!theWin3d->gl_win) return;
  if (!my_msglen) my_msg[0] = '\0';
  if (k>=' ' && k<'\177') {
    /* append printing characters */
    if (my_msglen>=94) my_msglen = 0;  /* crude overflow handling */
    my_msg[my_msglen++] = k;
    my_msg[my_msglen] = '\0';
  } else if (k=='\177' || k=='\010') {
    /* delete or backspace kills char */
    if (my_msglen)
      my_msg[--my_msglen] = '\0';
  } else if (k=='\027') {
    /* C-w kills word */
    int n = my_msglen;
    char c = n? my_msg[n-1] : '\0';
    if (c<'0' || (c>'9' && c<'A') || (c>'Z' && c<'a' && c!='_') || c>'z') {
      if (my_msglen)
        my_msg[--my_msglen] = '\0';
    } else {
      while (--n) {
        c = my_msg[n-1];
        if (c<'0' || (c>'9' && c<'A') || (c>'Z' && c<'a' && c!='_') ||
            c>'z') break;
      }
      my_msg[n] = '\0';
      my_msglen = n;
    }
  } else if (k=='\025') {
    /* C-u kills line */
    my_msg[0] = '\0';
    my_msglen = 0;
  } else if (k=='\012' || k=='\015') {
    /* linefeed or carriage return sends line to interpreter */
    int n = my_msglen;
    my_msg[n] = '\012';
    my_msg[n+1] = '\0';
    p_stdout(my_msg);
    my_msg[n] = '\0';
    if (g_on_keyline) g_on_keyline(my_msg);
    my_msg[0] = '\0';
    my_msglen = 0;
  }
  /* RedrawMessage(fxe); no feedback for now */
}

void yglInitlights(void)
{
  GLfloat ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat position[] = { 0.0, 0.0, 2.0, 1.0 };
  GLfloat mat_diffuse[] = { 0.6f, 0.6f, 0.6f, 1.0f };
  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_shininess[] = { 50.0 };

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_POSITION, position);

#if 0
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
#else
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
#endif
}

#if 0
void myReshape(int w, int h)
{
  if(w < 20) w= 20;
  if(h < 20) h= 20;
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if (w <= h)
    glOrtho(-4.0, 4.0, -4.0*(GLfloat)h/(GLfloat)w, 
	    4.0*(GLfloat)h/(GLfloat)w, -4.0, 4.0);
  else
    glOrtho(-4.0*(GLfloat)w/(GLfloat)h, 
	    4.0*(GLfloat)w/(GLfloat)h, -4.0, 4.0, -4.0, 4.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glCurrWin3d->width= w;
  glCurrWin3d->hite= h;
}
#endif
