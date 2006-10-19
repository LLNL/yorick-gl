/*
 * $Id: glcode.h,v 1.2 2006-10-19 14:48:19 dhmunro Exp $
 * Declarations for external OpenGL module.  These functions are
 * defined in glcode.c and are called appropriately by the CView
 * derived classes.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __GLCODE__
#define __GLCODE__

#ifdef WIN32
/* WARNING: It appears that <windows.h> must be included
   before <GL/gl.h>. !!!! */
#include <windows.h>
/* #include <GL/glaux.h> */
#endif

#include "dlist3d.h"
#include "play.h"
#include "playgl.h"

#ifdef MACINTOSH
#include <gl.h>
#else
#include <GL/gl.h>
#endif

#ifdef __cplusplus

extern "C" {

#endif

extern int alpha_pass;

typedef struct g_callbacks g_callbacks;
struct g_callbacks {
  char *name;
  void (*gl_on_expose)(void *c, int *xy);
  void (*gl_on_destroy)(void *c);
  void (*gl_on_resize)(void *c,int w,int h);
  void (*gl_on_focus)(void *c,int in);
  void (*gl_on_key)(void *c,int k,int md);
  void (*gl_on_click)(void *c,int b,int md,int x,int y, unsigned long ms);
  void (*gl_on_motion)(void *c,int md,int x,int y);
  void (*gl_on_deselect)(void *c);
} ;

typedef struct glInnerWinProp glInnerWinProp;
struct glInnerWinProp {
  g_callbacks *on;
  struct glWinProp *topwin;
} ;

#ifdef WIN32
  typedef void (APIENTRY *MYTEX3DFUNC) (GLenum target, GLint level, 
                 GLenum internalFormat, GLsizei width, 
                 GLsizei height, GLsizei depth, GLint border, GLenum format,
                 GLenum type, const GLvoid *pixels );
#else
  typedef int MYTEX3DFUNC;
#endif
typedef void (*glBindTextureProc)(GLenum target, GLint level);

typedef struct glWinProp glWinProp;
struct glWinProp {
  g_callbacks *on;
  struct glInnerWinProp *inner;
  char *name;
  glInnerWinProp innerWin;
  
  p_glwin *gl_win;
  p_win *top_win;
  p_scr *s;
  int dirty;
  float back_red, back_green, back_blue, back_alpha;
  /* The "cage" is the three coord. planes at the "back" of the scene 
   as viewed from the current viewpoint. The planes have one color
   and the grid lines another. The planes are chosen from (xmin,xmax),
   (ymin,ymax), and (zmin,zmax) pairs. If cage_style is negative,
   these planes will be the limits of the data, if positive, they will
   be values set by the user, if zero no cage will be drawn.
   */
  float cage_red, cage_green, cage_blue, cage_alpha;
  float grid_red, grid_green, grid_blue, grid_alpha;
  float cage_xmin, cage_xmax, cage_ymin, cage_ymax;
  float cage_zmin, cage_zmax;
  long cage_style;
  int num_xgrid, num_ygrid, num_zgrid;
  long cage_seq_num;
  long cage_state;

  /* Lighting components */
  GLfloat ambientLight[4];
  GLfloat diffuseLight[4];
  GLfloat specularLight[4];
  /* if the last parameter is zero, the light is at infinity */
  GLfloat positionLight[4];
  GLint light_model;
  GLfloat mat_spec;
  GLint shade_model;
  GLint cull_mode;
  GLint poly_sides;
  GLint poly_mode;
  GLint mat_color;

  /* parameters describing the color of the light */
  GLfloat curr_ambient[4];
  GLfloat curr_diffuse[4];
  GLfloat curr_specular[4];
  /* if the last parameter is zero, the light is at infinity */
  GLfloat curr_position[4];
  GLint curr_light_model;

  /* current specular color */
  GLfloat curr_mat_spec[4];
  GLint curr_shade_model;
  GLint curr_cull_mode;
  GLint curr_poly_sides;
  GLint curr_poly_mode;
  GLint curr_mat_color;

  double eye[3];
  double center[3];
  double up[3];
  double view[3];
  double viewdist;
  float fov;

  long width, hite;

  float curr_alpha;
  int have_gl_list, the_gl_list, always_show_obj, object_on;
  int cursor_action;
  long use_list;
  long use_array;
  long seq_num;
  long list_num;
  yBox3D boxAll;
  long BoxSeqNum;
  int hascubetex, hasTex3d, hasTex3dExt, hasTexExt, tex3dChecked;
  MYTEX3DFUNC myglTexImage3D_ptr;
  int myGL_TEXTURE_3D, myGL_PROXY_TEXTURE_3D;
  glBindTextureProc myglBindTexture3D;
};

  void yglInitRC(void *pData);
  void yglResize(glWinProp *theWin3d, GLsizei h, GLsizei w);
  void draw_plane(float pt1[3], float pt2[3], float pt3[3], int num1, int num2);
  void yglWireSphere(int list_num, double radius);
  void yglForceUpdateLight(void);
  void yglUpdateLight(void);
  void yglUpdateProperties(void);
  void yglForceUpdateProperties(void);
  void yglChekError(char *x);
  int yglNewList(void);
  int yglCloseList(void);
  int yglDrawList(int num);
  int yglDeleteList(int num);
  int yglArrlim3D(long nvert, double *arr, double *lims);
  
  extern glWinProp *glCurrWin3d;

  void yglInitWin3d(glWinProp *theWin);

  /* provide GLU wrappers (see glustub.c) */
#undef gluLookAt
#define gluLookAt my_gluLookAt
#undef gluPerspective
#define gluPerspective my_gluPerspective
#undef gluErrorString
#define gluErrorString my_gluErrorString
  void gluLookAt(GLdouble eyex, GLdouble eyey, GLdouble eyez,
                 GLdouble centerx, GLdouble centery, GLdouble centerz,
                 GLdouble upx, GLdouble upy, GLdouble upz);
  void gluPerspective(GLdouble fovy, GLdouble aspect,
                      GLdouble near, GLdouble far);
  const GLubyte *gluErrorString(GLenum errorCode);

#ifdef __cplusplus
}
#endif

/* Easier way to specify color for Windows bigots */
#define yglRGB(x, y, z)	glColor3ub((GLubyte)x, (GLubyte)y, (GLubyte)z)

#endif /* Include/Define */
