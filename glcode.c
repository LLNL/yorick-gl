/*
 * $Id: glcode.c,v 1.2 2006-10-19 14:48:19 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "glcode.h"
#include "glfunc.h"
#include <math.h>
#include <stdio.h>

#define USE_TRI

/*
#define CHEK_ERROR(x)	gl_chek_error(x)
*/
#define CHEK_ERROR(x)

extern void (*g_on_idle)(void);
extern void ygl_update_3d(void);

void yglChekError(char *x)
{
  GLenum val;
  char msg[120];

  /* check for any OpenGl errors and clear the error flag */
  val= glGetError();
  if(val == GL_NO_ERROR) {
    return;
  } else if(val == GL_INVALID_ENUM) {
    sprintf(msg, "OpenGL error GL_INVALID_ENUM in %s\n", x);
  } else if(val == GL_INVALID_VALUE) {
    sprintf(msg, "OpenGL error GL_INVALID_VALUE in %s\n", x);
  } else if(val == GL_INVALID_OPERATION) {
    sprintf(msg, "OpenGL error GL_INVALID_OPERATION in %s\n", x);
  } else if(val == GL_STACK_OVERFLOW) {
    sprintf(msg, "OpenGL error GL_STACK_OVERFLOW in %s\n", x);
  } else if(val == GL_STACK_UNDERFLOW) {
    sprintf(msg, "OpenGL error GL_STACK_UNDERFLOW in %s\n", x);
  } else if(val == GL_OUT_OF_MEMORY) {
    sprintf(msg, "OpenGL error GL_OUT_OF_MEMORY in %s\n", x);
  } else {
    sprintf(msg, "GLU error in %s\n", x);
  }
  p_stderr(msg);
}

/*  //////////////////////////////////////////////////////////////
   Do any initialization of the rendering context here, such as
   setting background colors, setting up lighting, or performing
   preliminary calculations.
*/
void yglInitRC(void *pData)
{
  glEnable(GL_DEPTH_TEST);	/* Hidden surface removal */
  /* The following should be GL defaults, but set them anyway */
  /* don't use alpha to decide whether to draw */
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_STENCIL_TEST);	/* disable stencil buffer */
  glDisable(GL_BLEND);	/* disable blending of incoming colors
                           and the color buffer */
  glDisable(GL_DITHER);	/* diasble dithering */

  glClearColor(glCurrWin3d->back_red, glCurrWin3d->back_green, glCurrWin3d->back_blue, glCurrWin3d->back_alpha);
  CHEK_ERROR("GLInitRC after setting clear color");

  glEnable(GL_LIGHTING);
  yglForceUpdateLight();
  glEnable(GL_LIGHT0);
  CHEK_ERROR("GLInitRC after enabling light");

  glEnable(GL_COLOR_MATERIAL);
  yglForceUpdateProperties();
  /* set pixel storage alignment for definiteness */
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  /* Draw in green */
  yglRGB(0, 255, 0);
  /* install the function that draws "dirty" 3D windows if one has not
     already been set */
  if(!g_on_idle) {
    g_on_idle= ygl_update_3d;
  }
  CHEK_ERROR("GLInitRC");
}

extern void p_glresize(p_glwin *w, int width, int height, int x, int y);

/*  //////////////////////////////////////////////////////////////
   Called by Windows when it receives the WM_SIZE message.
   Put any code needed here to recalc the viewing volume and
   viewport info.
*/
void yglResize(glWinProp *theWin3d, GLsizei w, GLsizei h)
{
#define USE_PERSPECTIVE
#ifdef USE_PERSPECTIVE
  GLdouble theAspectRatio;
#endif

  /* Prevent a divide by zero */
  if(h <= 20) h = 20;
  if(w <= 20) w = 20;
  theWin3d->width= w;
  theWin3d->hite= h;

  /* update the size of the "inner window" (someday the outer 2D
     window might be larger than the 3D window to provide room
     for buttons. */
  p_glresize(theWin3d->gl_win, w, h, 0, 0);
  /* Set Viewport to window dimensions */
  glViewport(0, 0, w, h);

  glMatrixMode(GL_PROJECTION);
  /* Reset coordinate system */
  glLoadIdentity();

  /* Establish clipping volume (left, right, bottom, top, near, far) */
  theAspectRatio= (GLdouble)theWin3d->width/(GLdouble)theWin3d->hite;
  /* Use perspective viewing with the requested field-of-view and
     an aspect ratio matching the window. Set the near and far clipping
     planes in z. Use the full opening angle instead of the half
     opening angle used internally in this package. */
  gluPerspective (2.0*theWin3d->fov, theAspectRatio, 
                  theWin3d->viewdist/25.0, 16.0*theWin3d->viewdist);
  /* Set up transformation. */ 
  glMatrixMode (GL_MODELVIEW);
  gluLookAt(theWin3d->eye[0], theWin3d->eye[1], theWin3d->eye[2], 
            theWin3d->center[0], theWin3d->center[1], 
            theWin3d->center[2], theWin3d->up[0], theWin3d->up[1], 
            theWin3d->up[2]);
  CHEK_ERROR("GLResize");
}

void draw_plane(float pt1[3], float pt2[3], float pt3[3], int num1, int num2)
{
  int i;
  float rgb[3], norm[3], lin1[3], lin2[3], pt4[3];
  float dx1, dy1, dz1, dx2, dy2, dz2;
  float ambi[]= {1.0f, 1.0f, 1.0f, 1.0f};
  double val;

  /* The three points are the first three corners of a rectangle.
     Draw the rectangle, num1 lines along the first edge, and
     num2 lines along the second edge.
  */

  if(alpha_pass) return;
  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();

  rgb[0]= glCurrWin3d->cage_red;
  rgb[1]= glCurrWin3d->cage_green;
  rgb[2]= glCurrWin3d->cage_blue;
  /* the normal is the cross product of the first two edges */
  dx1= pt2[0]-pt1[0];
  dy1= pt2[1]-pt1[1];
  dz1= pt2[2]-pt1[2];
  dx2= pt3[0]-pt2[0];
  dy2= pt3[1]-pt2[1];
  dz2= pt3[2]-pt2[2];
  pt4[0]= pt1[0]+dx2; pt4[1]= pt1[1]+dy2; pt4[2]= pt1[2]+dz2;
  norm[0]= dy1*dz2-dz1*dy2;
  norm[1]= dz1*dx2-dx1*dz2;
  norm[2]= dx1*dy2-dy1*dx2;
  val= 1.0/sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
  norm[0] *= (float)val; norm[1] *= (float)val; norm[2] *= (float)val;

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0, 1.0);
  glBegin(GL_POLYGON);
  glColor3fv(rgb);
  glNormal3fv(norm);
  glVertex3fv(pt1);
  glVertex3fv(pt2);
  glVertex3fv(pt3);
  glVertex3fv(pt4);
  glEnd();
  glDisable(GL_POLYGON_OFFSET_FILL);
  /* draw the grid lines */
  rgb[0]= glCurrWin3d->grid_red;
  rgb[1]= glCurrWin3d->grid_green;
  rgb[2]= glCurrWin3d->grid_blue;
  lin1[0]= pt1[0]; lin1[1]= pt1[1]; lin1[2]= pt1[2];
  lin2[0]= pt4[0]; lin2[1]= pt4[1]; lin2[2]= pt4[2];
  /*  this creates diffuse light not connected to any light source */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambi);
  /* turn off all light associated with the light source */
  glDisable(GL_LIGHT0);
  glBegin(GL_LINES);
  glColor3fv(rgb);
  for(i= 0; i <= num1+1; i++) {
    glVertex3fv(lin1);
    glVertex3fv(lin2);
    lin1[0] += dx1/(num1+1); lin1[1] += dy1/(num1+1); lin1[2] += dz1/(num1+1);
    lin2[0] += dx1/(num1+1); lin2[1] += dy1/(num1+1); lin2[2] += dz1/(num1+1);
  }
  glEnd();
  lin1[0]= pt1[0]; lin1[1]= pt1[1]; lin1[2]= pt1[2];
  lin2[0]= pt2[0]; lin2[1]= pt2[1]; lin2[2]= pt2[2];
  glBegin(GL_LINES);
  glColor3fv(rgb);
  for(i= 0; i <= num2+1; i++) {
    glVertex3fv(lin1);
    glVertex3fv(lin2);
    lin1[0] += dx2/(num2+1); lin1[1] += dy2/(num2+1); lin1[2] += dz2/(num2+1);
    lin2[0] += dx2/(num2+1); lin2[1] += dy2/(num2+1); lin2[2] += dz2/(num2+1);
  }
  glEnd();
  /* restore prior diffuse light not connected to any light source */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, glCurrWin3d->curr_ambient);
  glEnable(GL_LIGHT0);
}

void yglCells(long nx, long ny, float xyz[3][3], float *norm, 
              float *colr, long do_alpha)
{
  /* draw a cell array in 3D. The data is a grid of nx by ny cells.
     xyz is a 3 by 3 array. The first 3 is for x, y, and z.
     The second three is for the three corners of the cell array */
  long i, j, base, colrsiz;
  float x0, y0, z0, dx0, dy0, dz0, dx1, dy1, dz1;
  float xa, ya, za, xb, yb, zb, oldSpec;
  float blackvec[]= {0.0f, 0.0f, 0.0f, 0.0f};

  if(do_alpha && !alpha_pass) return;
  if(!do_alpha && alpha_pass) return;
  /* there are nx by ny cells to be drawn. the first two corners
     are along an "x-edge" and the 2nd and 3rd are along a "y-edge". */
  x0= xyz[0][0];
  y0= xyz[0][1];
  z0= xyz[0][2];
  dx0= ((xyz[1][0]-xyz[0][0])/nx);
  dy0= ((xyz[1][1]-xyz[0][1])/nx);
  dz0= ((xyz[1][2]-xyz[0][2])/nx);
  dx1= ((xyz[2][0]-xyz[1][0])/ny);
  dy1= ((xyz[2][1]-xyz[1][1])/ny);
  dz1= ((xyz[2][2]-xyz[1][2])/ny);

  oldSpec= yglGetMatSpec();
  yglSetMatSpec(0.0);  /* turn off specular highlights */
  yglUpdateProperties();
  /* set ambient and diffuse colors to black (they are tracking glColor) */
  glColor3f(0.0,0.0,0.0);
  if(do_alpha) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    colrsiz= 4;
  } else {
    glDisable(GL_BLEND);
    colrsiz= 3;
  }

  /* There are ny quad strips (one for each column). */
  for(j= 0; j < ny; j++) {
    xa= x0+j*dx1;
    xb= xa+dx1;
    ya= y0+j*dy1;
    yb= ya+dy1;
    za= z0+j*dz1;
    zb= za+dz1;
    glBegin(GL_QUAD_STRIP);
    glNormal3fv(norm);
    if(do_alpha) {
      /* there are nx quads in the strip. set nx coords in the looop
         and the last one after the loop */
      for(i= 0, base= colrsiz*j*nx; i < nx; i++) {
        /* NOTE: It makes no difference in local performance if the coords
           are converted to floats before being stored in the display list */
        glVertex3f(xa, ya, za);
        glVertex3f(xb, yb, zb);
        xa += dx0;
        xb += dx0;
        ya += dy0;
        yb += dy0;
        za += dz0;
        zb += dz0;
        glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,colr+colrsiz*i+base);
      }
    } else {
      /* there are nx quads in the strip. set nx coords in the looop
         and the last one after the loop */
      for(i= 0, base= colrsiz*j*nx; i < nx; i++) {
        /* NOTE: It makes no difference in local performance if the coords
           are converted to floats before being stored in the display list */
        glVertex3f(xa, ya, za);
        glVertex3f(xb, yb, zb);
        xa += dx0;
        xb += dx0;
        ya += dy0;
        yb += dy0;
        za += dz0;
        zb += dz0;
        glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,colr+colrsiz*i+base);
      }
    }
    glVertex3f(xa, ya, za);
    glVertex3f(xb, yb, zb);
    glEnd();
  }
  if(do_alpha) {
    glDisable(GL_BLEND);
  }
  yglSetMatSpec(oldSpec);
  /* turn emission color off */
  glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,blackvec);
  yglForceUpdateProperties();
}

void yglPlm(long nx, long ny, float *xyz, float *colr)
{
  long i, j;

  /* draw a logically 2D mesh in one color */
  if(nx <= 0 || ny <= 0) return;
  if(alpha_pass) return;

  /* NOTE: It makes no difference in local performance if the coords
     are converted to floats before being stored in the display list */
  /* bump the line thickness up so that it will show through
     the plane */
  for(i= 0; i < nx; i++) {
    glBegin(GL_LINE_STRIP);
    glColor3fv(colr);
    for(j= 0; j < ny; j++) {
      glVertex3fv(xyz+3*i+3*j*nx);
    }
    glEnd();
  }
  for(j= 0; j < ny; j++) {
    glBegin(GL_LINE_STRIP);
    glColor3fv(colr);
    for(i= 0; i < nx; i++) {
      glVertex3fv(xyz+3*i+3*j*nx);
    }
    glEnd();
  }
  CHEK_ERROR("yglPlm3d");
}

void yglPlf(long nx, long ny, float *xyz, float *colr)
{
  long i, j;
  float oldSpec;
  float blackvec[]= {0.0f, 0.0f, 0.0f, 0.0f};

  /* fill a logically 2D mesh */
  if(nx <= 0 || ny <= 0) return;
  if(alpha_pass) return;

  /* NOTE: It makes no difference in local performance if the coords
     are converted to floats before being stored in the display list */
  /* set ambient and diffuse colors to black (they are tracking glColor) */
  oldSpec= yglGetMatSpec();
  yglSetMatSpec(0.0);  /* turn off specular highlights */
  yglUpdateProperties();
  glColor4f(0.0,0.0,0.0, 0.0);
  for(j= 0; j < ny-1; j++) {
    for(i= 0; i < nx-1; i++) {
      glBegin(GL_POLYGON);
      glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,colr+4*i+4*j*(nx-1));
      glVertex3fv(xyz+3*i+3*j*nx);
      glVertex3fv(xyz+3*(i+1)+3*j*nx);
      glVertex3fv(xyz+3*(i+1)+3*(j+1)*nx);
      glVertex3fv(xyz+3*i+3*(j+1)*nx);
      glEnd();
    }
  }
  yglSetMatSpec(oldSpec);
  /* turn emission color off */
  glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,blackvec);
  yglForceUpdateProperties();
  CHEK_ERROR("yglPlf3d");
}

void yglSurf(long do_alpha, long nx, long ny, float *xyz, float *norm, float *colr)
{
  long i, j, base;

  /* fill a logically 2D mesh */
  if(nx <= 0 || ny <= 0) return;
  if(do_alpha && !alpha_pass) return;
  if(!do_alpha && alpha_pass) return;

  /* NOTE: It makes no difference in local performance if the coords
     are converted to floats before being stored in the display list */
  yglSetPolyMode(0);
  /* use smooth shading */
  yglSetShade(1);
  yglSetColorType(1);
  yglUpdateProperties();
  /* set ambient and diffuse colors to desired value */
  if(do_alpha) {
    glColor4fv(colr);
  } else {
    glColor3fv(colr);
  }
  for(j= 0; j < ny-1; j++) {
    glBegin(GL_QUAD_STRIP);
    for(i= 0; i < nx; i++) {
      base= 3*i+3*j*nx;
      glNormal3fv(norm+base);
      glVertex3fv(xyz +base);
      glNormal3fv(norm+base+3*nx);
      glVertex3fv(xyz +base+3*nx);
    }
    glEnd();
  }
  CHEK_ERROR("yglSurf3d");
}

void yglColrSurf(long do_alpha, long nx, long ny, float *xyz, float *norm, float *colr)
{
  long i, j, base;

  /* fill a logically 2D mesh with a color per vertex */
  if(nx <= 0 || ny <= 0) return;
  if(do_alpha && !alpha_pass) return;
  if(!do_alpha && alpha_pass) return;

  /* NOTE: It makes no difference in local performance if the coords
     are converted to floats before being stored in the display list */
  yglSetPolyMode(0);
  /* use smooth shading */
  yglSetShade(1);
/*
  yglSetColorType(1);
*/
  yglUpdateProperties();
  /* set ambient and diffuse colors to desired value */
  if(do_alpha) {
    for(j= 0; j < ny-1; j++) {
      glBegin(GL_QUAD_STRIP);
      for(i= 0; i < nx; i++) {
        base= 3*i+3*j*nx;
        glColor4fv(colr+4*i+4*j*nx);
        glNormal3fv(norm+base);
        glVertex3fv(xyz +base);
        glColor4fv(colr+4*i+4*(j+1)*nx);
        glNormal3fv(norm+base+3*nx);
        glVertex3fv(xyz +base+3*nx);
      }
      glEnd();
    }
  } else {
    for(j= 0; j < ny-1; j++) {
      glBegin(GL_QUAD_STRIP);
      for(i= 0; i < nx; i++) {
        base= 3*i+3*j*nx;
        glColor3fv(colr+base);
        glNormal3fv(norm+base);
        glVertex3fv(xyz +base);
        glColor3fv(colr+base+3*nx);
        glNormal3fv(norm+base+3*nx);
        glVertex3fv(xyz +base+3*nx);
      }
      glEnd();
    }
  }
  CHEK_ERROR("yglColrSurf3d");
}

void yglLineWidth(double width)
{
  ygl_fpemask(0);
  glLineWidth((float)width);
  ygl_fpemask(1);
}

void yglLines(long nvert, float *xyz, float *colr)
{
  long i;
  float ambi[]= {1.0f, 1.0f, 1.0f, 1.0f};

  /* draw the polyline */
  if(nvert <= 1) return;
  if(alpha_pass) return;

  /*  this creates diffuse light not connected to any light source */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambi);
  /* turn off all light associated with the light source */
  glDisable(GL_LIGHT0);
  glBegin(GL_LINE_STRIP);
  glColor3fv(colr);
  for(i= 0; i < nvert; i++) {
    glVertex3fv(xyz+3*i);
  }
  glEnd();
  /* restore prior diffuse light not connected to any light source */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, glCurrWin3d->curr_ambient);
  glEnable(GL_LIGHT0);
}

void yglPoints(long nvert, float *xyz, float *colr)
{
  long i;
  float ambi[]= {1.0f, 1.0f, 1.0f, 1.0f};

  /* draw the points */
  if(nvert <= 1) return;
  if(alpha_pass) return;

  /*  glPointSize(3.5); */

  /*  this creates diffuse light not connected to any light source */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambi);
  /* turn off all light associated with the light source */
  glDisable(GL_LIGHT0);
  glBegin(GL_POINTS);
  for(i= 0; i < nvert; i++) {
    glColor3fv(colr+3*i);
    glVertex3fv(xyz+3*i);
  }
  glEnd();
  /* restore prior diffuse light not connected to any light source */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, glCurrWin3d->curr_ambient);
  glEnable(GL_LIGHT0);
}

void yglFrontFace3d(long dir)
{
  ygl_fpemask(0);
  if(dir) glFrontFace(GL_CCW);
  else glFrontFace(GL_CW);
  ygl_fpemask(1);
}

void yglGetCenter3d(double *center)
{
  center[0]= glCurrWin3d->center[0];
  center[1]= glCurrWin3d->center[1];
  center[2]= glCurrWin3d->center[2];
}

void yglGetEye3d(double *eye)
{
  eye[0]= glCurrWin3d->eye[0];
  eye[1]= glCurrWin3d->eye[1];
  eye[2]= glCurrWin3d->eye[2];
}

void yglGetUp3d(double *up)
{
  up[0]= glCurrWin3d->up[0];
  up[1]= glCurrWin3d->up[1];
  up[2]= glCurrWin3d->up[2];
}

void yglSetLight3d(double ambient, double diffuse, double spec, 
                 double spower, double *sdir)
{
  if(glCurrWin3d->ambientLight[0] != (float)ambient) {
    glCurrWin3d->ambientLight[0] = (float)ambient;
    glCurrWin3d->ambientLight[1] = (float)ambient;
    glCurrWin3d->ambientLight[2] = (float)ambient;
  }
  if(glCurrWin3d->diffuseLight[0] != (float)diffuse) {
    glCurrWin3d->diffuseLight[0] = (float)diffuse;
    glCurrWin3d->diffuseLight[1] = (float)diffuse;
    glCurrWin3d->diffuseLight[2] = (float)diffuse;
  }
  if(glCurrWin3d->specularLight[0] != (float)spec) {
    glCurrWin3d->specularLight[0] = (float)spec;
    glCurrWin3d->specularLight[1] = (float)spec;
    glCurrWin3d->specularLight[2] = (float)spec;
  }
  glCurrWin3d->positionLight[0] = (float)sdir[0];
  glCurrWin3d->positionLight[1] = (float)sdir[1];
  glCurrWin3d->positionLight[2] = (float)sdir[2];
  CHEK_ERROR("set_gl_light");
}

void yglPrepContext(void)
{
  GLdouble theAspectRatio;

  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity();
  /* The light should be set immediately after loading
     the identity transform. All properties should
     be set because this is the start of a plot. */
  yglForceUpdateLight();

  /* Reset coordinate system */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* Establish clipping volume (left, right, bottom, top, near, far) */
  theAspectRatio= (GLdouble)glCurrWin3d->width/(GLdouble)glCurrWin3d->hite;
  /* Use perspective viewing with the requested field-of-view and
     an aspect ratio matching the window. Set the near and far clipping
     planes in z. Use the full opening angle instead of he half
     opening angle used internally in this package. */
  gluPerspective (2.0*glCurrWin3d->fov, theAspectRatio, 
                  glCurrWin3d->viewdist/25.0, 16.0*glCurrWin3d->viewdist);
  /* Set up transformation. */ 
  glMatrixMode (GL_MODELVIEW);
  gluLookAt(glCurrWin3d->eye[0], glCurrWin3d->eye[1], glCurrWin3d->eye[2], 
            glCurrWin3d->center[0], glCurrWin3d->center[1], 
            glCurrWin3d->center[2], glCurrWin3d->up[0], glCurrWin3d->up[1], 
            glCurrWin3d->up[2]);
  CHEK_ERROR("yglPrepContext");
  yglForceUpdateProperties();
}

void yglForceUpdateLight(void)
{
  glCurrWin3d->curr_ambient[0]= glCurrWin3d->ambientLight[0];
  glCurrWin3d->curr_ambient[1]= glCurrWin3d->ambientLight[1];
  glCurrWin3d->curr_ambient[2]= glCurrWin3d->ambientLight[2];
  glCurrWin3d->curr_diffuse[0]= glCurrWin3d->diffuseLight[0];
  glCurrWin3d->curr_diffuse[1]= glCurrWin3d->diffuseLight[1];
  glCurrWin3d->curr_diffuse[2]= glCurrWin3d->diffuseLight[2];
  glCurrWin3d->curr_specular[0]= glCurrWin3d->specularLight[0];
  glCurrWin3d->curr_specular[1]= glCurrWin3d->specularLight[1];
  glCurrWin3d->curr_specular[2]= glCurrWin3d->specularLight[2];
  glCurrWin3d->curr_position[0]= glCurrWin3d->positionLight[0];
  glCurrWin3d->curr_position[1]= glCurrWin3d->positionLight[1];
  glCurrWin3d->curr_position[2]= glCurrWin3d->positionLight[2];
  glCurrWin3d->curr_position[3]= glCurrWin3d->positionLight[3];
  glCurrWin3d->curr_light_model= glCurrWin3d->light_model;
  glLightfv(GL_LIGHT0,GL_AMBIENT,glCurrWin3d->curr_ambient);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,glCurrWin3d->curr_diffuse);
  glLightfv(GL_LIGHT0,GL_SPECULAR,glCurrWin3d->curr_specular);
  glLightfv(GL_LIGHT0,GL_POSITION,glCurrWin3d->curr_position);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, glCurrWin3d->curr_light_model);
  CHEK_ERROR("ForceUpdateLight");
}

void yglUpdateLight(void)
{
  if(glCurrWin3d->curr_ambient[0] != glCurrWin3d->ambientLight[0] || 
     glCurrWin3d->curr_ambient[1] != glCurrWin3d->ambientLight[1] || 
     glCurrWin3d->curr_ambient[2] != glCurrWin3d->ambientLight[2]) {
    glCurrWin3d->curr_ambient[0]= glCurrWin3d->ambientLight[0];
    glCurrWin3d->curr_ambient[1]= glCurrWin3d->ambientLight[1];
    glCurrWin3d->curr_ambient[2]= glCurrWin3d->ambientLight[2];
    glLightfv(GL_LIGHT0,GL_AMBIENT,glCurrWin3d->curr_ambient);
  }
  if(glCurrWin3d->curr_diffuse[0] != glCurrWin3d->diffuseLight[0] || 
     glCurrWin3d->curr_diffuse[1] != glCurrWin3d->diffuseLight[1] || 
     glCurrWin3d->curr_diffuse[2] != glCurrWin3d->diffuseLight[2]) {
    glCurrWin3d->curr_diffuse[0]= glCurrWin3d->diffuseLight[0];
    glCurrWin3d->curr_diffuse[1]= glCurrWin3d->diffuseLight[1];
    glCurrWin3d->curr_diffuse[2]= glCurrWin3d->diffuseLight[2];
    glLightfv(GL_LIGHT0,GL_DIFFUSE,glCurrWin3d->curr_diffuse);
  }
  if(glCurrWin3d->curr_specular[0] != glCurrWin3d->specularLight[0] || 
     glCurrWin3d->curr_specular[1] != glCurrWin3d->specularLight[1] || 
     glCurrWin3d->curr_specular[2] != glCurrWin3d->specularLight[2]) {
    glCurrWin3d->curr_specular[0]= glCurrWin3d->specularLight[0];
    glCurrWin3d->curr_specular[1]= glCurrWin3d->specularLight[1];
    glCurrWin3d->curr_specular[2]= glCurrWin3d->specularLight[2];
    glLightfv(GL_LIGHT0,GL_SPECULAR,glCurrWin3d->curr_specular);
  }
  /* not changing fourth element, so it remains "directional" or
     "positional" as initially set. */
  if(glCurrWin3d->curr_position[0] != glCurrWin3d->positionLight[0] || 
     glCurrWin3d->curr_position[1] != glCurrWin3d->positionLight[1] || 
     glCurrWin3d->curr_position[2] != glCurrWin3d->positionLight[2] || 
     glCurrWin3d->curr_position[3] != glCurrWin3d->positionLight[3]) {
    glCurrWin3d->curr_position[0]= glCurrWin3d->positionLight[0];
    glCurrWin3d->curr_position[1]= glCurrWin3d->positionLight[1];
    glCurrWin3d->curr_position[2]= glCurrWin3d->positionLight[2];
    glCurrWin3d->curr_position[3]= glCurrWin3d->positionLight[3];
    glLightfv(GL_LIGHT0,GL_POSITION,glCurrWin3d->curr_position);
  }
  if(glCurrWin3d->curr_light_model != glCurrWin3d->light_model) {
    glCurrWin3d->curr_light_model= glCurrWin3d->light_model;
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, glCurrWin3d->curr_light_model);
  }
  CHEK_ERROR("UpdateLight");
}

void yglSetPolyMode(int mode)
{
  if(!mode) glCurrWin3d->poly_mode= GL_FILL;
  else glCurrWin3d->poly_mode= GL_LINE;
}

void yglSetPolySides(int sides)
{
  if(sides) glCurrWin3d->poly_sides= GL_FRONT_AND_BACK;
  else glCurrWin3d->poly_sides= GL_FRONT;
}

void yglSetMatSpec(float spec)
{
  if(spec >= 0.0 && spec <= 1.0) glCurrWin3d->mat_spec= spec;
}

float yglGetMatSpec(void)
{
  return glCurrWin3d->mat_spec;
}

void yglSetShade(int smooth)
{
  if(smooth) glCurrWin3d->shade_model= GL_SMOOTH;
  else glCurrWin3d->shade_model= GL_FLAT;
}

void yglSetColorType(int colorType)
{
  /* specifies which sort of color is set by glColorMaterial */
  if(colorType) glCurrWin3d->mat_color= GL_AMBIENT_AND_DIFFUSE;
  else glCurrWin3d->mat_color= GL_EMISSION;
}

void yglUpdateProperties(void)
{
  int chg_sides;

  if(glCurrWin3d->curr_poly_sides != glCurrWin3d->poly_sides) {
    chg_sides= 1;
    glCurrWin3d->curr_poly_sides= glCurrWin3d->poly_sides;
  } else {
    chg_sides= 0;
  }
  if(chg_sides || glCurrWin3d->curr_poly_mode != glCurrWin3d->poly_mode) {
    glCurrWin3d->curr_poly_mode= glCurrWin3d->poly_mode;
    glPolygonMode(glCurrWin3d->curr_poly_sides, glCurrWin3d->curr_poly_mode);
  }
  if(chg_sides || glCurrWin3d->curr_mat_spec[0] != glCurrWin3d->mat_spec) {
    glCurrWin3d->curr_mat_spec[0]= glCurrWin3d->mat_spec;
    glCurrWin3d->curr_mat_spec[1]= glCurrWin3d->mat_spec;
    glCurrWin3d->curr_mat_spec[2]= glCurrWin3d->mat_spec;
    glMaterialfv(glCurrWin3d->curr_poly_sides, GL_SPECULAR, glCurrWin3d->curr_mat_spec);
  }
  if(glCurrWin3d->curr_cull_mode != glCurrWin3d->cull_mode) {
    glCurrWin3d->curr_cull_mode= glCurrWin3d->cull_mode;
    if(glCurrWin3d->curr_cull_mode) {
      glEnable(GL_CULL_FACE);
    } else {
      glDisable(GL_CULL_FACE);
    }
  }
  if(chg_sides || glCurrWin3d->curr_mat_color != glCurrWin3d->mat_color) {
    glCurrWin3d->curr_mat_color= glCurrWin3d->mat_color;
    glColorMaterial(glCurrWin3d->curr_poly_sides, glCurrWin3d->curr_mat_color);
    glEnable(GL_COLOR_MATERIAL);
  }
  if(chg_sides) {
    glMateriali(glCurrWin3d->curr_poly_sides, GL_SHININESS, 100);
  }
  if(glCurrWin3d->curr_shade_model != glCurrWin3d->shade_model) {
    glCurrWin3d->curr_shade_model= glCurrWin3d->shade_model;
    glShadeModel(glCurrWin3d->curr_shade_model);
  }
}

void yglForceUpdateProperties(void)
{
  glPolygonMode(glCurrWin3d->curr_poly_sides, glCurrWin3d->curr_poly_mode);
  glMaterialfv(glCurrWin3d->curr_poly_sides, GL_SPECULAR, glCurrWin3d->curr_mat_spec);
  if(glCurrWin3d->curr_cull_mode) {
    glEnable(GL_CULL_FACE);
  } else {
    glDisable(GL_CULL_FACE);
  }
  glMateriali(glCurrWin3d->curr_poly_sides, GL_SHININESS, 100);
  glColorMaterial(glCurrWin3d->curr_poly_sides, glCurrWin3d->curr_mat_color);
  /*  glColorMaterial(curr_poly_sides, GL_AMBIENT_AND_DIFFUSE); */
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(glCurrWin3d->curr_shade_model);
}

void yglWireSphere(int list_num, double radius)
{
  /* draw a wire frame sphere with radius 1.0 centered at
     the origin. Draw lines every 30 degrees in latitude and
     longitude. Put a marker at the north pole.
     longitude lines in the western hemisphere are green
     and in the eastern hemisphere cyan.
	*/
  int i, j, numi, numj, numli, numlj, num2i;
  double lat, lon, x, y, z, csth, snth, csph, snph, pi;

  if(alpha_pass) return;
#define ang_step 30
#define lin_step 5
  pi= 4.0*atan(1.0);
  numi= (int)(360.0/ang_step);
  numj= (int)(180.0/ang_step);
  numli= (int)(360.0/lin_step);
  numlj= (int)(180.0/lin_step);
  glNewList(list_num, GL_COMPILE);
  CHEK_ERROR("yglWireSphere glNewList");
  /* draw latitude lines */
  for(j= 1; j < numj; j++) {
    lat= j*pi/numj;
    csth= cos(lat);
    snth= sqrt(1.0-csth*csth);
    glBegin(GL_LINE_STRIP);
    glColor3d(1.0,1.0,1.0);
    for(i= 0; i <= numli; i++) {
      lon= i*2.0*pi/numli;
      csph= cos(lon);
      snph= sin(lon);
      x= radius*csph*snth;
      y= radius*snph*snth;
      z= radius*csth;
      glVertex3d(x,y,z);
    }
    glEnd();
  }
  CHEK_ERROR("yglWireSphere after latitude");
  /* draw longitude lines */
  num2i= numi/2;
  for(i= 0; i < numi; i++) {
    lon= i*2.0*pi/numi;
    csph= cos(lon);
    snph= sin(lon);
    glBegin(GL_LINE_STRIP);
    if(i < num2i) glColor3d(0.0, 0.0, 0.0);
    else glColor3d(0.0, 1.0, 0.0);
    for(j= 0; j <= numlj; j++) {
      lat= j*pi/numlj;
      csth= cos(lat);
      snth= sqrt(1.0-csth*csth);
      x= radius*csph*snth;
      y= radius*snph*snth;
      z= radius*csth;
      glVertex3d(x,y,z);
    }
    glEnd();
  }
  CHEK_ERROR("yglWireSphere after longitude");
  /* draw polar cap */
  glBegin(GL_POLYGON);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glColor3d(1.0,1.0,1.0);
  for(i= 0; i <= numi; i++) {
    lon= i*2.0*pi/numi;
    csph= cos(lon);
    snph= sin(lon);
    csth= cos(pi/18.0);
    snth= sqrt(1.0-csth*csth);
    x= radius*csph*snth;
    y= radius*snph*snth;
    z= radius*csth;
    glVertex3d(x,y,z);
  }
  glEnd();
  CHEK_ERROR("yglWireSphere after polar cap");
  glEndList();
  CHEK_ERROR("yglWireSphere");
}

int yglNewList(void)
{
  int num;
  /* create one new OpenGL display list */
  num= glGenLists(1);
  glNewList(num, GL_COMPILE);
  return num;
}

int yglCloseList(void)
{
  /* closes the current OpenGL display list (only one 
     can be open at a time) */
  glEndList();
  return 0;
}

int yglDrawList(int num)
{
  /* draw the specified OpenGL display list */
  glCallList(num);
  return 0;
}

int yglDeleteList(int num)
{
  /* free the specified OpenGL display list */
  /* !!! WARNING  !! FIX FIX FIX !!
     OpenGL display lists are destroyed when their OpenGL 
     context is destroyed. Yorick must thus manage the
     list of active display lists. This should be per
     window information */
  return 0;
}

int yglArrlim3d(long nvert, double *arr, double *lims)
{
  long base;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double xval, yval, zval;

  xmin= 1.0e100;
  xmax= -1.0e100;
  ymin= 1.0e100;
  ymax= -1.0e100;
  zmin= 1.0e100;
  zmax= -1.0e100;
  /* the x-array is 3-by-nvert */
  for(base= 0; base < nvert; base += 3) {
    xval= arr[base];
    yval= arr[base+1];
    zval= arr[base+2];
    if(xmin > xval) xmin= xval;
    if(xmax < xval) xmax= xval;
    if(ymin > yval) ymin= yval;
    if(ymax < yval) ymax= yval;
    if(zmin > zval) zmin= zval;
    if(zmax < zval) zmax= zval;
  }
  lims[0]= xmin;
  lims[1]= xmax;
  lims[2]= ymin;
  lims[3]= ymax;
  lims[4]= zmin;
  lims[5]= zmax;
  return 0;
}
