/*
 * $Id: glPolys.c,v 1.3 2006-10-19 14:48:19 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "glcode.h"
#include "glfunc.h"
#include "glPolys.h"
#include "pstdlib.h"
#include <math.h>
#include <stdlib.h>

static void yglPolysArrNoLite(long npoly, long *nverts, float *xyz, 
                              float *colr);
static void yglPolysNoArrNoLite(long npoly, long *nverts, float *xyz, 
                                float *colr);
static void yglPolysSmArr(long npoly, long *nverts, float *xyz, float *norm, 
			  float *colr);
static void yglPolysSmNoArr(long npoly, long *nverts, float *xyz, float *norm, 
                            float *colr);
static void yglPolysArr(long npoly, long *nverts, float *xyz, float *norm, 
                        float *colr);
static void yglPolysNoArr(long npoly, long *nverts, float *xyz, float *norm, 
			  float *colr);

/*
#define CHEK_ERROR(x)	gl_chek_error(x)
*/
#define CHEK_ERROR(x)

void yglPolys(long npolys, long *len, float *xyz, float *norm, 
              float *colr, long edge, long smooth, long do_light)
{
  /* The input is a list of polygons. There is one color per polygon.
     There is one normal per vertex or one normal per polygon
     (depending on smooth), or none for no lighting.
  */
  float oldSpec;

  if(alpha_pass) return;
  yglSetPolyMode(edge);
  oldSpec= yglGetMatSpec();
  if(!do_light) {
    yglSetMatSpec(0.0);  /* turn off specular highlights */
    /* polys are not lit */
    yglSetColorType(0);
    if(glCurrWin3d->use_array) {
      yglPolysArrNoLite(npolys, len, xyz, colr);
    } else {
      yglPolysNoArrNoLite(npolys, len, xyz, colr);
    }
  } else if(smooth) {
    /* smooth shading */
    yglSetColorType(1);
    if(glCurrWin3d->use_array) {
      yglPolysSmArr(npolys, len, xyz, norm, colr);
    } else {
      yglPolysSmNoArr(npolys, len, xyz, norm, colr);
    }
    yglSetMatSpec(oldSpec);
  } else {
    /* flat shading */
    yglSetColorType(1);
    if(glCurrWin3d->use_array) {
      yglPolysArr(npolys, len, xyz, norm, colr);
    } else {
      yglPolysNoArr(npolys, len, xyz, norm, colr);
    }
  }
}

void yglPolysArr(long npoly, long *nverts, float *xyz, float *norm, 
                 float *colr)
{
  /* This function expands each polygon out into a series of
     triangles and will INCREASE the amount of data sent
     to the server. */
  long i, j;
  float *arr_colr, *arr_norm, *arr_vert;
  long ntri, ind, vrt;
#ifdef USE_ALPHA
  long indc;
#endif

  /* draw the polygon list now */
  if(npoly <= 0) return;
#ifdef USE_ALPHA
  if(!alpha_pass) return;
#else
  if(alpha_pass) return;
#endif

#ifdef USE_ALPHA
  glEnable(GL_BLEND);
#endif
  /* use flat shading (no per-vertex normals supplied) */
  yglSetShade(0);
  yglUpdateProperties();

  /* NOTE: cannot use an index list, because the indices must be the 
     same for all active arrays. This means that I could only do
     one polygon at a time, and that does not seem worthwhile.
  */
  /* compute the number of triangles in the set */
  for(i= 0, ntri= 0; i < npoly; i++) {
    ntri += nverts[i]-2;
  }
  /* make an array big enough to hold all the colors, normals, and
     vertices for the polygon list */
#ifdef USE_ALPHA
  arr_colr= (float *)p_malloc(12*ntri*sizeof(float));
#else
  arr_colr= (float *)p_malloc(9*ntri*sizeof(float));
#endif
  arr_norm= (float *)p_malloc(9*ntri*sizeof(float));
  arr_vert= (float *)p_malloc(9*ntri*sizeof(float));
  ind= 0;
#ifdef USE_ALPHA
  indc= 0;
#endif
  vrt= 0;
  for(i= 0; i < npoly; i++) {
    long nv= nverts[i];
    long i3= 3*i;
    /* loop through all triangles for this polygon */
    for(j= 1; j < nv-1; j++) {
#ifdef USE_ALPHA
      arr_colr[indc]= colr[i3];
      arr_colr[indc+1]= colr[i3+1];
      arr_colr[indc+2]= colr[i3+2];
      arr_colr[indc+3]= curr_alpha;
      arr_colr[indc+4]= colr[i3];
      arr_colr[indc+5]= colr[i3+1];
      arr_colr[indc+6]= colr[i3+2];
      arr_colr[indc+7]= curr_alpha;
      arr_colr[indc+8]= colr[i3];
      arr_colr[indc+9]= colr[i3+1];
      arr_colr[indc+10]= colr[i3+2];
      arr_colr[indc+11]= curr_alpha;
      indc += 12;
#else
      arr_colr[ind]= colr[i3];
      arr_colr[ind+1]= colr[i3+1];
      arr_colr[ind+2]= colr[i3+2];
      arr_colr[ind+3]= colr[i3];
      arr_colr[ind+4]= colr[i3+1];
      arr_colr[ind+5]= colr[i3+2];
      arr_colr[ind+6]= colr[i3];
      arr_colr[ind+7]= colr[i3+1];
      arr_colr[ind+8]= colr[i3+2];
#endif

      arr_norm[ind]= norm[i3];
      arr_norm[ind+1]= norm[i3+1];
      arr_norm[ind+2]= norm[i3+2];
      arr_norm[ind+3]= norm[i3];
      arr_norm[ind+4]= norm[i3+1];
      arr_norm[ind+5]= norm[i3+2];
      arr_norm[ind+6]= norm[i3];
      arr_norm[ind+7]= norm[i3+1];
      arr_norm[ind+8]= norm[i3+2];

      arr_vert[ind]= xyz[vrt];
      arr_vert[ind+1]= xyz[vrt+1];
      arr_vert[ind+2]= xyz[vrt+2];
      arr_vert[ind+3]= xyz[vrt+3*j];
      arr_vert[ind+4]= xyz[vrt+3*j+1];
      arr_vert[ind+5]= xyz[vrt+3*j+2];
      arr_vert[ind+6]= xyz[vrt+3*j+3];
      arr_vert[ind+7]= xyz[vrt+3*j+4];
      arr_vert[ind+8]= xyz[vrt+3*j+5];

      ind += 9;
    }
    vrt += 3*nv;
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
#ifdef USE_ALPHA
  glColorPointer(4, GL_FLOAT, 0, arr_colr);
#else
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
#endif
  glNormalPointer(GL_FLOAT, 0, arr_norm);
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles using previously specified colors, vertices,
     and normals */
  glDrawArrays(GL_TRIANGLES, 0, 3*ntri);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_norm);
  p_free(arr_vert);
#ifdef USE_ALPHA
  glDisable(GL_BLEND);
#endif
  CHEK_ERROR("gl_polys_arr");
}

void yglPolysArrNoLite(long npoly, long *nverts, float *xyz, 
                       float *colr)
{
  /* This function expands each polygon out into a series of
     triangles and will INCREASE the amount of data sent
     to the server. */
  long i, j;
  float *arr_colr, *arr_vert;
  long ntri, ind, vrt;
#ifdef USE_ALPHA
  long indc;
#endif

  /* draw the polygon list now */
  if(npoly <= 0) return;
#ifdef USE_ALPHA
  if(!alpha_pass) return;
#else
  if(alpha_pass) return;
#endif

#ifdef USE_ALPHA
  glEnable(GL_BLEND);
#endif
  /* use flat shading (no per-vertex normals supplied) */
  yglSetShade(0);
  yglUpdateProperties();

  /* NOTE: cannot use an index list, because the indices must be the 
     same for all active arrays. This means that I could only do
     one polygon at a time, and that does not seem worthwhile.
  */
  /* compute the number of triangles in the set */
  for(i= 0, ntri= 0; i < npoly; i++) {
    ntri += nverts[i]-2;
  }
  /* make an array big enough to hold all the colors, normals, and
     vertices for the polygon list */
#ifdef USE_ALPHA
  arr_colr= (float *)p_malloc(12*ntri*sizeof(float));
#else
  arr_colr= (float *)p_malloc(9*ntri*sizeof(float));
#endif
  arr_vert= (float *)p_malloc(9*ntri*sizeof(float));
  ind= 0;
#ifdef USE_ALPHA
  indc= 0;
#endif
  vrt= 0;
  for(i= 0; i < npoly; i++) {
    long nv= nverts[i];
    long i3= 3*i;
    /* loop through all triangles for this polygon */
    for(j= 1; j < nv-1; j++) {
#ifdef USE_ALPHA
      arr_colr[indc]= colr[i3];
      arr_colr[indc+1]= colr[i3+1];
      arr_colr[indc+2]= colr[i3+2];
      arr_colr[indc+3]= curr_alpha;
      arr_colr[indc+4]= colr[i3];
      arr_colr[indc+5]= colr[i3+1];
      arr_colr[indc+6]= colr[i3+2];
      arr_colr[indc+7]= curr_alpha;
      arr_colr[indc+8]= colr[i3];
      arr_colr[indc+9]= colr[i3+1];
      arr_colr[indc+10]= colr[i3+2];
      arr_colr[indc+11]= curr_alpha;
      indc += 12;
#else
      arr_colr[ind]= colr[i3];
      arr_colr[ind+1]= colr[i3+1];
      arr_colr[ind+2]= colr[i3+2];
      arr_colr[ind+3]= colr[i3];
      arr_colr[ind+4]= colr[i3+1];
      arr_colr[ind+5]= colr[i3+2];
      arr_colr[ind+6]= colr[i3];
      arr_colr[ind+7]= colr[i3+1];
      arr_colr[ind+8]= colr[i3+2];
#endif

      arr_vert[ind]= xyz[vrt];
      arr_vert[ind+1]= xyz[vrt+1];
      arr_vert[ind+2]= xyz[vrt+2];
      arr_vert[ind+3]= xyz[vrt+3*j];
      arr_vert[ind+4]= xyz[vrt+3*j+1];
      arr_vert[ind+5]= xyz[vrt+3*j+2];
      arr_vert[ind+6]= xyz[vrt+3*j+3];
      arr_vert[ind+7]= xyz[vrt+3*j+4];
      arr_vert[ind+8]= xyz[vrt+3*j+5];

      ind += 9;
    }
    vrt += 3*nv;
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
#ifdef USE_ALPHA
  glColorPointer(4, GL_FLOAT, 0, arr_colr);
#else
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
#endif
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles using previously specified colors, vertices,
     and normals */
  glDrawArrays(GL_TRIANGLES, 0, 3*ntri);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_vert);
#ifdef USE_ALPHA
  glDisable(GL_BLEND);
#endif
  CHEK_ERROR("gl_polys_arr");
}

void yglPolysNoArr(long npoly, long *nverts, float *xyz, float *norm, 
                   float *colr)
{
  long i, j, base;

  /* draw the polygon list now */
  if(npoly <= 0) return;
#ifdef USE_ALPHA
  if(!alpha_pass) return;
#else
  if(alpha_pass) return;
#endif

#ifdef USE_ALPHA
  glEnable(GL_BLEND);
#endif
  /* use flat shading (no per-vertex normals supplied) */
  yglSetShade(0);
  yglUpdateProperties();

  base= 0;
  for(i= 0; i < npoly; i++) {
#ifdef USE_TRI
    glBegin(GL_TRIANGLE_FAN);
#else
    glBegin(GL_POLYGON);
#endif
#ifdef USE_ALPHA
    glColor4f(*(colr+3*i), *(colr+3*i+1), *(colr+3*i+2), curr_alpha);
#else
    glColor3fv(colr+3*i);
#endif
    glNormal3fv(norm+3*i);
    for(j= 0; j < nverts[i]; j++) {
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      glVertex3fv(xyz+base);
      base += 3;
    }
    glEnd();
  }
#ifdef USE_ALPHA
  glDisable(GL_BLEND);
#endif
  CHEK_ERROR("gl_polys_no_arr");
}

void yglPolysNoArrNoLite(long npoly, long *nverts, float *xyz, 
                         float *colr)
{
  long i, j, base;

  /* draw the polygon list now */
  if(npoly <= 0) return;
#ifdef USE_ALPHA
  if(!alpha_pass) return;
#else
  if(alpha_pass) return;
#endif

#ifdef USE_ALPHA
  glEnable(GL_BLEND);
#endif
  /* use flat shading */
  yglSetShade(0);
  yglUpdateProperties();

  base= 0;
  for(i= 0; i < npoly; i++) {
#ifdef USE_TRI
    glBegin(GL_TRIANGLE_FAN);
#else
    glBegin(GL_POLYGON);
#endif
#ifdef USE_ALPHA
    glColor4f(*(colr+3*i), *(colr+3*i+1), *(colr+3*i+2), curr_alpha);
#else
    glColor3fv(colr+3*i);
#endif
    for(j= 0; j < nverts[i]; j++) {
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      glVertex3fv(xyz+base);
      base += 3;
    }
    glEnd();
  }
#ifdef USE_ALPHA
  glDisable(GL_BLEND);
#endif
  CHEK_ERROR("gl_polys_no_arr");
}

void yglPolysSmArr(long npoly, long *nverts, float *xyz, float *norm, 
                   float *colr)
{
  /* This function expands each polygon out into a series of
     triangles and will INCREASE the amount of data sent
     to the server. */
  long i, j;
  float *arr_colr, *arr_norm, *arr_vert;
  long ntri, ind, vrt;
#ifdef USE_ALPHA
  long indc;
#endif

  /* draw the polygon list, with a normal per vertex */
  if(npoly <= 0) return;
#ifdef USE_ALPHA
  if(!alpha_pass) return;
#else
  if(alpha_pass) return;
#endif

#ifdef USE_ALPHA
  glEnable(GL_BLEND);
#endif

  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();

  /* NOTE: cannot use an index list, because the indices must be the 
     same for all active arrays. This means that I could only do
     one polygon at a time, and that does not seem worthwhile.
  */
  /* compute the number of triangles in the set */
  for(i= 0, ntri= 0; i < npoly; i++) {
    ntri += nverts[i]-2;
  }
  /* make an array big enough to hold all the colors, normals, and
     vertices for the polygon list */
#ifdef USE_ALPHA
  arr_colr= (float *)p_malloc(12*ntri*sizeof(float));
#else
  arr_colr= (float *)p_malloc(9*ntri*sizeof(float));
#endif
  arr_norm= (float *)p_malloc(9*ntri*sizeof(float));
  arr_vert= (float *)p_malloc(9*ntri*sizeof(float));
  ind= 0;
#ifdef USE_ALPHA
  indc= 0;
#endif
  vrt= 0;
  for(i= 0; i < npoly; i++) {
    long nv= nverts[i];
    long i3= 3*i;
    /* loop through all triangles for this polygon */
    for(j= 1; j < nv-1; j++) {
#ifdef USE_ALPHA
      arr_colr[indc]= colr[i3];
      arr_colr[indc+1]= colr[i3+1];
      arr_colr[indc+2]= colr[i3+2];
      arr_colr[indc+3]= curr_alpha;
      arr_colr[indc+4]= colr[i3];
      arr_colr[indc+5]= colr[i3+1];
      arr_colr[indc+6]= colr[i3+2];
      arr_colr[indc+7]= curr_alpha;
      arr_colr[indc+8]= colr[i3];
      arr_colr[indc+9]= colr[i3+1];
      arr_colr[indc+10]= colr[i3+2];
      arr_colr[indc+11]= curr_alpha;
      indc += 12;
#else
      arr_colr[ind]= colr[i3];
      arr_colr[ind+1]= colr[i3+1];
      arr_colr[ind+2]= colr[i3+2];
      arr_colr[ind+3]= colr[i3];
      arr_colr[ind+4]= colr[i3+1];
      arr_colr[ind+5]= colr[i3+2];
      arr_colr[ind+6]= colr[i3];
      arr_colr[ind+7]= colr[i3+1];
      arr_colr[ind+8]= colr[i3+2];
#endif

      arr_norm[ind]= norm[vrt];
      arr_norm[ind+1]= norm[vrt+1];
      arr_norm[ind+2]= norm[vrt+2];
      arr_norm[ind+3]= norm[vrt+3*j];
      arr_norm[ind+4]= norm[vrt+3*j+1];
      arr_norm[ind+5]= norm[vrt+3*j+2];
      arr_norm[ind+6]= norm[vrt+3*j+3];
      arr_norm[ind+7]= norm[vrt+3*j+4];
      arr_norm[ind+8]= norm[vrt+3*j+5];

      arr_vert[ind]= xyz[vrt];
      arr_vert[ind+1]= xyz[vrt+1];
      arr_vert[ind+2]= xyz[vrt+2];
      arr_vert[ind+3]= xyz[vrt+3*j];
      arr_vert[ind+4]= xyz[vrt+3*j+1];
      arr_vert[ind+5]= xyz[vrt+3*j+2];
      arr_vert[ind+6]= xyz[vrt+3*j+3];
      arr_vert[ind+7]= xyz[vrt+3*j+4];
      arr_vert[ind+8]= xyz[vrt+3*j+5];

      ind += 9;
    }
    vrt += 3*nv;
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
#ifdef USE_ALPHA
  glColorPointer(4, GL_FLOAT, 0, arr_colr);
#else
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
#endif
  glNormalPointer(GL_FLOAT, 0, arr_norm);
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles using previously specified colors, vertices,
     and normals */
  glDrawArrays(GL_TRIANGLES, 0, 3*ntri);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_norm);
  p_free(arr_vert);
#ifdef USE_ALPHA
  glDisable(GL_BLEND);
#endif
  CHEK_ERROR("gl_polys_sm_arr");
}

void yglPolysSmNoArr(long npoly, long *nverts, float *xyz, float *norm, 
                     float *colr)
{
  long i, j, base;

  /* draw the polygon list, with a normal per vertex */
  if(npoly <= 0) return;
#ifdef USE_ALPHA
  if(!alpha_pass) return;
#else
  if(alpha_pass) return;
#endif

#ifdef USE_ALPHA
  glEnable(GL_BLEND);
#endif

  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();

  base= 0;
  for(i= 0; i < npoly; i++) {
#ifdef USE_TRI
    glBegin(GL_TRIANGLE_FAN);
#else
    glBegin(GL_POLYGON);
#endif
#ifdef USE_ALPHA
    glColor4f(*(colr+3*i), *(colr+3*i+1), *(colr+3*i+2), curr_alpha);
#else
    glColor3fv(colr+3*i);
#endif
    for(j= 0; j < nverts[i]; j++) {
      glNormal3fv(norm+base);
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      glVertex3fv(xyz+base);
      base += 3;
    }
    glEnd();
  }
#ifdef USE_ALPHA
  glDisable(GL_BLEND);
#endif
  CHEK_ERROR("gl_polys_sm_no_arr");
}
