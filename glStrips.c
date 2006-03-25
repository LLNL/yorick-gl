/*
 * $Id: glStrips.c,v 1.2 2006-03-25 03:12:29 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "glcode.h"
#include "glfunc.h"
#include "glStrips.h"
#include "glWrappers.h"
#include "pstdlib.h"
#include <math.h>
#include <stdlib.h>

static void yglTstripsSmNoArrAlpha(long nstrip, long *len, float *xyz, 
                                   float *norm, float *colr);
static void yglTstripSmNoArr(long nvert, float *xyz, float *norm, float *colr);
static void yglQstripsSmNoArrAlpha(long nstrip, long *len, float *xyz, 
                                   float *norm, float *colr);
static void yglQstripSmNoArr(long nvert, float *xyz, float *norm, float *colr);
static void yglTstripNoArrNoLiteAlpha(long nvert, float *xyz, float *colr);
static void yglTstripArrNoLite(long nvert, float *xyz, float *colr);
static void yglTstripNoArrNoLite(long nvert, float *xyz, float *colr);
static void yglTstripSmArr(long nvert, float *xyz, float *norm, float *colr);
static void yglTstripsSmNoArr(long nstrip, long *len, float *xyz, 
                              float *norm, float *colr);
static void yglTstripArr(long nvert, float *xyz, float *norm, float *colr);
static void yglTstripNoArr(long nvert, float *xyz, float *norm, float *colr);
static void yglQstripArrNoLite(long nvert, float *xyz, float *colr);
static void yglQstripSmArr(long nvert, float *xyz, float *norm, float *colr);
static void yglQstripsSmNoArr(long nstrip, long *len, float *xyz, float *norm,
                              float *colr);
static void yglQstripArr(long nvert, float *xyz, float *norm, float *colr);
static void yglQstripNoArr(long nvert, float *xyz, float *norm, float *colr);
static void yglQstripNoArrNoLite(long nvert, float *xyz, float *colr);

/*
#define CHEK_ERROR(x)	yygl_chek_error(x)
*/
#define CHEK_ERROR(x)

void yglTstrips(long nstrip, long *len, float *xyz, float *norm, 
                float *colr, long edge, long smooth, long do_light)
{
  /* The input is a list of triangle strips. There is one color per triangle.
     There is one normal per vertex or one normal per triangle
     (depending on smooth), or none for no lighting.
  */
  long i, num, vert, ntri;
  float oldSpec;

  yglSetPolyMode(edge);
  if(smooth) {
    /* use smooth shading */
    yglSetShade(1);
  } else {
    /* use flat shading */
    yglSetShade(0);
  }
  yglUpdateProperties();
  if(!do_light) {
    oldSpec= yglGetMatSpec();
    yglSetMatSpec(0.0);  /* turn off specular highlights */
    /* triangles are not lit */
    yglSetColorType(0);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripArrNoLite(num, xyz+vert*3, colr+3*ntri);
        vert += num;
	ntri += num-2;
      }
    } else {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripNoArrNoLite(num, xyz+vert*3, colr+3*ntri);
        vert += num;
        ntri += num-2;
      }
    }
    yglSetMatSpec(oldSpec);
  } else if(smooth) {
    /* smooth shading */
    yglSetColorType(1);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripSmArr(num, xyz+vert*3, norm+vert*3, colr+3*ntri);
        vert += num;
        ntri += num-2;
      }
    } else {
      yglTstripsSmNoArr(nstrip, len, xyz, norm, colr);
    }
  } else {
    /* flat shading */
    yglSetColorType(1);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripArr(num, xyz+vert*3, norm+ntri*3, colr+3*ntri);
        vert += num;
        ntri += num-2;
      }
    } else {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripNoArr(num, xyz+vert*3, norm+ntri*3, colr+3*ntri);
        vert += num;
        ntri += num-2;
      }
    }
  }
}

void yglTstripsAlpha(long nstrip, long *len, float *xyz, float *norm, 
                     float *colr, long edge, long smooth, long do_light)
{
  /* The input is a list of triangle strips. There is one color per triangle.
     There is one normal per vertex or one normal per triangle
     (depending on smooth), or none for no lighting.
  */
  long i, num, vert, ntri;
  float oldSpec;
#if 0
  GLfloat mat_white[]= {1.0, 1.0, 1.0, 1.0};
  GLfloat mat_transparent[]= {0.0f, 0.0f, 0.0f, 0.3f};
  GLfloat mat_emission[]= {0.0f, 1.0f, 0.0f, 0.3f};
  GLfloat mat_nul[]= {0.0, 0.0, 0.0, 1.0};
#endif
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);
  yglSetPolyMode(edge);
  if(smooth) {
    /* use smooth shading */
    yglSetShade(1);
  } else {
    /* use flat shading */
    yglSetShade(0);
  }

  if(!do_light) {
    oldSpec= yglGetMatSpec();
    yglSetMatSpec(0.0);  /* turn off specular highlights */
    /* triangles are not lit */
    yglSetColorType(0);
    yglUpdateProperties();
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripArrNoLite(num, xyz+vert*3, colr+4*ntri);
        vert += num;
        ntri += num-2;
      }
    } else {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripNoArrNoLiteAlpha(num, xyz+vert*3, colr+4*ntri);
        vert += num;
        ntri += num-2;
      }
    }
    yglSetMatSpec(oldSpec);
  } else if(smooth) {
    /* smooth shading */
    yglSetColorType(1);  /* set ambient and diffuse color */
    yglUpdateProperties();
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripSmArr(num, xyz+vert*3, norm+vert*3, colr+4*ntri);
        vert += num;
        ntri += num-2;
      }
    } else {
      yglTstripsSmNoArrAlpha(nstrip, len, xyz, norm, colr);
    }
  } else {
    /* flat shading */
    yglSetColorType(1);
    yglUpdateProperties();
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripArr(num, xyz+vert*3, norm+ntri*3, colr+4*ntri);
        vert += num;
        ntri += num-2;
      }
    } else {
      for(i= 0, vert= 0, ntri= 0; i < nstrip; i++) {
        num= len[i];
        yglTstripNoArr(num, xyz+vert*3, norm+ntri*3, colr+4*ntri);
        vert += num;
        ntri += num-2;
      }
    }
  }
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
}

void yglQstrips(long nstrip, long *len, float *xyz, float *norm, 
                float *colr, long edge, long smooth, long do_light)
{
  /* The input is a list of quad strips. There is one color per quad.
     There is one normal per vertex or one normal per quad
     (depending on smooth).
  */
  long i, num, vert, nquad;
  float oldSpec;

  yglSetPolyMode(edge);
  if(smooth) {
    /* use smooth shading */
    yglSetShade(1);
  } else {
    /* use flat shading */
    yglSetShade(0);
  }
  yglUpdateProperties();
  if(!do_light) {
    oldSpec= yglGetMatSpec();
    yglSetMatSpec(0.0);  /* turn off specular highlights */
    /* quads are not lit */
    yglSetColorType(0);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripArrNoLite(num, xyz+vert*3, colr+3*nquad);
        vert += num;
        nquad += num-1;
      }
    } else {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripNoArrNoLite(num, xyz+vert*3, colr+3*nquad);
        vert += num;
        nquad += num-1;
      }
    }
    yglSetMatSpec(oldSpec);
  } else if(smooth) {
    /* smooth shading */
    yglSetColorType(1);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripSmArr(num, xyz+vert*3, norm+vert*3, colr+3*nquad);
        vert += num;
        nquad += num-1;
      }
    } else {
      yglQstripsSmNoArr(nstrip, len, xyz, norm, colr);
    }
  } else {
    /* flat shading */
    yglSetColorType(1);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripArr(num, xyz+vert*3, norm+nquad*3, colr+3*nquad);
        vert += num;
        nquad += num-1;
      }
    } else {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripNoArr(num, xyz+vert*3, norm+nquad*3, colr+3*nquad);
        vert += num;
        nquad += num-1;
      }
    }
  }
}

void yglQstripsAlpha(long nstrip, long *len, float *xyz, float *norm, 
                     float *colr, long edge, long smooth, long do_light)
{
  /* The input is a list of quad strips. There is one color per quad.
     There is one normal per vertex or one normal per quad
     (depending on smooth).
  */
  long i, num, vert, nquad;
  float oldSpec;

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  yglSetPolyMode(edge);
  if(smooth) {
    /* use smooth shading */
    yglSetShade(1);
  } else {
    /* use flat shading */
    yglSetShade(0);
  }
  yglUpdateProperties();
  if(!do_light) {
    oldSpec= yglGetMatSpec();
    yglSetMatSpec(0.0);  /* turn off specular highlights */
    /* quads are not lit */
    yglSetColorType(0);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripArrNoLite(num, xyz+vert*3, colr+4*nquad);
        vert += num;
        nquad += num-1;
      }
    } else {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripNoArrNoLite(num, xyz+vert*3, colr+4*nquad);
        vert += num;
        nquad += num-1;
      }
    }
    yglSetMatSpec(oldSpec);
  } else if(smooth) {
    /* smooth shading */
    yglSetColorType(1);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripSmArr(num, xyz+vert*3, norm+vert*3, colr+4*nquad);
        vert += num;
        nquad += num-1;
      }
    } else {
      yglQstripsSmNoArrAlpha(nstrip, len, xyz, norm, colr);
    }
  } else {
    /* flat shading */
    yglSetColorType(1);
    if(glCurrWin3d->use_array) {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripArr(num, xyz+vert*3, norm+nquad*3, colr+4*nquad);
        vert += num;
        nquad += num-1;
      }
    } else {
      for(i= 0, vert= 0, nquad= 0; i < nstrip; i++) {
        num= len[i];
        yglQstripNoArr(num, xyz+vert*3, norm+nquad*3, colr+4*nquad);
        vert += num;
        nquad += num-1;
      }
    }
  }
  glDisable(GL_BLEND);
}

void yglTstripsSmNoArr(long nstrip, long *len, float *xyz, 
                       float *norm, float *colr)
{
  long i, i3, num;
  float old_red= -1.0, old_green= -1.0, old_blue= -1.0;

  for(i= 0; i < nstrip; i++) {
    num= len[i];  /* number of vertices in this strip */
    /* draw the strip */
    DEMAND(num > 2, "triangle strip with less than 3 vertices in yglTstripsSmNoArr")
    glBegin(GL_TRIANGLE_STRIP);
    if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue) {
      old_red= colr[0];
      old_green= colr[1];
      old_blue= colr[2];
      glColor3fv(colr);
    }
    glNormal3fv(norm);
    glVertex3fv(xyz);
    glNormal3fv(norm+3);
    glVertex3fv(xyz+3);
    norm += 6;
    xyz += 6;
    for(i3= 0; i3 < num-2; i3++) {
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue) {
        old_red= colr[0];
        old_green= colr[1];
        old_blue= colr[2];
        glColor3fv(colr);
      }
      glNormal3fv(norm);
      glVertex3fv(xyz);
      colr += 3;
      xyz += 3;
      norm += 3;
    }
    glEnd();
  }
  CHEK_ERROR("yglTstripsSmNoArr");
}

void yglTstripsSmNoArrAlpha(long nstrip, long *len, float *xyz, 
                            float *norm, float *colr)
{
  long i, i3, num;
  float old_red= -1.0, old_green= -1.0, old_blue= -1.0, old_alpha= -1.0;

  for(i= 0; i < nstrip; i++) {
    num= len[i];  /* number of vertices in this strip */
    /* draw the strip */
    DEMAND(num > 2, "triangle strip with less than 3 vertices in yglTstripsSmNoArrAlpha")
    glBegin(GL_TRIANGLE_STRIP);
    if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue || 
       colr[3] != old_alpha) {
      old_red= colr[0];
      old_green= colr[1];
      old_blue= colr[2];
      old_alpha= colr[3];
      glColor4fv(colr);
    }
    glNormal3fv(norm);
    glVertex3fv(xyz);
    glNormal3fv(norm+3);
    glVertex3fv(xyz+3);
    norm += 6;
    xyz += 6;
    for(i3= 0; i3 < num-2; i3++) {
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue || 
         colr[3] != old_alpha) {
        old_red= colr[0];
        old_green= colr[1];
        old_blue= colr[2];
        old_alpha= colr[3];
        glColor4fv(colr);
      }
      glNormal3fv(norm);
      glVertex3fv(xyz);
      colr += 4;
      xyz += 3;
      norm += 3;
    }
    glEnd();
  }
  CHEK_ERROR("yglTstripsSmNoArrAlpha");
}

void yglTstripSmNoArr(long nvert, float *xyz, float *norm, float *colr)
{
  long i3;
  float old_red, old_green, old_blue;

  /* draw the strip */
  if(nvert <= 2) return;

  glBegin(GL_TRIANGLE_STRIP);
  old_red= colr[0];
  old_green= colr[1];
  old_blue= colr[2];
  glColor3fv(colr);
  glNormal3fv(norm);
  glVertex3fv(xyz);
  glNormal3fv(norm+3);
  glVertex3fv(xyz+3);
  for(i3= 0; i3 < 3*(nvert-2); i3 += 3) {
    /* NOTE: It makes no difference in local performance if the coords
       are converted to floats before being stored in the display list */
    if(colr[i3] != old_red || colr[i3+1] != old_green || 
       colr[i3+2] != old_blue) {
      old_red= colr[i3];
      old_green= colr[i3+1];
      old_blue= colr[i3+2];
      glColor3fv(colr+i3);
    }
    glNormal3fv(norm+i3+6);
    glVertex3fv(xyz+i3+6);
  }
  glEnd();
  CHEK_ERROR("yglTstripSmNoArr");
}

void yglQstripsSmNoArr(long nstrip, long *len, float *xyz, float *norm,
                       float *colr)
{
  long i, i3, num;
  float old_red= -1.0, old_green= -1.0, old_blue= -1.0;

  for(i= 0; i < nstrip; i++) {
    num= len[i];  /* number of vertices on one edge of this strip */
    DEMAND(num > 1, "quad strip with less than 2 vertices in yglQstripsSmNoArr");
    glBegin(GL_QUAD_STRIP);
    if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue) {
      old_red= colr[0];
      old_green= colr[1];
      old_blue= colr[2];
      glColor3fv(colr);
    }
    glNormal3fv(norm);
    glVertex3fv(xyz);
    glNormal3fv(norm+3);
    glVertex3fv(xyz+3);
    norm += 6;
    xyz += 6;
    for(i3= 0; i3 < num-1; i3++) {
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue) {
        old_red= colr[0];
        old_green= colr[1];
        old_blue= colr[2];
        glColor3fv(colr);
      }
      glNormal3fv(norm);
      glVertex3fv(xyz);
      glNormal3fv(norm+3);
      glVertex3fv(xyz+3);
      colr += 3;
      xyz += 6;
      norm += 6;
    }
    glEnd();
  }
  CHEK_ERROR("yglQstripsSmNoArr");
}

void yglQstripsSmNoArrAlpha(long nstrip, long *len, float *xyz, 
                            float *norm, float *colr)
{
  long i, i3, num;
  float old_red= -1.0, old_green= -1.0, old_blue= -1.0, old_alpha= -1.0;

  for(i= 0; i < nstrip; i++) {
    num= len[i];  /* number of vertices on one edge of this strip */
    DEMAND(num > 1,"quad strip with less than 2 vertices in yglQstripsSmNoArrAlpha")
    glBegin(GL_QUAD_STRIP);
    if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue || 
       colr[3] != old_alpha) {
      old_red= colr[0];
      old_green= colr[1];
      old_blue= colr[2];
      old_alpha= colr[3];
      glColor4fv(colr);
    }
    glNormal3fv(norm);
    glVertex3fv(xyz);
    glNormal3fv(norm+3);
    glVertex3fv(xyz+3);
    norm += 6;
    xyz += 6;
    for(i3= 0; i3 < num-1; i3++) {
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue || 
         colr[3] != old_alpha) {
        old_red= colr[0];
        old_green= colr[1];
        old_blue= colr[2];
        old_alpha= colr[3];
        glColor4fv(colr);
      }
      glNormal3fv(norm);
      glVertex3fv(xyz);
      glNormal3fv(norm+3);
      glVertex3fv(xyz+3);
      colr += 4;
      xyz += 6;
      norm += 6;
    }
    glEnd();
  }
  CHEK_ERROR("yglQstripsSmNoArrAlpha");
}

void yglQstripSmNoArr(long nvert, float *xyz, float *norm, float *colr)
{
  long i3, i6;

  /* draw the strip */
  if(nvert <= 1) return;

  glBegin(GL_QUAD_STRIP);
  glColor3fv(colr);
  glNormal3fv(norm);
  glVertex3fv(xyz);
  glNormal3fv(norm+3);
  glVertex3fv(xyz+3);
  for(i3= 0, i6= 6; i3 < 3*(nvert-1); i3 += 3, i6 += 6) {
    /* NOTE: It makes no difference in local performance if the coords
       are converted to floats before being stored in the display list */
    glColor3fv(colr+i3);
    glNormal3fv(norm+i6);
    glVertex3fv(xyz+i6);
    glNormal3fv(norm+i6+3);
    glVertex3fv(xyz+i6+3);
  }
  glEnd();
  CHEK_ERROR("yglQstripSmNoArr");
}

void yglTstripSmArr(long nvert, float *xyz, float *norm, float *colr)
{
  long i3;
  float *arr_colr, *arr_norm, *arr_vert;
  long ind_colr, ind_vert;

  /* draw the strip */
  if(nvert <= 2) return;

  /* make an array big enough to hold all the colors, normals, and
     vertices for the quad strip */
  ind_colr= ind_vert= 0;
  arr_colr= (float *)p_malloc(3*nvert*sizeof(float));
  arr_norm= (float *)p_malloc(3*nvert*sizeof(float));
  arr_vert= (float *)p_malloc(3*nvert*sizeof(float));
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  for(i3= 0; i3 < 3*(nvert-2); i3 += 3) {
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
  }
  for(i3= 0; i3 < 3*nvert; i3 += 3) {
    arr_norm[ind_vert]= norm[i3];
    arr_vert[ind_vert++]= xyz[i3];
    arr_norm[ind_vert]= norm[i3+1];
    arr_vert[ind_vert++]= xyz[i3+1];
    arr_norm[ind_vert]= norm[i3+2];
    arr_vert[ind_vert++]= xyz[i3+2];
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
  glNormalPointer(GL_FLOAT, 0, arr_norm);
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles/quads using previously specified colors, vertices,
     and normals */
  glDrawArrays(GL_TRIANGLE_STRIP, 0, nvert);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_norm);
  p_free(arr_vert);
  CHEK_ERROR("yglTstripSmArr");
}

void yglQstripSmArr(long nvert, float *xyz, float *norm, float *colr)
{
  long i3, i6;
  float *arr_colr, *arr_norm, *arr_vert;
  long ind_colr, ind_vert;

  /* draw the strip */
  if(nvert <= 1) return;

  /* make an array big enough to hold all the colors, normals, and
     vertices for the quad strip */
  ind_colr= ind_vert= 0;
  arr_colr= (float *)p_malloc(6*nvert*sizeof(float));
  arr_norm= (float *)p_malloc(6*nvert*sizeof(float));
  arr_vert= (float *)p_malloc(6*nvert*sizeof(float));
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  for(i3= 0; i3 < 3*(nvert-1); i3 += 3) {
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
  }
  for(i6= 0; i6 < 6*nvert; i6 += 6) {
    arr_norm[ind_vert]= norm[i6];
    arr_vert[ind_vert++]= xyz[i6];
    arr_norm[ind_vert]= norm[i6+1];
    arr_vert[ind_vert++]= xyz[i6+1];
    arr_norm[ind_vert]= norm[i6+2];
    arr_vert[ind_vert++]= xyz[i6+2];
    arr_norm[ind_vert]= norm[i6+3];
    arr_vert[ind_vert++]= xyz[i6+3];
    arr_norm[ind_vert]= norm[i6+4];
    arr_vert[ind_vert++]= xyz[i6+4];
    arr_norm[ind_vert]= norm[i6+5];
    arr_vert[ind_vert++]= xyz[i6+5];
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
  glNormalPointer(GL_FLOAT, 0, arr_norm);
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles/quads using previously specified colors, vertices,
     and normals */
  glDrawArrays(GL_QUAD_STRIP, 0, 2*nvert);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_norm);
  p_free(arr_vert);
  CHEK_ERROR("yglQstripSmArr");
}

void yglTstripArr(long nvert, float *xyz, float *norm, float *colr)
{
  long i3;
  float *arr_colr, *arr_norm, *arr_vert;
  long ind_colr, ind_norm, ind_vert;

  /* draw the strip */
  if(nvert <= 2) return;

  /* make an array big enough to hold all the colors, normals, and
     vertices for the quad strip */
  ind_colr= ind_norm= ind_vert= 0;
  arr_colr= (float *)p_malloc(3*nvert*sizeof(float));
  arr_norm= (float *)p_malloc(3*nvert*sizeof(float));
  arr_vert= (float *)p_malloc(3*nvert*sizeof(float));
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  arr_norm[ind_norm++]= norm[0];
  arr_norm[ind_norm++]= norm[1];
  arr_norm[ind_norm++]= norm[2];
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  arr_norm[ind_norm++]= norm[0];
  arr_norm[ind_norm++]= norm[1];
  arr_norm[ind_norm++]= norm[2];
  for(i3= 0; i3 < 3*(nvert-2); i3 += 3) {
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
    arr_norm[ind_norm++]= norm[i3];
    arr_norm[ind_norm++]= norm[i3+1];
    arr_norm[ind_norm++]= norm[i3+2];
  }
  for(i3= 0; i3 < 3*nvert; i3 += 3) {
    arr_vert[ind_vert++]= xyz[i3];
    arr_vert[ind_vert++]= xyz[i3+1];
    arr_vert[ind_vert++]= xyz[i3+2];
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
  glNormalPointer(GL_FLOAT, 0, arr_norm);
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles/quads using previously specified colors, vertices,
     and normals */
  glDrawArrays(GL_TRIANGLE_STRIP, 0, nvert);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_norm);
  p_free(arr_vert);
  CHEK_ERROR("yglTstripArr");
}

void yglQstripArr(long nvert, float *xyz, float *norm, float *colr)
{
  long i3, i6;
  float *arr_colr, *arr_norm, *arr_vert;
  long ind_colr, ind_norm, ind_vert;

  /* draw the strip */
  if(nvert <= 1) return;

  /* make an array big enough to hold all the colors, normals, and
     vertices for the quad strip */
  ind_colr= ind_norm= ind_vert= 0;
  arr_colr= (float *)p_malloc(6*nvert*sizeof(float));
  arr_norm= (float *)p_malloc(6*nvert*sizeof(float));
  arr_vert= (float *)p_malloc(6*nvert*sizeof(float));
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  arr_norm[ind_norm++]= norm[0];
  arr_norm[ind_norm++]= norm[1];
  arr_norm[ind_norm++]= norm[2];
  arr_norm[ind_norm++]= norm[0];
  arr_norm[ind_norm++]= norm[1];
  arr_norm[ind_norm++]= norm[2];
  for(i3= 0; i3 < 3*(nvert-1); i3 += 3) {
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
    arr_norm[ind_norm++]= norm[i3];
    arr_norm[ind_norm++]= norm[i3+1];
    arr_norm[ind_norm++]= norm[i3+2];
    arr_norm[ind_norm++]= norm[i3];
    arr_norm[ind_norm++]= norm[i3+1];
    arr_norm[ind_norm++]= norm[i3+2];
  }
  for(i6= 0; i6 < 6*nvert; i6 += 6) {
    arr_vert[ind_vert++]= xyz[i6];
    arr_vert[ind_vert++]= xyz[i6+1];
    arr_vert[ind_vert++]= xyz[i6+2];
    arr_vert[ind_vert++]= xyz[i6+3];
    arr_vert[ind_vert++]= xyz[i6+4];
    arr_vert[ind_vert++]= xyz[i6+5];
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
  glNormalPointer(GL_FLOAT, 0, arr_norm);
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles/quads using previously specified colors, vertices,
     and normals */
  glDrawArrays(GL_QUAD_STRIP, 0, 2*nvert);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_norm);
  p_free(arr_vert);
  CHEK_ERROR("yglQstripArr");
}

void yglTstripNoArr(long nvert, float *xyz, float *norm, float *colr)
{
  long i3;
  double old_red, old_green, old_blue;

  /* draw the strip */
  if(nvert <= 2) return;

  glBegin(GL_TRIANGLE_STRIP);
  glColor3fv(colr);
  old_red= colr[0];
  old_green= colr[1];
  old_blue= colr[2];
  glVertex3fv(xyz);
  glVertex3fv(xyz+3);
  for(i3= 0; i3 < 3*(nvert-2); i3 += 3) {
    if(colr[i3] != old_red || colr[i3+1] != old_green || 
       colr[i3+2] != old_blue) {
      glColor3fv(colr+i3);
      old_red= colr[i3];
      old_green= colr[i3+1];
      old_blue= colr[i3+2];
    }
    glNormal3fv(norm+i3);
    /* NOTE: It makes no difference in local performance if the coords
       are converted to floats before being stored in the display list */
    glVertex3fv(xyz+i3+6);
  }
  glEnd();
  CHEK_ERROR("yglTstripNoArr");
}

void yglQstripNoArr(long nvert, float *xyz, float *norm, float *colr)
{
  long i3, i6;
  double old_red, old_green, old_blue;

  /* draw the strip */
  if(nvert <= 2) return;

  glBegin(GL_QUAD_STRIP);
  glColor3fv(colr);
  old_red= colr[0];
  old_green= colr[1];
  old_blue= colr[2];
  glColor3fv(colr);
  glVertex3fv(xyz);
  glVertex3fv(xyz+3);
  for(i3= 0, i6= 6; i3 < 3*(nvert-1); i3 += 3, i6 += 6) {
    if(colr[i3] != old_red || colr[i3+1] != old_green || 
       colr[i3+2] != old_blue) {
      glColor3fv(colr+i3);
      old_red= colr[i3];
      old_green= colr[i3+1];
      old_blue= colr[i3+2];
    }
    glNormal3fv(norm+i3);
    /* NOTE: It makes no difference in local performance if the coords
       are converted to floats before being stored in the display list */
    glVertex3fv(xyz+i6);
    glVertex3fv(xyz+i6+3);
  }
  glEnd();
  CHEK_ERROR("yglQstripNoArr");
}

void yglTstripArrNoLite(long nvert, float *xyz, float *colr)
{
  long i3;
  float *arr_colr, *arr_vert;
  long ind_colr, ind_vert;

  /* draw the strip */
  if(nvert <= 2) return;

  /* make an array big enough to hold all the colors and
     vertices for the quad strip */
  ind_colr= ind_vert= 0;
  arr_colr= (float *)p_malloc(3*nvert*sizeof(float));
  arr_vert= (float *)p_malloc(3*nvert*sizeof(float));
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  for(i3= 0; i3 < 3*(nvert-2); i3 += 3) {
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
  }
  for(i3= 0; i3 < 3*nvert; i3 += 3) {
    arr_vert[ind_vert++]= xyz[i3];
    arr_vert[ind_vert++]= xyz[i3+1];
    arr_vert[ind_vert++]= xyz[i3+2];
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles/quads using previously specified colors, vertices */
  glDrawArrays(GL_TRIANGLE_STRIP, 0, nvert);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_vert);
  CHEK_ERROR("yglTstripArrNoLite");
}

void yglQstripArrNoLite(long nvert, float *xyz, float *colr)
{
  long i3, i6;
  float *arr_colr, *arr_vert;
  long ind_colr, ind_vert;

  /* draw the strip */
  if(nvert <= 1) return;

  /* make an array big enough to hold all the colors and
     vertices for the quad strip */
  ind_colr= ind_vert= 0;
  arr_colr= (float *)p_malloc(6*nvert*sizeof(float));
  arr_vert= (float *)p_malloc(6*nvert*sizeof(float));
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  arr_colr[ind_colr++]= colr[0];
  arr_colr[ind_colr++]= colr[1];
  arr_colr[ind_colr++]= colr[2];
  for(i3= 0; i3 < 3*(nvert-1); i3 += 3) {
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
    arr_colr[ind_colr++]= colr[i3];
    arr_colr[ind_colr++]= colr[i3+1];
    arr_colr[ind_colr++]= colr[i3+2];
  }
  for(i6= 0; i6 < 6*nvert; i6 += 6) {
    arr_vert[ind_vert++]= xyz[i6];
    arr_vert[ind_vert++]= xyz[i6+1];
    arr_vert[ind_vert++]= xyz[i6+2];
    arr_vert[ind_vert++]= xyz[i6+3];
    arr_vert[ind_vert++]= xyz[i6+4];
    arr_vert[ind_vert++]= xyz[i6+5];
  }
  /* NOTE: for each enabled array, a corresponding array must be
     defined and have data filled in before a call to glDrawArrays,
     because that call will use data from each enabled array */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_FLOAT, 0, arr_colr);
  glVertexPointer(3, GL_FLOAT, 0, arr_vert);
  /* draw all triangles/quads using previously specified colors, vertices */
  glDrawArrays(GL_QUAD_STRIP, 0, 2*nvert);
  /* free up the storage */
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  p_free(arr_colr);
  p_free(arr_vert);
  CHEK_ERROR("yglQstripArrNoLite");
}

void yglTstripNoArrNoLite(long nvert, float *xyz, float *colr)
{
  long i3;
  double old_red, old_green, old_blue;

  /* draw the strip */
  if(nvert <= 2) return;

  glBegin(GL_TRIANGLE_STRIP);
  glColor3fv(colr);
  old_red= colr[0];
  old_green= colr[1];
  old_blue= colr[2];
  glVertex3fv(xyz);
  glVertex3fv(xyz+3);
  for(i3= 0; i3 < 3*(nvert-2); i3 += 3) {
    if(colr[0] != old_red || colr[1] != old_green || 
       colr[2] != old_blue) {
      glColor3fv(colr);
      old_red= colr[0];
      old_green= colr[1];
      old_blue= colr[2];
    }
    colr += 3;
    /* NOTE: It makes no difference in local performance if the coords
       are converted to floats before being stored in the display list */
    glVertex3fv(xyz+i3+6);
  }
  glEnd();
  CHEK_ERROR("yglTstripNoArrNoLite");
}

void yglTstripNoArrNoLiteAlpha(long nvert, float *xyz, float *colr)
{
  long i3;
  float old_red, old_green, old_blue, old_alpha;

  /* draw the strip */
  if(nvert <= 2) return;

  glBegin(GL_TRIANGLE_STRIP);
  old_red= colr[0];
  old_green= colr[1];
  old_blue= colr[2];
  old_alpha= colr[3];
  glColor4fv(colr);
  glVertex3fv(xyz);
  glVertex3fv(xyz+3);
  for(i3= 0; i3 < 3*(nvert-2); i3 += 3) {
    if(colr[0] != old_red || colr[1] != old_green || 
       colr[2] != old_blue || colr[3] != old_alpha) {
      glColor4fv(colr);
      old_red= colr[0];
      old_green= colr[1];
      old_blue= colr[2];
      old_alpha= colr[3];
    }
    colr += 4;
    /* NOTE: It makes no difference in local performance if the coords
       are converted to floats before being stored in the display list */
    glVertex3fv(xyz+i3+6);
  }
  glEnd();
  CHEK_ERROR("yglTstripNoArrNoLiteAlpha");
}

void yglQstripNoArrNoLite(long nvert, float *xyz, float *colr)
{
  long i3, i6;
  double old_red, old_green, old_blue;

  /* draw the strip */
  if(nvert <= 2) return;

  glBegin(GL_QUAD_STRIP);
  glColor3fv(colr);
  old_red= colr[0];
  old_green= colr[1];
  old_blue= colr[2];
  glVertex3fv(xyz);
  glVertex3fv(xyz+3);
  for(i3= 0, i6= 6; i3 < 3*(nvert-1); i3 += 3, i6 += 6) {
    if(colr[i3] != old_red || colr[i3+1] != old_green || 
       colr[i3+2] != old_blue) {
      glColor3fv(colr+i3);
      old_red= colr[i3];
      old_green= colr[i3+1];
      old_blue= colr[i3+2];
    }
    /* NOTE: It makes no difference in local performance if the coords
       are converted to floats before being stored in the display list */
    glVertex3fv(xyz+i6);
    glVertex3fv(xyz+i6+3);
  }
  glEnd();
  CHEK_ERROR("yglQstripNoArrNoLite");
}

void yglTstripsNdx(long nstrip, long numedg, long ntri, long *len, long *ndx,
                   float *xyz, float *norm, float *colr, long edge)
{
  /* The input is a list of triangle strips. There is one color per triangle.
     There is one normal per vertex.
     The coordinates and normals are accessed via ndx, an array with
     one index into xyz and norm per vertex in the tri strip.
  */
  long i, num, i3, indv;
  float old_red= -1, old_green= -1, old_blue= -1;

  yglSetPolyMode(edge);
  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();
  /* smooth shading */
  yglSetColorType(1);

  for(i= 0, indv= 0; i < nstrip; i++) {
    num= len[i];  /* number of vertices in this strip */
    /* draw the strip */
    DEMAND(num > 2, "triangle strip with less than 3 vertices in yglTstripsNdx")
    glBegin(GL_TRIANGLE_STRIP);
    if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue) {
      old_red= colr[0];
      old_green= colr[1];
      old_blue= colr[2];
      glColor3fv(colr);
    }
    glNormal3fv(norm+3*ndx[indv]);
    glVertex3fv(xyz+3*ndx[indv]);
    glNormal3fv(norm+3*ndx[indv+1]);
    glVertex3fv(xyz+3*ndx[indv+1]);
    indv += 2;
    for(i3= 0; i3 < num-2; i3++) {
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue) {
        old_red= colr[0];
        old_green= colr[1];
        old_blue= colr[2];
        glColor3fv(colr);
      }
      glNormal3fv(norm+3*ndx[indv]);
      glVertex3fv(xyz+3*ndx[indv]);
      indv++;
      colr += 3;
    }
    glEnd();
  }
  CHEK_ERROR("yglTstripsNdx");
}

void yglTstripsAlphaNdx(long nstrip, long numedg, long ntri, long *len, 
                        long *ndx, float *xyz, float *norm, float *colr,
                        long edge)
{
  /* The input is a list of triangle strips. There is one color per triangle.
     There is one normal per vertex.
     The coordinates and normals are accessed via ndx, an array with
     one index into xyz and norm per vertex in the tri strip.
  */
  long i, num, i3, indv;
  float old_red= -1, old_green= -1, old_blue= -1, old_alpha= -1;

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);
  yglSetPolyMode(edge);
  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();
  /* smooth shading */
  yglSetColorType(1);

  for(i= 0, indv= 0; i < nstrip; i++) {
    num= len[i];  /* number of vertices in this strip */
    /* draw the strip */
    DEMAND(num > 2, "triangle strip with less than 3 vertices in yglTstripsNdx")
    glBegin(GL_TRIANGLE_STRIP);
    if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue) {
      old_red= colr[0];
      old_green= colr[1];
      old_blue= colr[2];
      old_alpha= colr[3];
      glColor4fv(colr);
    }
    glNormal3fv(norm+ndx[indv]);
    glVertex3fv(xyz+ndx[indv]);
    glNormal3fv(norm+ndx[indv+1]);
    glVertex3fv(xyz+ndx[indv+1]);
    indv += 2;
    for(i3= 0; i3 < num-2; i3++) {
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      if(colr[0] != old_red || colr[1] != old_green || colr[2] != old_blue 
         || colr[3] != old_alpha) {
        old_red= colr[0];
        old_green= colr[1];
        old_blue= colr[2];
        old_alpha= colr[3];
        glColor4fv(colr);
      }
      glNormal3fv(norm+ndx[indv]);
      glVertex3fv(xyz+ndx[indv]);
      indv++;
      colr += 4;
    }
    glEnd();
  }
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  CHEK_ERROR("yglTstripsAlphaNdx");
}
