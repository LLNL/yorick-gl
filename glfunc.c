/*
 * $Id: glfunc.c,v 1.1 2005-09-18 22:07:49 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include "ydata.h"

#include "glfunc.h"
#include "glcode.h"
#include "glBasic.h"
#include "glStrips.h"
#include "glMouse.h"
#include "glWrappers.h"
#include <math.h>

/* Yorick tends to work with doubles, but floats are better for OpenGL.
   Use floats here.
*/

void yglForceWin3d(void)
{
  /* force the creation of a 3D window if none is currently open */
  if(!glCurrWin3d) {
    yglWin3d(0, 500, 500);
   }
}

float yglGetFov3d(void)
{
  /* Return the field of view in degrees. */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  return glCurrWin3d->fov;
}

void yglSetFov3d(float fov)
{
  /* require the field-of-view to be somewhat reasonable */
  ASSERT(fov >= 1.0e-4 && fov <= 90.0, "gl_set_fov: field of view must be between 1.0e-4 and 90 degrees")
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  glCurrWin3d->fov= fov;
}

long yglGetWidth3d(void)
{
  /* Return the width in pixels of the OpenGL window. */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  return glCurrWin3d->width;
}

long yglGetHite3d(void)
{
  /* Return the height in pixels of the OpenGL window. */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  return glCurrWin3d->hite;
}

int yglHasTex3d(void)
{
  int res;

  /* Non-zero indicates that this OpenGL supports 3D textures. */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  ygl_fpemask(0);
  res= yglQueryTex3d(glCurrWin3d);
  ygl_fpemask(1);
  return res;
}

int yglHasTexcube3d(void)
{
  int res;

  /* Non-zero indicates that this OpenGL supports cube map textures. */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  ygl_fpemask(0);
  res= yglQueryTexCube();
  ygl_fpemask(1);
  return res;
}

void yglAlwaysShowObj3d(long flag)
{
  /* Object will be drawn even when rotating if flag is non-zero */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(flag) glCurrWin3d->always_show_obj= 1;
  else glCurrWin3d->always_show_obj= 0;
}

void yglUseList3d(long flag)
{
  /* Turn on use of an OpenGL display list on if non-zero */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(flag) glCurrWin3d->use_list= 1;
  else glCurrWin3d->use_list= 0;
}

void yglUseArray3d(long flag)
{
  /* Turn use of OpenGL arrays on if non-zero */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(flag) glCurrWin3d->use_array= 1;
  else glCurrWin3d->use_array= 0;
}

void yglMsmovVal3d(double val)
{
  ASSERT(val >= 2.0 && val <= 200.0, "gl_mov_val must be between 2 and 200")
  ygl_ms_mov_val= val;
}

void yglGetBackRGB3d(double *rgb)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  rgb[0]= glCurrWin3d->back_red;
  rgb[1]= glCurrWin3d->back_green;
  rgb[2]= glCurrWin3d->back_blue;
}

void yglBackRGB3d(double *rgb)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(rgb[0] >= 0.0 && rgb[0] <= 1.0) glCurrWin3d->back_red= (float) rgb[0];
  if(rgb[1] >= 0.0 && rgb[1] <= 1.0) glCurrWin3d->back_green= (float) rgb[1];
  if(rgb[2] >= 0.0 && rgb[2] <= 1.0) glCurrWin3d->back_blue= (float) rgb[2];
}

void yglBackRGBA3d(double *rgba)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(rgba[0] >= 0.0 && rgba[0] <= 1.0) glCurrWin3d->back_red= (float) rgba[0];
  if(rgba[1] >= 0.0 && rgba[1] <= 1.0) glCurrWin3d->back_green= (float) rgba[1];
  if(rgba[2] >= 0.0 && rgba[2] <= 1.0) glCurrWin3d->back_blue= (float) rgba[2];
  if(rgba[3] >= 0.0 && rgba[3] <= 1.0) glCurrWin3d->back_alpha= (float) rgba[3];
}

void yglCageStyle3d(long style)
{
  /* zero means no cage, negative means use data bounds, positive means 
     use limits supplied by user in the most recent call to gl_cage_limits. */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(style != glCurrWin3d->cage_style) {
    glCurrWin3d->cage_style= style;
    glCurrWin3d->cage_seq_num++;  /* changing the cage style changes the plot */
  }
}

void yglGetCageRGB3d(double *rgb)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  rgb[0]= glCurrWin3d->cage_red;
  rgb[1]= glCurrWin3d->cage_green;
  rgb[2]= glCurrWin3d->cage_blue;
}

void yglCageRGB3d(double *rgb)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(rgb[0] >= 0.0 && rgb[0] <= 1.0) glCurrWin3d->cage_red= (float) rgb[0];
  if(rgb[1] >= 0.0 && rgb[1] <= 1.0) glCurrWin3d->cage_green= (float) rgb[1];
  if(rgb[2] >= 0.0 && rgb[2] <= 1.0) glCurrWin3d->cage_blue= (float) rgb[2];
}

void yglCageRGBA3d(double *rgba)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(rgba[0] >= 0.0 && rgba[0] <= 1.0) glCurrWin3d->cage_red= (float) rgba[0];
  if(rgba[1] >= 0.0 && rgba[1] <= 1.0) glCurrWin3d->cage_green= (float) rgba[1];
  if(rgba[2] >= 0.0 && rgba[2] <= 1.0) glCurrWin3d->cage_blue= (float) rgba[2];
  if(rgba[3] >= 0.0 && rgba[3] <= 1.0) glCurrWin3d->cage_alpha= (float) rgba[3];
}

void yglGetGridRGB3d(double *rgb)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  rgb[0]= glCurrWin3d->grid_red;
  rgb[1]= glCurrWin3d->grid_green;
  rgb[2]= glCurrWin3d->grid_blue;
}

void yglGridRGB3d(double *rgb)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(rgb[0] >= 0.0 && rgb[0] <= 1.0) glCurrWin3d->grid_red= (float) rgb[0];
  if(rgb[1] >= 0.0 && rgb[1] <= 1.0) glCurrWin3d->grid_green= (float) rgb[1];
  if(rgb[2] >= 0.0 && rgb[2] <= 1.0) glCurrWin3d->grid_blue= (float) rgb[2];
}

void yglGridRGBA3d(double *rgba)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(rgba[0] >= 0.0 && rgba[0] <= 1.0) glCurrWin3d->grid_red= (float) rgba[0];
  if(rgba[1] >= 0.0 && rgba[1] <= 1.0) glCurrWin3d->grid_green= (float) rgba[1];
  if(rgba[2] >= 0.0 && rgba[2] <= 1.0) glCurrWin3d->grid_blue= (float) rgba[2];
  if(rgba[3] >= 0.0 && rgba[3] <= 1.0) glCurrWin3d->grid_alpha= (float) rgba[3];
}

void yglGetCageLimits3d(double *lims)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  lims[0]= glCurrWin3d->cage_xmin;
  lims[1]= glCurrWin3d->cage_xmax;
  lims[2]= glCurrWin3d->cage_ymin;
  lims[3]= glCurrWin3d->cage_ymax;
  lims[4]= glCurrWin3d->cage_zmin;
  lims[5]= glCurrWin3d->cage_zmax;
}

void yglCageLimits3d(double *lims)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  if(lims[0] <= lims[1]) {
    glCurrWin3d->cage_xmin= (float)lims[0]; glCurrWin3d->cage_xmax= (float)lims[1];
  } else {
    glCurrWin3d->cage_xmin= (float)lims[1]; glCurrWin3d->cage_xmax= (float)lims[0];
  }
  if(lims[2] <= lims[3]) {
    glCurrWin3d->cage_ymin= (float)lims[2]; glCurrWin3d->cage_ymax= (float)lims[3];
  } else {
    glCurrWin3d->cage_ymin= (float)lims[3]; glCurrWin3d->cage_ymax= (float)lims[2];
  }
  if(lims[4] <= lims[5]) {
    glCurrWin3d->cage_zmin= (float)lims[4]; glCurrWin3d->cage_zmax= (float)lims[5];
  } else {
    glCurrWin3d->cage_zmin= (float)lims[5]; glCurrWin3d->cage_zmax= (float)lims[4];
  }
  /* increase the sequence number so that this change will
     be reflected in the next plot (only if active) */
  if(glCurrWin3d->cage_style > 0) glCurrWin3d->cage_seq_num++;
}

void yglLookat3d(double *eye, double *center, double *up)
{
  double dot, len;

  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  glCurrWin3d->eye[0]= eye[0];
  glCurrWin3d->eye[1]= eye[1];
  glCurrWin3d->eye[2]= eye[2];
  glCurrWin3d->center[0]= center[0];
  glCurrWin3d->center[1]= center[1];
  glCurrWin3d->center[2]= center[2];
  /* compute the unit view vector and viewing distance */
  glCurrWin3d->view[0]= glCurrWin3d->eye[0]-glCurrWin3d->center[0];
  glCurrWin3d->view[1]= glCurrWin3d->eye[1]-glCurrWin3d->center[1];
  glCurrWin3d->view[2]= glCurrWin3d->eye[2]-glCurrWin3d->center[2];
  glCurrWin3d->viewdist= sqrt(glCurrWin3d->view[0]*glCurrWin3d->view[0]
                         +glCurrWin3d->view[1]*glCurrWin3d->view[1]
                         +glCurrWin3d->view[2]*glCurrWin3d->view[2]);
  glCurrWin3d->view[0] /= glCurrWin3d->viewdist;
  glCurrWin3d->view[1] /= glCurrWin3d->viewdist;
  glCurrWin3d->view[2] /= glCurrWin3d->viewdist;
  /* force the up vector to have unit length and be normal
     to the viewing direction. */
  dot= up[0]*glCurrWin3d->view[0]+up[1]*glCurrWin3d->view[1]+up[2]*glCurrWin3d->view[2];
  glCurrWin3d->up[0]= up[0]-dot*glCurrWin3d->view[0];
  glCurrWin3d->up[1]= up[1]-dot*glCurrWin3d->view[1];
  glCurrWin3d->up[2]= up[2]-dot*glCurrWin3d->view[2];
  len= sqrt(glCurrWin3d->up[0]*glCurrWin3d->up[0]+glCurrWin3d->up[1]*glCurrWin3d->up[1]+
            glCurrWin3d->up[2]*glCurrWin3d->up[2]);
  if(len < 1.0e-4) {
    /* a bad up vector was supplied - an impossible error if the interpreted
       wrapper function was called.
       Try a "random" direction. 
    */
    glCurrWin3d->up[0]= 0.5;
    glCurrWin3d->up[1]= 1.0/sqrt(2.0);
    glCurrWin3d->up[2]= 0.5;
  } else {
    glCurrWin3d->up[0] /= len;
    glCurrWin3d->up[1] /= len;
    glCurrWin3d->up[2] /= len;
  }
}

int yglUpdateList3d(void)
{
  /* Return a flag indicating whether objects need to be
     added to the OpenGL display list. */
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  return (glCurrWin3d->seq_num > glCurrWin3d->list_num);
}

void yglIncSeq3d(void)
{
  if(!glCurrWin3d) {
    yglForceWin3d();
  }
  glCurrWin3d->seq_num++;
}

void yglArsum3d(long nx0, long ny0, long nz0, long nsx, long nsy, long nsz, 
              double *vin, double *vout)
{
  long nxout, nyout, nzout, i, j, k, iout, jout, kout, ndin, ndout;
  
  /* sum up nsx-by-nsy-by-nsz zone chunks of vin and
     store the result in vout.
     This function thus makes an integral preserving,
     lower resolution version of the input array.
     If the caller wants a smoothed version of vin,
     he/she should divide vout by nsx*nsy*nsz. */

  nxout= nx0/nsx;  nyout= ny0/nsy;  nzout= nz0/nsz;
  /* clear the output array */
  for(k= 0; k < nzout; k++) {
    for(j= 0; j < nyout; j++) {
      for(i= 0; i < nxout; i++) {
        ndout= i/nsx+j/nsy*nxout+k/nsz*nxout*nyout;
        vout[ndout]= 0.0;
      }
    }
  }
  /* the following nested loops sums the input array in memory
     order into the output array. */
  for(k= 0; k < nz0; k++) {
    for(j= 0; j < ny0; j++) {
      for(i= 0; i < nx0; i++) {
        iout= i/nsx; jout= j/nsy; kout= k/nsz;
        ndout= iout+jout*nxout+kout*nxout*nyout;
        ndin= i+j*nx0+k*nx0*ny0;
        vout[ndout] += vin[ndin];
      }
    }
  }
}
