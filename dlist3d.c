/*
 * $Id: dlist3d.c,v 1.3 2006-10-19 14:48:19 dhmunro Exp $
 * Implement functions used for manipulating 3D display lists.
 * The functions in this file maintains the display list for the 3D graphics 
 * package in Yorick.
 * The display routines MUST all be thread safe to permit this package to work
 * on Power Walls and other "multi-graphics-card" displays.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include "glcode.h"
#include "glBasic.h"
#include "dlist3d.h"
#include "glfunc.h"
#include "glStrips.h"
#include "glWrappers.h"
#include "ydata.h"
#include "pstdlib.h"
#include <string.h>
#include <stdlib.h>

extern long yglGetBoundsList3d(yBox3D *box, yList3d_Elem *elem);

yList3d_Elem *yListDirectHead= 0;
yList3d_Elem *yListCachedHead= 0;
int yDrawBBox3d= 0;
int alpha_pass= 0;

void yglDrawCurr3d(void)
{
  yglDraw3d(glCurrWin3d);
}

void yglDraw3d(glWinProp *theWin3d)
{
  glWinProp *oldWin3d= glCurrWin3d;

  /* prepare to draw, draw both cached and direct display lists,
     and finish drawing.
  */
  if(!theWin3d) return;
  ygl_fpemask(0);
  glCurrWin3d= theWin3d;
  yglPrepDraw(theWin3d);
  yglUpdateLight();
  yglPrepContext();
  /* draw all items in the "cached" display list */
  yglDrawListCache3d();
  /* draw all items in the "direct" display list */
  yglDrawListDirect3d();
  yglDrawCage();
  yglDrawGnomon();
  yglFinFrame();
  theWin3d->dirty= 0;
  glCurrWin3d= oldWin3d;
  ygl_fpemask(1);
}

void yglClearList3d(void)
{
  ygl_fpemask(0);
  yglClearCachedList3d();
  yglClearDirectList3d();
  ygl_fpemask(1);
}

void yglClearCachedList3d(void)
{
  yList3d_Elem *elem;
  
  /* empty the display list and free all storage associated with it */
  while(yListCachedHead) {
    elem= yListCachedHead;
    yListCachedHead= yListCachedHead->next;
    p_free(elem->data);
    p_free(elem);
  }
  /* the OpenGL display list is now out-of-date */
  if(glCurrWin3d && (glCurrWin3d->seq_num <= glCurrWin3d->list_num) ) {
    glCurrWin3d->seq_num++;
  }
}

void yglClearDirectList3d(void)
{
  yList3d_Elem *elem;
  
  /* empty the display list and free all storage associated with it */
  while(yListDirectHead) {
    elem= yListDirectHead;
    yListDirectHead= yListDirectHead->next;
    p_free(elem->data);
    p_free(elem);
  }
}

void yglDrawList3d(void)
{
  /* draw all items in the "cached" display list */
  yglDrawListCache3d();
  /* draw all items in the "direct" display list */
  yglDrawListDirect3d();
}

void yglDrawListCache3d(void)
{
  yList3d_Elem *elem;
  
  if(glCurrWin3d && (glCurrWin3d->seq_num > glCurrWin3d->list_num) ) {
    /* need a new display list */
    /* draw all items in the "cached" display list */
    yglPrepList();
    /* in the first pass draw all opaque objects */
    alpha_pass= 0;
    elem= yListCachedHead;
    while(elem) {
      (*(elem->func))(0, elem->data);
      elem= elem->next;
    }
    /* In the second pass draw all translucent objects.
       NOTE: objects with an alpha of 1.0 are still "translucent"
       and need not produce the same effect as an opaque object. */
    alpha_pass= 1;
    elem= yListCachedHead;
    while(elem) {
      (*(elem->func))(0, elem->data);
      elem= elem->next;
    }
    alpha_pass= 0;
  }
  /* if necessary, finish the display list, then draw it */
  yglFinCache();
}

void yglDrawListDirect3d(void)
{
  yList3d_Elem *elem;
  
  /* draw all items in the "direct" display list */
  elem= yListDirectHead;
  while(elem) {
    (*(elem->func))(0, elem->data);
    elem= elem->next;
  }
}

long yglGetBounds3d(yBox3D *box)
{
  yBox3D boxCached, boxDirect;
  long hasDirect, hasCached;

  /* if the bounding box for this scene has already been
     computed, return it. */
  if(!glCurrWin3d) return 0;  /* indicate failure */
  if((glCurrWin3d->BoxSeqNum >= glCurrWin3d->seq_num) &&
     (glCurrWin3d->cage_seq_num <= glCurrWin3d->cage_state) ) {
    *box= glCurrWin3d->boxAll;
    return 1;
  }
  /* Get the bounds of the objects in the scene */
  hasCached= yglGetBoundsCachedList3d(&boxCached);
  hasDirect= yglGetBoundsDirectList3d(&boxDirect);
  if(hasCached) {
    if(hasDirect) {
      /* merge the two bounding boxes */
      glCurrWin3d->boxAll= boxCached;
      if(glCurrWin3d->boxAll.xmin > boxDirect.xmin) glCurrWin3d->boxAll.xmin= boxDirect.xmin;
      if(glCurrWin3d->boxAll.xmax < boxDirect.xmax) glCurrWin3d->boxAll.xmax= boxDirect.xmax;
      if(glCurrWin3d->boxAll.ymin > boxDirect.ymin) glCurrWin3d->boxAll.ymin= boxDirect.ymin;
      if(glCurrWin3d->boxAll.ymax < boxDirect.ymax) glCurrWin3d->boxAll.ymax= boxDirect.ymax;
      if(glCurrWin3d->boxAll.zmin > boxDirect.zmin) glCurrWin3d->boxAll.zmin= boxDirect.zmin;
      if(glCurrWin3d->boxAll.zmax < boxDirect.zmax) glCurrWin3d->boxAll.zmax= boxDirect.zmax;
    } else {
      /* only cached objects */
      glCurrWin3d->boxAll= boxCached;
    }
  } else {
    if(hasDirect) {
      /* only direct objects */
      glCurrWin3d->boxAll= boxDirect;
    } else {
      /* no objects, so return "failure" and an empty box */
      glCurrWin3d->boxAll.xmin= 0.0; glCurrWin3d->boxAll.xmax= 0.0;
      glCurrWin3d->boxAll.ymin= 0.0; glCurrWin3d->boxAll.ymax= 0.0;
      glCurrWin3d->boxAll.zmin= 0.0; glCurrWin3d->boxAll.zmax= 0.0;
      *box= glCurrWin3d->boxAll;
      return 0;
    }
  }
  if(glCurrWin3d->cage_style > 0) {
    /* modify the box to include the "cage" set by the user */
    if(glCurrWin3d->boxAll.xmin > glCurrWin3d->cage_xmin) glCurrWin3d->boxAll.xmin= glCurrWin3d->cage_xmin;
    if(glCurrWin3d->boxAll.xmax < glCurrWin3d->cage_xmax) glCurrWin3d->boxAll.xmax= glCurrWin3d->cage_xmax;
    if(glCurrWin3d->boxAll.ymin > glCurrWin3d->cage_ymin) glCurrWin3d->boxAll.ymin= glCurrWin3d->cage_ymin;
    if(glCurrWin3d->boxAll.ymax < glCurrWin3d->cage_ymax) glCurrWin3d->boxAll.ymax= glCurrWin3d->cage_ymax;
    if(glCurrWin3d->boxAll.zmin > glCurrWin3d->cage_zmin) glCurrWin3d->boxAll.zmin= glCurrWin3d->cage_zmin;
    if(glCurrWin3d->boxAll.zmax < glCurrWin3d->cage_zmax) glCurrWin3d->boxAll.zmax= glCurrWin3d->cage_zmax;
  }
  *box= glCurrWin3d->boxAll;
  glCurrWin3d->BoxSeqNum= glCurrWin3d->seq_num;
  glCurrWin3d->cage_state= glCurrWin3d->cage_seq_num;
  return 1;
}

long yglGetBoundsDirectList3d(yBox3D *box)
{
  long i;
  ygl_fpemask(0);
  i = yglGetBoundsList3d(box, yListDirectHead);
  ygl_fpemask(1);
  return i;
}

long yglGetBoundsCachedList3d(yBox3D *box)
{
  long i;
  ygl_fpemask(0);
  i = yglGetBoundsList3d(box, yListCachedHead);
  ygl_fpemask(1);
  return i;
}

long yglGetBoundsList3d(yBox3D *box, yList3d_Elem *elem)
{
  double xmin, xmax, ymin, ymax, zmin, zmax;
  
  /* find the bounding box containing all items in the specified display list */
  if(elem) {
    xmin= (elem->box).xmin;
    xmax= (elem->box).xmax;
    ymin= (elem->box).ymin;
    ymax= (elem->box).ymax;
    zmin= (elem->box).zmin;
    zmax= (elem->box).zmax;
    elem= elem->next;
    while(elem) {
      if(xmin > (elem->box).xmin) xmin= (elem->box).xmin;
      if(xmax < (elem->box).xmax) xmax= (elem->box).xmax;
      if(ymin > (elem->box).ymin) ymin= (elem->box).ymin;
      if(ymax < (elem->box).ymax) ymax= (elem->box).ymax;
      if(zmin > (elem->box).zmin) zmin= (elem->box).zmin;
      if(zmax < (elem->box).zmax) zmax= (elem->box).zmax;
      elem= elem->next;
    }
    box->xmin= xmin;
    box->xmax= xmax;
    box->ymin= ymin;
    box->ymax= ymax;
    box->zmin= zmin;
    box->zmax= zmax;
    return 1;
  } else {
    box->xmin= 0.0;
    box->xmax= 0.0;
    box->ymin= 0.0;
    box->ymax= 0.0;
    box->zmin= 0.0;
    box->zmax= 0.0;
    return 0;
  }
}

yList3d_Elem *yglNewCachedList3dElem(void)
{
  yList3d_Elem *elem;
  
  /* add another element at the head of the display list (it is
     a linked list) */
  elem= (yList3d_Elem *) p_malloc(sizeof(yList3d_Elem));
  elem->next= yListCachedHead;
  yListCachedHead= elem;
  return elem;
}

yList3d_Elem *yglNewDirectList3dElem(void)
{
  yList3d_Elem *elem;
  
  /* add another element at the head of the display list (it is
     a linked list) */
  elem= (yList3d_Elem *) p_malloc(sizeof(yList3d_Elem));
  elem->next= yListDirectHead;
  yListDirectHead= elem;
  return elem;
}

void yglSetLims3d(yList3d_Elem *elem, long nvert, float *xyz)
{
  long i;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  float xx, yy, zz;
  
  if(nvert <= 0) return;
  xmin= xmax= xyz[0];
  ymin= ymax= xyz[1];
  zmin= zmax= xyz[2];
  for(i= 1; i < nvert; i++) {
    xx= xyz[3*i];
    yy= xyz[3*i+1];
    zz= xyz[3*i+2];
    if(xx < xmin) xmin= xx;
    if(xx > xmax) xmax= xx;
    if(yy < ymin) ymin= yy;
    if(yy > ymax) ymax= yy;
    if(zz < zmin) zmin= zz;
    if(zz > zmax) zmax= zz;
  }
  elem->box.xmin= xmin;
  elem->box.xmax= xmax;
  elem->box.ymin= ymin;
  elem->box.ymax= ymax;
  elem->box.zmin= zmin;
  elem->box.zmax= zmax;
}

void yglDrawGnomon(void)
{
  /* draw an indicator of axis orientation */
}

void yglDrawCage(void)
{
  float pt1[3], pt2[3], pt3[3];
  yBox3D boxAll;

  /* Determine the three furthest planes, given the current
     viewpoint. Draw them. If the "cage" uses the data bounds,
     put them in place now.
   */
  if(!glCurrWin3d) return;
  if(!glCurrWin3d->cage_style) {
    /* do not draw a cage */
    return;
  }
  if(glCurrWin3d->cage_style < 0) {
    /* Use the bounds of the objects in the scene. 
       If no objects, just return */
    if(!yglGetBounds3d(&boxAll)) {
      return;
    }
    glCurrWin3d->cage_xmin= (float)boxAll.xmin; glCurrWin3d->cage_xmax= (float)boxAll.xmax;
    glCurrWin3d->cage_ymin= (float)boxAll.ymin; glCurrWin3d->cage_ymax= (float)boxAll.ymax;
    glCurrWin3d->cage_zmin= (float)boxAll.zmin; glCurrWin3d->cage_zmax= (float)boxAll.zmax;
  }
  if(glCurrWin3d->view[0] >= 0.0) {
    /* draw the plane at cage_xmin because the eye
       is at a larger x than the center */
    pt1[0]= glCurrWin3d->cage_xmin; pt1[1]= glCurrWin3d->cage_ymin; pt1[2]= glCurrWin3d->cage_zmin;
    pt2[0]= glCurrWin3d->cage_xmin; pt2[1]= glCurrWin3d->cage_ymin; pt2[2]= glCurrWin3d->cage_zmax;
    pt3[0]= glCurrWin3d->cage_xmin; pt3[1]= glCurrWin3d->cage_ymax; pt3[2]= glCurrWin3d->cage_zmax;
  } else {
    pt1[0]= glCurrWin3d->cage_xmax; pt1[1]= glCurrWin3d->cage_ymin; pt1[2]= glCurrWin3d->cage_zmin;
    pt2[0]= glCurrWin3d->cage_xmax; pt2[1]= glCurrWin3d->cage_ymin; pt2[2]= glCurrWin3d->cage_zmax;
    pt3[0]= glCurrWin3d->cage_xmax; pt3[1]= glCurrWin3d->cage_ymax; pt3[2]= glCurrWin3d->cage_zmax;
  }
  draw_plane(pt1, pt2, pt3, glCurrWin3d->num_zgrid, glCurrWin3d->num_ygrid);
  if(glCurrWin3d->view[1] >= 0.0) {
    /* draw the plane at cage_ymin because the eye
       is at a larger y than the center */
    pt1[0]= glCurrWin3d->cage_xmin; pt1[1]= glCurrWin3d->cage_ymin; pt1[2]= glCurrWin3d->cage_zmin;
    pt2[0]= glCurrWin3d->cage_xmin; pt2[1]= glCurrWin3d->cage_ymin; pt2[2]= glCurrWin3d->cage_zmax;
    pt3[0]= glCurrWin3d->cage_xmax; pt3[1]= glCurrWin3d->cage_ymin; pt3[2]= glCurrWin3d->cage_zmax;
  } else {
    pt1[0]= glCurrWin3d->cage_xmin; pt1[1]= glCurrWin3d->cage_ymax; pt1[2]= glCurrWin3d->cage_zmin;
    pt2[0]= glCurrWin3d->cage_xmin; pt2[1]= glCurrWin3d->cage_ymax; pt2[2]= glCurrWin3d->cage_zmax;
    pt3[0]= glCurrWin3d->cage_xmax; pt3[1]= glCurrWin3d->cage_ymax; pt3[2]= glCurrWin3d->cage_zmax;
  }
  draw_plane(pt1, pt2, pt3, glCurrWin3d->num_zgrid, glCurrWin3d->num_xgrid);
  if(glCurrWin3d->view[2] >= 0.0) {
    /* draw the plane at cage_zmin because the eye
       is at a larger z than the center */
    pt1[0]= glCurrWin3d->cage_xmin; pt1[1]= glCurrWin3d->cage_ymin; pt1[2]= glCurrWin3d->cage_zmin;
    pt2[0]= glCurrWin3d->cage_xmin; pt2[1]= glCurrWin3d->cage_ymax; pt2[2]= glCurrWin3d->cage_zmin;
    pt3[0]= glCurrWin3d->cage_xmax; pt3[1]= glCurrWin3d->cage_ymax; pt3[2]= glCurrWin3d->cage_zmin;
  } else {
    pt1[0]= glCurrWin3d->cage_xmin; pt1[1]= glCurrWin3d->cage_ymin; pt1[2]= glCurrWin3d->cage_zmax;
    pt2[0]= glCurrWin3d->cage_xmin; pt2[1]= glCurrWin3d->cage_ymax; pt2[2]= glCurrWin3d->cage_zmax;
    pt3[0]= glCurrWin3d->cage_xmax; pt3[1]= glCurrWin3d->cage_ymax; pt3[2]= glCurrWin3d->cage_zmax;
  }
  draw_plane(pt1, pt2, pt3, glCurrWin3d->num_ygrid, glCurrWin3d->num_xgrid);
}

void yglPolys3d(long npolys, long *len, double *xyz, double *norm,
                double *colr, long edge, long smooth, long do_light)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yPoly3dData *data;
  long size, nvert, i, nval;
  float *dxyz, *dnorm, *dcolr;
  
  /* There is one color per polygon (i.e. npolys colors)
     The number of polygon vertices is the sum of the npolys
     elements of len[].
     There is one coordinate and one normal ???????????????? per vertex */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawPolys3d;
  /* compute the size of the data structure plus all data */
  nvert= 0;
  for(i= 0; i < npolys; i++) {
    nvert += len[i];
  }
  size= sizeof(yPoly3dData)+sizeof(long)*npolys+sizeof(float)*(2*3*nvert+3*npolys);
  data= (yPoly3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->npolys= npolys;
  data->edge= edge;
  data->smooth= smooth;
  data->do_light= do_light;
  data->len= (long *) ((char *)data+sizeof(yPoly3dData));
  data->xyz= dxyz= (float *) (data->len+npolys);
  data->norm= dnorm= data->xyz+3*nvert;
  data->colr= dcolr= data->norm+3*nvert;
  /* copy data into the new storage space */
  memcpy(data->len, len, sizeof(long)*npolys);
  /* copy data into the new storage space */
  nval= 3*nvert;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
    dnorm[i]= (float) norm[i];
  }
  nval= 3*npolys;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, nvert, data->xyz);
  ygl_fpemask(1);
}

void yglGlyphs3d(long nglyph, double *origin, double *scal,
                  double *theta, double *phi, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yGlyph3dData *data;
  float *dorigin, *dscal, *dtheta, *dphi, *dcolr;
  long i, size, nval;
  
  /* There is an origin, radius, ellipticity, color, and 
     direction (theta and phi) per glyph */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawGlyphs3d;
  /* compute the size of the data structure plus all data */
  size= sizeof(yGlyph3dData)+sizeof(float)*nglyph*(3+1+1+1+3);
  data= (yGlyph3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nglyph= nglyph;
  data->origin= dorigin= (float *) ((char *)data+sizeof(yGlyph3dData));
  data->scal= dscal= data->origin+3*nglyph;
  data->theta=  dtheta= data->scal+nglyph;
  data->phi=    dphi= data->theta+nglyph;
  data->colr=   dcolr= data->phi+nglyph;
  /* copy data into the new storage space with type conversion */
  for(i= 0; i < nglyph; i++) {
    dscal[i]=  (float) scal[i];
    dtheta[i]= (float) theta[i];
    dphi[i]=   (float) phi[i];
  }
  nval= 3*nglyph;
  for(i= 0; i < nval; i++) {
    dorigin[i]= (float) origin[i];
    dcolr[i]=   (float) colr[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, nglyph, data->origin);
  ygl_fpemask(1);
}

void yglCells3d(long nx, long ny, double *corners,
                double *norm, double *colr, long do_alpha)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yCell3dData *data;
  long i, size, nval;
  float *dxyz, *dnorm, *dcolr;
  
  /* there are nx-by-ny colors, three 3D coordinates and one 3D normal.
     NOTE: given three coordinates, the normal is redundant.
     !!!!!! shouldn't use a normal, just "emissive" color !!!! */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawCells3d;
  /* compute the size of the data structure plus all data */
  size= sizeof(yCell3dData)+sizeof(float)*3*nx*ny+sizeof(float)*(3*3+3);
  data= (yCell3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nx= nx;
  data->ny= ny;
  data->do_alpha= do_alpha;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(yCell3dData));
  data->norm= dnorm= data->xyz+3*3;
  data->colr= dcolr= data->norm+3;
  /* copy data into the new storage space */
  nval= 3;
  for(i= 0; i < nval; i++) {
    dnorm[i]= (float) norm[i];
  }
  nval= 3*3;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) corners[i];
  }
  nval= 3*nx*ny;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, 2, data->xyz);
  ygl_fpemask(1);
}

void yglPlm3d(long nx, long ny, double *xyz, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yPlm3dData *data;
  long i, size, nval;
  float *dxyz, *dcolr;
  
  /* there are nx-by-ny coordinates and one color */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawPlm3d;
  
  /* compute the size of the data structure plus all data */
  size= sizeof(yPlm3dData)+sizeof(float)*3*nx*ny+sizeof(float)*3;
  data= (yPlm3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nx= nx;
  data->ny= ny;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(yPlm3dData));
  data->colr= dcolr= data->xyz+3*nx*ny;
  /* copy data into the new storage space */
  nval= 3;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  nval= 3*nx*ny;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, nx*ny, data->xyz);
  ygl_fpemask(1);
}

void yglPlf3d(long nx, long ny, double *xyz, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yPlf3dData *data;
  long i, size, nval;
  float *dxyz, *dcolr;
  
  /* there are nx-by-ny coordinates and (nx-1)*(ny-1) colors ???????? */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawPlf3d;

  /* compute the size of the data structure plus all data */
  size= sizeof(yPlf3dData)+sizeof(float)*3*nx*ny+sizeof(float)*3*(nx-1)*(ny-1);
  data= (yPlf3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nx= nx;
  data->ny= ny;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(yPlf3dData));
  data->colr= dcolr= data->xyz+3*nx*ny;
  /* copy data into the new storage space */
  nval= 4*(nx-1)*(ny-1);
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  nval= 3*nx*ny;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, nx*ny, data->xyz);
  ygl_fpemask(1);
}

void yglSurf3d(long do_alpha, long nx, long ny, double *xyz, 
                double *norm, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  ySurf3dData *data;
  long i, size, nval;
  float *dxyz, *dnorm, *dcolr;
  
  /* there are nx-by-ny coordinates and normals and one color ????????? */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawSurf3d;
  
  /* compute the size of the data structure plus all data */
  size= sizeof(ySurf3dData)+sizeof(float)*2*3*nx*ny+sizeof(float)*3;
  data= (ySurf3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->do_alpha= do_alpha;
  data->nx= nx;
  data->ny= ny;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(ySurf3dData));
  data->norm= dnorm= data->xyz+3*nx*ny;
  data->colr= dcolr= data->norm+3*nx*ny;
  /* copy data into the new storage space */
  nval= 3*nx*ny;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
    dnorm[i]= (float) norm[i];
  }
  dcolr[0]= (float) colr[0];
  dcolr[1]= (float) colr[1];
  dcolr[2]= (float) colr[2];
  /* compute and save the bounding box */
  yglSetLims3d(elem, nx*ny, data->xyz);
  ygl_fpemask(1);
}

void yglColrsurf3d(long do_alpha, long nx, long ny, double *xyz, 
                    double *norm, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  ySurf3dData *data;
  long i, size, nval;
  float *dxyz, *dnorm, *dcolr;
  
  /* there are nx-by-ny coordinates and normals and one color ????????? */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawColrSurf3d;
  
  /* compute the size of the data structure plus all data */
  size= sizeof(ySurf3dData)+sizeof(float)*2*3*nx*ny+sizeof(float)*3*nx*ny;
  data= (ySurf3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->do_alpha= do_alpha;
  data->nx= nx;
  data->ny= ny;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(ySurf3dData));
  data->norm= dnorm= data->xyz+3*nx*ny;
  data->colr= dcolr= data->norm+3*nx*ny;
  /* copy data into the new storage space */
  nval= 3*nx*ny;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
    dnorm[i]= (float) norm[i];
  }
  if(do_alpha) nval= 4*nx*ny;
  else nval= 3*nx*ny;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, nx*ny, data->xyz);
  ygl_fpemask(1);
}

void yglLines3d(long nvert, double *xyz, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yLines3dData *data;
  long i, size, nval;
  float *dxyz, *dcolr;
  
  /* there are nvert coordinates and one color */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawLines3d;
  
  /* compute the size of the data structure plus all data */
  size= sizeof(yLines3dData)+sizeof(float)*3*nvert+sizeof(float)*3;
  data= (yLines3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nvert= nvert;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(yLines3dData));
  data->colr= dcolr= data->xyz+3*nvert;
  /* copy data into the new storage space */
  nval= 3;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  nval= 3*nvert;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, nvert, data->xyz);
  ygl_fpemask(1);
}

void yglPoints3d(long nvert, double *xyz, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yPoints3dData *data;
  long i, size, nval;
  float *dxyz, *dcolr;
  
  /* there are nvert coordinates and nvert colors */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawPoints3d;
  
  /* compute the size of the data structure plus all data */
  size= sizeof(yPoints3dData)+sizeof(float)*3*nvert+sizeof(float)*3*nvert;
  data= (yPoints3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nvert= nvert;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(yPoints3dData));
  data->colr= dcolr= data->xyz+3*nvert;
  /* copy data into the new storage space */
  nval= 3*nvert;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  nval= 3*nvert;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, nvert, data->xyz);
  ygl_fpemask(1);
}

void yglTstrips3d(long nstrips, long *len, double *xyz,
                   double *norm, double *colr, long edge,
                   long smooth, long do_light, long do_alpha)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yTstrips3dData *data;
  long size, i, ntri, nvert, colrsiz, nval, nrmsiz, *dlen;
  float *dxyz, *dnorm, *dcolr;
  
  /* The number of triangles is the sum of the nstrips elements of len[].
     There is one color per triangle. ??????
     There are ntri+2*nstrips vertices and normals. ???????
     edge non-zero means outline the polygon, not fill it.
     smooth non-zero means normals have been provided and smooth shading
     should be used.           ?????? smooth only ???????
     do_light means enable lighting ????? maybe means emissive lighting if zero ?????
     do_alpha means colors contain an alpha component
  */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawTstrips3d;
  if(do_alpha) {
    colrsiz= 4;
  } else {
    colrsiz= 3;
  }
  nvert= 0;
  for(i= 0; i < nstrips; i++) {
    nvert += len[i];
  }
  ntri= nvert-2*nstrips;
  /* compute the size of the data structure plus all data */
  size= sizeof(yTstrips3dData)+sizeof(float)*3*nvert+sizeof(float)*colrsiz*ntri;
  size += sizeof(long)*nstrips;
  if(smooth) {
    nrmsiz= 3*nvert;
  } else {
    if(do_light) {
      nrmsiz= 3*ntri;
	} else {
      /* no normals supplied, so don't allocate storage for them. */
      nrmsiz= 0;
	}
  }
  size += sizeof(float)*nrmsiz;
  data= (yTstrips3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nstrips= nstrips;
  data->edge= edge;
  data->smooth= smooth;
  data->do_light= do_light;
  data->do_alpha= do_alpha;
  data->len= dlen= (long *) ((char *)data+sizeof(yTstrips3dData));
  data->xyz= dxyz= (float *) (data->len+nstrips);
  data->norm= dnorm= data->xyz+3*nvert;
  data->colr= dcolr= data->norm+nrmsiz;
  /* copy data into the new storage space */
  for(i= 0; i < nstrips; i++) {
    dlen[i]= len[i];
  }
  nval= colrsiz*ntri;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  if(smooth) {
    nval= 3*nvert;
    for(i= 0; i < nval; i++) {
      dxyz[i]= (float) xyz[i];
      dnorm[i]= (float) norm[i];
    }
  } else {
    nval= 3*nvert;
    for(i= 0; i < nval; i++) {
      dxyz[i]= (float) xyz[i];
    }
    if(do_light) {
      nval= 3*ntri;
      for(i= 0; i < nval; i++) {
        dnorm[i]= (float) norm[i];
      }
    }
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, nvert, data->xyz);
  ygl_fpemask(1);
}

void yglQstrips3d(long nstrips, long *len, double *xyz,
                   double *norm, double *colr, long edge,
                   long smooth, long do_light, long do_alpha)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yQstrips3dData *data;
  long size, i, nquad, nvert, ncoord, colrsiz, nval, nrmsiz, *dlen;
  float *dxyz, *dnorm, *dcolr;
  
  /* The number of quads is the sum of the nstrips elements of len[].
     There is one color per quad. ??????
     There are 2*nquad+2*nstrips vertices and normals. ???????
     edge non-zero means outline the polygon, not fill it.
     smooth non-zero means normals have been provided and smooth shading
     should be used.           ?????? smooth only ???????
     do_light means enable lighting ????? maybe means emissive lighting if zero ?????
     do_alpha means colors contain an alpha component
  */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawQstrips3d;
  if(do_alpha) {
    colrsiz= 4;
  } else {
    colrsiz= 3;
  }
  /* first compute the number of vertices down one edge of the strip */
  nvert= 0;
  for(i= 0; i < nstrips; i++) {
    nvert += len[i];
  }
  nquad= nvert-nstrips;
  /* compute the number of vertices in the strip on both edges */
  ncoord= 2*nvert;
  /* compute the size of the data structure plus all data */
  size= sizeof(yQstrips3dData)+sizeof(float)*3*ncoord+sizeof(float)*colrsiz*nquad;
  size += sizeof(long)*nstrips;
  if(smooth) {
    nrmsiz= 3*ncoord;
  } else {
    if(do_light) {
      nrmsiz= 3*nquad;
    } else {
      nrmsiz= 0;
    }
  }
  size += sizeof(float)*nrmsiz;
  data= (yQstrips3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nstrips= nstrips;
  data->edge= edge;
  data->smooth= smooth;
  data->do_light= do_light;
  data->do_alpha= do_alpha;
  data->len= dlen= (long *) ((char *)data+sizeof(yQstrips3dData));
  data->xyz= dxyz= (float *) (data->len+nstrips);
  data->norm= dnorm= data->xyz+3*ncoord;
  data->colr= dcolr= data->norm+nrmsiz;
  /* copy data into the new storage space */
  for(i= 0; i < nstrips; i++) {
    dlen[i]= len[i];
  }
  nval= colrsiz*nquad;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  if(smooth) {
    nval= 3*ncoord;
    for(i= 0; i < nval; i++) {
      dxyz[i]= (float) xyz[i];
      dnorm[i]= (float) norm[i];
    }
  } else {
    nval= 3*ncoord;
    for(i= 0; i < nval; i++) {
      dxyz[i]= (float) xyz[i];
    }
    if(do_light) {
      for(i= 0; i < nrmsiz; i++) {
        dnorm[i]= (float) norm[i];
      }
    }
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, ncoord, data->xyz);
  ygl_fpemask(1);
}

void yglTstripsndx3d(long nstrips, long numedg, long ntri,
                      long *len, long *ndx, double *xyz,
                      double *norm, double *colr, long edge, long do_alpha)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yTstripsNdx3dData *data;
  long size, i, nvert, colrsiz, nval, *dlen, *dndx;
  float *dxyz, *dnorm, *dcolr;
  
  /* The number of triangles is the sum of the nstrips elements of len[].
     There is one color per triangle. ??????
     There are ntri+2*nstrips vertices and normals. ???????
     edge non-zero means outline the polygon, not fill it.
     do_alpha means colors contain an alpha component
  */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawTstripsNdx3d;
  if(do_alpha) {
    colrsiz= 4;
  } else {
    colrsiz= 3;
  }
  nvert= 0;
  for(i= 0; i < nstrips; i++) {
    nvert += len[i];
  }
  ntri= nvert-2*nstrips;
  /* compute the size of the data structure plus all data */
  size= sizeof(yTstripsNdx3dData)+sizeof(float)*2*3*numedg+sizeof(float)*colrsiz*ntri;
  size += sizeof(long)*nvert+sizeof(long)*nstrips;
  data= (yTstripsNdx3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nstrips= nstrips;
  data->numedg= numedg;
  data->nvert= nvert;
  data->ntri= ntri;
  data->edge= edge;
  data->do_alpha= do_alpha;
  data->len= dlen= (long *) ((char *)data+sizeof(yTstripsNdx3dData));
  data->ndx= dndx= data->len+nstrips;
  data->xyz= dxyz= (float *) (data->ndx+nvert);
  data->norm= dnorm= data->xyz+3*numedg;
  data->colr= dcolr= data->norm+3*numedg;
  /* copy data into the new storage space */
  for(i= 0; i < nstrips; i++) {
    dlen[i]= len[i];
  }
  for(i= 0; i < nvert; i++) {
    dndx[i]= ndx[i];
  }
  for(i= 0; i < 3*numedg; i++) {
    dxyz[i]= (float) xyz[i];
    dnorm[i]= (float) norm[i];
  }
  nval= colrsiz*ntri;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  /* compute and save the bounding box (NOTE: asumes everything in
     xyz is actually used) */
  yglSetLims3d(elem, numedg, data->xyz);
  ygl_fpemask(1);
}

void yglTarray3d(long ntri, double *xyz, double *norm, double *colr, 
                  long edge, long smooth, long do_light, long do_alpha,
                  long cpervrt, long cubemap, long emit)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yTarray3dData *data;
  long size, colrsiz, ncolr, i, nval;
  float *dxyz, *dnorm, *dcolr;
  
  /* There are ntri triangles.
     There is one color per triangle or one color per triangle vertex.
     There are 3*ntri vertices and normals.
     edge non-zero means outline the polygon, not fill it.
     smooth non-zero means normals have been provided and smooth shading
     should be used.           ?????? smooth only ???????
     do_light means enable lighting ????? maybe means emissive lighting if zero ?????
     do_alpha means colors contain an alpha component
     colr_per_vert non-zero means one color per vertex
  */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawTarray3d;
  if(do_alpha) {
    colrsiz= 4;
  } else {
    colrsiz= 3;
  }
  if(cpervrt) {
    ncolr= 3*ntri;
  } else {
    ncolr= ntri;
  }
  /* compute the size of the data structure plus all data */
  size= sizeof(yTarray3dData)+sizeof(float)*9*2*ntri+sizeof(float)*colrsiz*ncolr;
  data= (yTarray3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->ntri= ntri;
  data->edge= edge;
  data->smooth= smooth;
  data->do_light= do_light;
  data->do_alpha= do_alpha;
  data->cpervrt= cpervrt;
  data->cubemap= cubemap;
  data->emit= emit;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(yTarray3dData));
  data->norm= dnorm= data->xyz+9*ntri;
  data->colr= dcolr= data->norm+9*ntri;
  /* copy data into the new storage space */
  nval= colrsiz*ncolr;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  nval= 3*3*ntri;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
    dnorm[i]= (float) norm[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, 3*ntri, data->xyz);
  ygl_fpemask(1);
}

void yglQarray3d(long nquad, double *xyz, double *norm, double *colr, 
                  long edge, long smooth, long do_light, long do_alpha,
                  long cpervrt)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yQarray3dData *data;
  long size, colrsiz, ncolr, i, nval;
  float *dxyz, *dnorm, *dcolr;
  
  /* There are nquad quadrangless.
     There is one color per quad or one color per quad vertex.
     There are 4*nquad vertices and normals.
     edge non-zero means outline the polygon, not fill it.
     smooth non-zero means normals have been provided and smooth shading
     should be used.           ?????? smooth only ???????
     do_light means enable lighting ????? maybe means emissive lighting if zero ?????
     do_alpha means colors contain an alpha component
     colr_per_vert non-zero means one color per vertex
  */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawQarray3d;
  if(do_alpha) {
    colrsiz= 4;
  } else {
    colrsiz= 3;
  }
  if(cpervrt) {
    ncolr= 4*nquad;
  } else {
    ncolr= nquad;
  }
  /* compute the size of the data structure plus all data */
  size= sizeof(yQarray3dData)+sizeof(float)*12*2*nquad+sizeof(float)*colrsiz*ncolr;
  data= (yQarray3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nquad= nquad;
  data->edge= edge;
  data->smooth= smooth;
  data->do_light= do_light;
  data->do_alpha= do_alpha;
  data->cpervrt= cpervrt;
  data->xyz= dxyz= (float *) ((char *)data+sizeof(yQarray3dData));
  data->norm= dnorm= data->xyz+12*nquad;
  data->colr= dcolr= data->norm+12*nquad;
  /* copy data into the new storage space */
  nval= colrsiz*ncolr;
  for(i= 0; i < nval; i++) {
    dcolr[i]= (float) colr[i];
  }
  nval= 3*4*nquad;
  for(i= 0; i < nval; i++) {
    dxyz[i]= (float) xyz[i];
    dnorm[i]= (float) norm[i];
  }
  /* compute and save the bounding box */
  yglSetLims3d(elem, 4*nquad, data->xyz);
  ygl_fpemask(1);
}

void yglTivarray3d(long ntri, long nvert, long *ptndx, double *xyz, 
                    double *norm, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yTivarray3dData *data;
  long size, i, numndx;
  unsigned int *dptndx;
  float *dileave;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  
  /* There are ntri triangles and nvert vertices.
     Each triangle is specified by three indices into the 
     xyz, norm, and colr arrays.
     The colors are RGBA.
  */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawTivarray3d;
  numndx= 3*ntri;
  /* compute the size of the data structure plus all data */
  size= sizeof(yTivarray3dData)+sizeof(float)*(3+3+4)*nvert+sizeof(dptndx[0])*numndx;
  data= (yTivarray3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->ntri= ntri;
  data->nvert= nvert;
  data->ptndx= dptndx= (unsigned int *) ((char *)data+sizeof(yTivarray3dData));
  data->ileave= dileave= (float *) (data->ptndx+numndx);
  /* copy data into the new storage space */
  for(i= 0; i < numndx; i++) {
    dptndx[i]= (unsigned int) ptndx[i];
  }
  for(i= 0; i < nvert; i++) {
    dileave[10*i  ]= (float) colr[4*i];
    dileave[10*i+1]= (float) colr[4*i+1];
    dileave[10*i+2]= (float) colr[4*i+2];
    dileave[10*i+3]= (float) colr[4*i+3];
    dileave[10*i+4]= (float) norm[3*i];
    dileave[10*i+5]= (float) norm[3*i+1];
    dileave[10*i+6]= (float) norm[3*i+2];
    dileave[10*i+7]= (float) xyz[3*i];
    dileave[10*i+8]= (float) xyz[3*i+1];
    dileave[10*i+9]= (float) xyz[3*i+2];
  }
  /* Compute and save the bounding box.
     WARNING: this could be WRONG if there are points in the arrays
     that are never referenced. */
  if(nvert > 0) {
    xmin= xmax= xyz[0];
    ymin= ymax= xyz[1];
    zmin= zmax= xyz[2];
    for(i= 1; i < nvert; i++) {
      if(xyz[3*i] < xmin)   xmin= xyz[3*i];
      if(xyz[3*i] > xmax)   xmax= xyz[3*i];
      if(xyz[3*i+1] < ymin) ymin= xyz[3*i+1];
      if(xyz[3*i+1] > ymax) ymax= xyz[3*i+1];
      if(xyz[3*i+2] < zmin) zmin= xyz[3*i+2];
      if(xyz[3*i+2] > zmax) zmax= xyz[3*i+2];
    }
    elem->box.xmin= (float) xmin;
    elem->box.xmax= (float) xmax;
    elem->box.ymin= (float) ymin;
    elem->box.ymax= (float) ymax;
    elem->box.zmin= (float) zmin;
    elem->box.zmax= (float) zmax;
  }
  ygl_fpemask(1);
}

void yglTvarray3d(long ntri, long nvert, long do_alpha, long cpervrt, long *ptndx,
                  double *xyz, double *norm, double *colr)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yTvarray3dData *data;
  long size, i, ncopy, clrsiz;
  unsigned int *dptndx;
  float *dxyz, *dnorm, *dcolr;
  double xmin, xmax, ymin, ymax, zmin, zmax;
#undef TVDBG
#ifdef TVDBG
  long ndmin=2000000000, ndmax=-2000000000;
  double nrmin=1.0e150, nrmax=-1.0e150, ntmp;
  double nrmin2=1.0e150, nrmax2=-1.0e150;
  xmin=1.0e150; xmax=-1.0e150;
#endif
  
  /* There are ntri triangles and nvert vertices.
     Each triangle is specified by three indices into the 
     xyz, norm, and colr arrays. An option is also provided to use a single
     color for all triangles.
     The colors are RGB or RGBA.
     This is like yglTivarray3d except it uses three separate arrays
     instead of one interleaved array.
  */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawTvarray3d;
  /* compute the size of the data structure plus all data */
  size= sizeof(yTvarray3dData)+sizeof(float)*(3+3)*nvert+sizeof(dptndx[0])*3*ntri;
  if(do_alpha) clrsiz= 4;
  else clrsiz= 3;
  if(cpervrt) {
    size += sizeof(float)*clrsiz*nvert;
  } else {
    size += sizeof(float)*clrsiz;
  }
  data= (yTvarray3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->ntri= ntri;
  data->nvert= nvert;
  data->cpervrt= cpervrt;
  data->do_alpha= do_alpha;
  data->ptndx= dptndx= (unsigned int *) ((char *)data+sizeof(yTvarray3dData));
  data->xyz= dxyz= (float *) (data->ptndx+3*ntri);
  data->norm= dnorm= data->xyz+3*nvert;
  data->colr= dcolr= data->norm+3*nvert;
#ifdef TVDBG
  printf("xyz pointer is %x, norm pointer is %x, colr pointer is %x\n", data->xyz, data->norm, data->colr);
#endif
  /* copy data into the new storage space */
  for(i= 0; i < 3*ntri; i++) {
#ifdef TVDBG
    if(ptndx[i] < ndmin) ndmin= ptndx[i];
    if(ptndx[i] > ndmax) ndmax= ptndx[i];
#endif
    dptndx[i]= (unsigned int) ptndx[i];
  }
#ifdef TVDBG
/*
  printf("elements for xyz is %ld\n", data->norm-data->xyz);
  printf("elements for norm is %ld\n", data->colr-data->norm);
    printf("struct uses %ld bytes, %ld bytes allocated\n", (char *)(data->colr+clrsiz*nvert)-(char *)data, size);
*/
  printf("minimum index is %ld and max is %ld, nvert is %ld\n", ndmin, ndmax, nvert);
#endif
  /* there are the same number of coords and normals in all cases */
  for(i= 0; i < nvert; i++) {
#ifdef TVDBG
    if(xyz[3*i] < xmin) xmin= xyz[3*i];
    if(xyz[3*i] > xmax) xmax= xyz[3*i];
    if(xyz[3*i+1] < xmin) xmin= xyz[3*i+1];
    if(xyz[3*i+1] > xmax) xmax= xyz[3*i+1];
    if(xyz[3*i+2] < xmin) xmin= xyz[3*i+2];
    if(xyz[3*i+2] > xmax) xmax= xyz[3*i+2];
    ntmp= norm[3*i]*norm[3*i]+norm[3*i+1]*norm[3*i+1]+norm[3*i+2]*norm[3*i+2];
    if(ntmp < nrmin) nrmin= ntmp;
    if(ntmp > nrmax) nrmax= ntmp;
#endif
    dxyz[3*i  ]= (float) xyz[3*i];
    dxyz[3*i+1]= (float) xyz[3*i+1];
    dxyz[3*i+2]= (float) xyz[3*i+2];
    dnorm[3*i  ]= (float) norm[3*i];
    dnorm[3*i+1]= (float) norm[3*i+1];
    dnorm[3*i+2]= (float) norm[3*i+2];
#ifdef TVDBG
    ntmp= dnorm[3*i]*dnorm[3*i]+dnorm[3*i+1]*dnorm[3*i+1]+dnorm[3*i+2]*dnorm[3*i+2];
    if(ntmp < nrmin2) nrmin2= ntmp;
    if(ntmp > nrmax2) nrmax2= ntmp;
#endif
  }
#ifdef TVDBG
  printf("minimum x is %e and max is %e\n", xmin, xmax);
  printf("minimum norm is %e and max is %e\n", nrmin, nrmax);
  printf("minimum norm in struct is %e and max is %e\n", nrmin2, nrmax2);
#endif
  if(cpervrt) {
    if(do_alpha) {
      ncopy= 4*nvert;
    } else {
      ncopy= 3*nvert;
    }
  } else {
    if(do_alpha) {
      ncopy= 4;
    } else {
      ncopy= 3;
    }
  }
  for(i= 0; i < ncopy; i++) {
    dcolr[i]= (float) colr[i];
  }
  /* Compute and save the bounding box.
     WARNING: this could be WRONG if there are points in the arrays
     that are never referenced. */
  if(nvert > 0) {
    xmin= xmax= xyz[0];
    ymin= ymax= xyz[1];
    zmin= zmax= xyz[2];
    for(i= 1; i < nvert; i++) {
      if(xyz[3*i] < xmin)   xmin= xyz[3*i];
      if(xyz[3*i] > xmax)   xmax= xyz[3*i];
      if(xyz[3*i+1] < ymin) ymin= xyz[3*i+1];
      if(xyz[3*i+1] > ymax) ymax= xyz[3*i+1];
      if(xyz[3*i+2] < zmin) zmin= xyz[3*i+2];
      if(xyz[3*i+2] > zmax) zmax= xyz[3*i+2];
    }
    elem->box.xmin= (float) xmin;
    elem->box.xmax= (float) xmax;
    elem->box.ymin= (float) ymin;
    elem->box.ymax= (float) ymax;
    elem->box.zmin= (float) zmin;
    elem->box.zmax= (float) zmax;
  }
  ygl_fpemask(1);
}

void yglTex3d(float ds, double *origin, double *boxsiz)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yTex3dData *data;
  long i, size, nval;
  double *dorigin, *dboxsiz;
  
  /* there are nx-by-ny coordinates and normals and one color ????????? */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  ygl_fpemask(0);
  elem= yglNewDirectList3dElem();
  elem->func= yglDrawTex3d;
  
  /* compute the size of the data structure plus all data */
  size= sizeof(yTex3dData)+sizeof(double)*2*3;
  data= (yTex3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->ds= ds;
  data->origin= dorigin= (double *) ((char *)data+sizeof(yTex3dData));
  data->boxsiz= dboxsiz= data->origin+3;
  /* copy data into the new storage space */
  nval= 3;
  for(i= 0; i < nval; i++) {
    dorigin[i]= origin[i];
    dboxsiz[i]=  boxsiz[i];
  }
  /* compute and save the bounding box */
  elem->box.xmin= origin[0];
  elem->box.xmax= origin[0]+boxsiz[0];
  elem->box.ymin= origin[1];
  elem->box.ymax= origin[1]+boxsiz[1];
  elem->box.zmin= origin[2];
  elem->box.zmax= origin[2]+boxsiz[2];
  ygl_fpemask(1);
}

void yglTexcell2d(long nx, long ny, long nz, double *znsiz, char *texval)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yTexcell2dData *data;
  long i, size, nval;
  unsigned char *dtexval;
  double *dznsiz;
  
  /* there are nx-by-ny coordinates and normals and one color ????????? */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawTexcell2d;
  
  /* compute the size of the data structure plus all data */
  size= sizeof(yTexcell2dData)+sizeof(double)*3+sizeof(char)*4*nx*ny*nz;
  data= (yTexcell2dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nx= nx;
  data->ny= ny;
  data->nz= nz;
  data->znsiz= dznsiz= (double *) ((char *)data+sizeof(yTexcell2dData));
  data->texval= dtexval= (unsigned char *) (data->znsiz+3);
  /* copy data into the new storage space */
  nval= 3;
  for(i= 0; i < nval; i++) {
    dznsiz[i]= znsiz[i];
  }
  nval= 4*nx*ny*nz;
  for(i= 0; i < nval; i++) {
    dtexval[i]= (unsigned char) texval[i];
  }
  /* compute and save the bounding box */
  elem->box.xmin= 0.0;
  elem->box.xmax= (nx-1)*znsiz[0];
  elem->box.ymin= 0.0;
  elem->box.ymax= (ny-1)*znsiz[1];
  elem->box.zmin= 0.0;
  elem->box.zmax= (nz-1)*znsiz[2];
  ygl_fpemask(1);
}

void yglPlpix3d(long nx, long ny, char *pix)
{
  /* variables for the display list */
  yList3d_Elem *elem;
  yPix3dData *data;
  long i, size, nval;
  unsigned char *dpix;
  
  /* there are nx-by-ny pixels */
  /* add a new entry to the 3D display list (it will be in the list
     when it is returned). */
  if(!glCurrWin3d) return;
  ygl_fpemask(0);
  if(glCurrWin3d->use_list) {
    elem= yglNewCachedList3dElem();
  } else {
    elem= yglNewDirectList3dElem();
  }
  elem->func= yglDrawPix3d;
  /* compute the size of the data structure plus all data */
  size= sizeof(yPix3dData)+sizeof(char)*3*nx*ny;
  data= (yPix3dData *) p_malloc(size);
  elem->data= data;
  /* set all elements of the data structure, including pointers into the data block */
  data->nx= nx;
  data->ny= ny;
  data->pix= dpix= (unsigned char *)data+sizeof(yPix3dData);
  /* copy data into the new storage space */
  nval= 3*nx*ny;
  for(i= 0; i < nval; i++) {
    dpix[i]= (char) pix[i];
  }
  /* A bounding box doesn't make sense for this function because
     it uses pixel coords. */
  ygl_fpemask(1);
}

void yglDrawPolys3d(int mode, void *vdata)
{
  yPoly3dData *data= (yPoly3dData *)vdata;

  yglPolys(data->npolys, data->len, data->xyz, data->norm, data->colr, 
           data->edge, data->smooth, data->do_light);
}

void yglDrawGlyphs3d(int mode, void *vdata)
{
  yGlyph3dData *data= (yGlyph3dData *)vdata;

  yglGlyphs(data->nglyph, data->origin, data->scal, data->theta, 
            data->phi, data->colr);
}

void yglDrawCells3d(int mode, void *vdata)
{
  yCell3dData *data= (yCell3dData *)vdata;

  yglCells(data->nx, data->ny, (float (*)[3])(data->xyz), data->norm, 
           data->colr, data->do_alpha);
}

void yglDrawPlm3d(int mode, void *vdata)
{
  yPlm3dData *data= (yPlm3dData *)vdata;

  yglPlm(data->nx, data->ny, data->xyz, data->colr);
}

void yglDrawPlf3d(int mode, void *vdata)
{
  yPlf3dData *data= (yPlf3dData *)vdata;

  yglPlf(data->nx, data->ny, data->xyz, data->colr);
}

void yglDrawSurf3d(int mode, void *vdata)
{
  ySurf3dData *data= (ySurf3dData *)vdata;

  yglSurf(data->do_alpha, data->nx, data->ny, data->xyz, data->norm, data->colr);
}

void yglDrawColrSurf3d(int mode, void *vdata)
{
  yColrSurf3dData *data= (yColrSurf3dData *)vdata;

  yglColrSurf(data->do_alpha, data->nx, data->ny, data->xyz, data->norm, data->colr);
}

void yglDrawLines3d(int mode, void *vdata)
{
  yLines3dData *data= (yLines3dData *)vdata;

  yglLines(data->nvert, data->xyz, data->colr);
}

void yglDrawPoints3d(int mode, void *vdata)
{
  yPoints3dData *data= (yPoints3dData *)vdata;

  yglPoints(data->nvert, data->xyz, data->colr);
}

void yglDrawTstrips3d(int mode, void *vdata)
{
  yTstrips3dData *data= (yTstrips3dData *)vdata;

  if(data->do_alpha) {
    yglTstripsAlpha(data->nstrips, data->len, data->xyz, data->norm, 
             data->colr, data->edge, data->smooth, data->do_light);
  } else {
    yglTstrips(data->nstrips, data->len, data->xyz, data->norm, 
             data->colr, data->edge, data->smooth, data->do_light);
  }
}

void yglDrawTstripsNdx3d(int mode, void *vdata)
{
  yTstripsNdx3dData *data= (yTstripsNdx3dData *)vdata;

  if(data->do_alpha) {
    yglTstripsAlphaNdx(data->nstrips, data->numedg, data->ntri, data->len, 
                       data->ndx, data->xyz, data->norm, data->colr, data->edge);
  } else {
    yglTstripsNdx(data->nstrips, data->numedg, data->ntri, data->len, 
                  data->ndx, data->xyz, data->norm, data->colr, data->edge);
  }
}

void yglDrawQstrips3d(int mode, void *vdata)
{
  yQstrips3dData *data= (yQstrips3dData *)vdata;

  if(data->do_alpha) {
    yglQstripsAlpha(data->nstrips, data->len, data->xyz, data->norm, 
             data->colr, data->edge, data->smooth, data->do_light);
  } else {
    yglQstrips(data->nstrips, data->len, data->xyz, data->norm, 
             data->colr, data->edge, data->smooth, data->do_light);
  }
}

void yglDrawTarray3d(int mode, void *vdata)
{
  yTarray3dData *data= (yTarray3dData *)vdata;

  if(data->do_alpha) {
    yglTarrayAlpha(data->smooth, data->ntri, data->xyz, data->norm, data->colr, 
             data->edge, data->cpervrt, data->emit);
  } else {
    yglTarray(data->smooth, data->ntri, data->xyz, data->norm, data->colr, 
             data->edge, data->cpervrt, data->emit);
  }
}

void yglDrawQarray3d(int mode, void *vdata)
{
  yQarray3dData *data= (yQarray3dData *)vdata;

  if(data->do_alpha) {
    yglQarrayAlpha(data->smooth, data->nquad, data->xyz, data->norm, data->colr, 
             data->edge, data->cpervrt);
  } else {
    yglQarray(data->smooth, data->nquad, data->xyz, data->norm, data->colr, 
             data->edge, data->cpervrt);
  }
}

void yglDrawTivarray3d(int mode, void *vdata)
{
  yTivarray3dData *data= (yTivarray3dData *)vdata;

  yglTivarray(data->ntri, data->ptndx, data->ileave);
}

void yglDrawTvarray3d(int mode, void *vdata)
{
  yTvarray3dData *data= (yTvarray3dData *)vdata;

  yglTvarray(data->do_alpha, data->cpervrt, data->ntri, data->ptndx, data->xyz,
             data->norm, data->colr);
}

void yglDrawTex3d(int mode, void *vdata)
{
  yTex3dData *data= (yTex3dData *)vdata;

  yglTex3dbox(data->ds, data->origin, data->boxsiz);
}

void yglDrawTexcell2d(int mode, void *vdata)
{
  yTexcell2dData *data= (yTexcell2dData *)vdata;

  yglTexcell2(data->nx, data->ny, data->nz, data->znsiz, data->texval);
}

void yglDrawPix3d(int mode, void *vdata)
{
  yPix3dData *data= (yPix3dData *)vdata;

  yglPutPixels(data->nx, data->ny, data->pix);
}
