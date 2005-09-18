/*
 * $Id: TriUtil.c,v 1.1 2005-09-18 22:07:41 dhmunro Exp $
 * This file contains utility functions for manipulating lists
 * of triangle arrays and strips
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include <stdlib.h>
#include <math.h>
#include "TriStruct.h"
#include "glWrappers.h"
#include "glcode.h"
#include "glfunc.h"

typedef struct entry entry;
struct entry {
  double val;
  long ndx;
} ;

extern void myqsort(entry v[], long left, long right);

long yglSizeTriArraysNdx3d(TriVertexGrp *list)
{
  long numtri;

  /* 
	 This function returns the number of triangles in a list of triangle arrays.
  */
  numtri= 0;
  while(list) {
	numtri += list->numTri;
	list= list->next;
  }
  /* Return the number of triangles summed over everything in the list. */
  return numtri;
}

void yglCopyTriArrayNdx3d( TriVertexGrp *list, TriVertexGrp *nlist)
{
  long i, numtri, numpt, *oldIDs, *newIDs, *oldNdx, *newNdx;
  yPoint3D *oldVert, *newVert, *oldNorm, *newNorm;
  double *oldVar2, *newVar2;

  /* 
	 This function copies a triangle array into another 
     triangle array.
	 WARNING!!! it does not copy coolors!!!!!.
     The caller guarantees that there is enough room in
     nlist to store all the triangles in list 
  */
  newVert= nlist->xyzverts;
  newNorm= nlist->normals;
  newIDs=  nlist->cellIDs;
  newVar2= nlist->var2;
  newNdx=  nlist->ptndx;
  oldVert= list->xyzverts;
  oldNorm= list->normals;
  oldIDs=  list->cellIDs;
  oldVar2= list->var2;
  oldNdx=  list->ptndx;
  numtri=  list->numTri;
  numpt=   list->numEdg;
  for(i=0; i < numtri; i++) {
	newNdx[0]= oldNdx[0];
	newNdx[1]= oldNdx[1];
	newNdx[2]= oldNdx[2];
	newNdx += 3;
	oldNdx += 3;
  }
  for(i=0; i < numpt; i++) {
	newVert[0]= oldVert[0];
	newVert[1]= oldVert[1];
	newVert[2]= oldVert[2];
	newVert += 3;
	oldVert += 3;
	newNorm[0]= oldNorm[0];
	newNorm[1]= oldNorm[1];
	newNorm[2]= oldNorm[2];
	newNorm += 3;
	oldNorm += 3;
	if(oldVar2) {
      newVar2[0] = oldVar2[0];
      newVar2[1] = oldVar2[1];
      newVar2[2] = oldVar2[2];
	  oldVar2 += 3;
	  newVar2 += 3;
	}
	*newIDs++  = *oldIDs++;
  }
}

void yglCollapseTriArraysNdx3d(long colrtyp, TriVertexGrp *list, TriVertexGrp *nlist)
{
  long i, count, nedg, comp4, colrskip, numtri, ptbase;
  long *oldIDs, *newIDs, *oldNdx, *newNdx;
  yPoint3D *oldVert, *newVert, *oldNorm, *newNorm;
  float *oldColr, *newColr;
  double *oldVar2, *newVar2;

  /* 
	 This function merges a list of triangle arrays into a 
     single triangle array.
	 abs(colrtyp) is 3 for RGB colors and 4 for RGBA.
	 if colrtyp<0, the input list has scalar colors.
     The caller guarantees that there is enough room in
     nlist to store all the triangles in list 
  */
  newVert= nlist->xyzverts;
  newNorm= nlist->normals;
  newColr= nlist->colors;
  newIDs=  nlist->cellIDs;
  newVar2= nlist->var2;
  newNdx=  nlist->ptndx;
  if(colrtyp < 0) {
    colrskip= 0;
	if(colrtyp == -4) comp4= 1;
	else comp4= 0;
  } else {
    colrskip= colrtyp;
	if(colrtyp == 4) comp4= 1;
	else comp4= 0;
  }
  numtri= 0;
  ptbase= 0;
  
  while(list) {
    oldVert= list->xyzverts;
    oldNorm= list->normals;
    oldColr= list->colors;
    oldIDs=  list->cellIDs;
    oldVar2= list->var2;
    oldNdx=  list->ptndx;
	count= list->numTri;
    nedg=   list->numEdg;
	numtri += count;
    for(i=0; i < count; i++) {
      /* the vertices for this tri array will be stored after the ones
         for all previous arrays, so the indices must be adjusted. */
	  newNdx[0]= oldNdx[0]+ptbase;
	  newNdx[1]= oldNdx[1]+ptbase;
	  newNdx[2]= oldNdx[2]+ptbase;
	  newNdx += 3;
	  oldNdx += 3;
	}
    for(i=0; i < nedg; i++) {
	  newVert[0]= oldVert[0];
	  newVert[1]= oldVert[1];
	  newVert[2]= oldVert[2];
	  newVert += 3;
	  oldVert += 3;
	  newNorm[0]= oldNorm[0];
	  newNorm[1]= oldNorm[1];
	  newNorm[2]= oldNorm[2];
	  newNorm += 3;
	  oldNorm += 3;
	  newColr[0]= oldColr[0];
	  newColr[1]= oldColr[1];
	  newColr[2]= oldColr[2];
	  if(comp4) {
        newColr[3]= oldColr[3];
		newColr += 4;
	  } else {
	    newColr += 3;
	  }
	  oldColr += colrskip;
      if(oldVar2) {
        newVar2[0] = oldVar2[0];
        newVar2[1] = oldVar2[1];
        newVar2[2] = oldVar2[2];
	    oldVar2 += 3;
	    newVar2 += 3;
	  }
	  *newIDs++  = *oldIDs++;
	}
    ptbase += nedg;
	list= list->next;
  }
  nlist->numTri= numtri;
  nlist->numEdg= ptbase;
}

void yglDoSortTriNdx3d(TriVertexGrp *oldtri, long *newptndx)
{
  entry *dist;
  double view[3], viewlen;
  yPoint3D *oldxyz;
  long *oldndx, i, ntri, ndx, base, pt0, pt1, pt2;;

  /* 
	 This function performs a depth sort of a triangle
	 array. It uses the current viewpoint.
	 The CALLER is responsible for inserting the new point indices
     into the triangle array (because the input tri array uses
     indices into coord, normal, etc. arrays, only the indices must be changed).
	 The resulting triangle array should display correctly
	 with translucency when viewed with this viewing transform.
	 Rotating the scene may lead to incorrect results
	 if the triangles are translucent.
	 newptndx must be of length at least 3*ntri.
  */
  view[0]= glCurrWin3d->eye[0] - glCurrWin3d->center[0];
  view[1]= glCurrWin3d->eye[1] - glCurrWin3d->center[1];
  view[2]= glCurrWin3d->eye[2] - glCurrWin3d->center[2];
  viewlen= 1.0e-80+sqrt(view[0]*view[0]+view[1]*view[1]+view[2]*view[2]);
  view[0] /= viewlen;
  view[1] /= viewlen;
  view[2] /= viewlen;
  /* compute the average z-depth of the three vertices of each triangle */
  ntri= oldtri->numTri;
  oldxyz= oldtri->xyzverts;
  oldndx= oldtri->ptndx;
  dist= (entry *) malloc(ntri*sizeof(entry));
  for(i= 0; i < ntri; i++) {
    dist[i].ndx= i;
    base= 3*i;
    pt0= oldndx[base];
    pt1= oldndx[base+1];
    pt2= oldndx[base+2];
    /* this uses coords that are on average a factor of 3 too large, 
       but that doesn't affect the sort order */
    dist[i].val= view[0]*(oldxyz[pt0].x+oldxyz[pt1].x+oldxyz[pt2].x)
		    +view[1]*(oldxyz[pt0].y+oldxyz[pt1].y+oldxyz[pt2].y)
			+view[2]*(oldxyz[pt0].z+oldxyz[pt1].z+oldxyz[pt2].z);
  }
  /* save new triangle order that sorts the polygons from smallest z to largest z */
  for(i= 0; i < ntri; i++) {
    ndx= 3*(dist[i].ndx);
	base= 3*i;
	newptndx[base]=   oldndx[ndx];
	newptndx[base+1]= oldndx[ndx+1];
	newptndx[base+2]= oldndx[ndx+2];
  }
  free(dist);
}

long yglSizeTriArrays3d(TriArrayGrp *list)
{
  long numtri;

  /* 
	 This function returns the number of triangles in a list of triangle arrays.
  */
  numtri= 0;
  while(list) {
	numtri += list->numTri;
	list= list->next;
  }
  /* Return the number of triangles summed over everything in the list. */
  return numtri;
}

void yglCopyTriArray3d(long numtri, TriArrayGrp *list, TriArrayGrp *nlist)
{
  long i, *oldIDs, *newIDs;
  yPoint3D *oldVert, *newVert, *oldNorm, *newNorm;
  double *oldVar2, *newVar2;

  /*
	 This function copies a triangle array into another 
     triangle array.
	 WARNING!!! it does not copy coolors!!!!!.
     The caller guarantees that there is enough room in
     nlist to store all the triangles in list 
  */
  newVert= nlist->xyzverts;
  newNorm= nlist->normals;
  newIDs=  nlist->cellIDs;
  newVar2= nlist->var2;
  oldVert= list->xyzverts;
  oldNorm= list->normals;
  oldIDs=  list->cellIDs;
  oldVar2= list->var2;
  for(i=0; i < numtri; i++) {
	newVert[0]= oldVert[0];
	newVert[1]= oldVert[1];
	newVert[2]= oldVert[2];
	newVert += 3;
	oldVert += 3;
	newNorm[0]= oldNorm[0];
	newNorm[1]= oldNorm[1];
	newNorm[2]= oldNorm[2];
	newNorm += 3;
	oldNorm += 3;
	if(oldVar2) {
      newVar2[0] = oldVar2[0];
      newVar2[1] = oldVar2[1];
      newVar2[2] = oldVar2[2];
	  oldVar2 += 3;
	  newVar2 += 3;
	}
	*newIDs++  = *oldIDs++;
  }
}

void yglCollapseTriArrays3d(long colrtyp, TriArrayGrp *list, TriArrayGrp *nlist)
{
  long i, count, comp4, colrskip, numtri, *oldIDs, *newIDs, color_loop, il;
  yPoint3D *oldVert, *newVert, *oldNorm, *newNorm;
  float *oldColr, *newColr;
  double *oldVar2, *newVar2;

  /* 
     This function merges a list of triangle arrays into a 
     single triangle array.
     abs(colrtyp) is 3 for RGB colors and 4 for RGBA.
     if colrtyp<0, the input list has scalar colors.
     The caller guarantees that there is enough room in
     nlist to store all the triangles in list 
  */
  newVert= nlist->xyzverts;
  newNorm= nlist->normals;
  newColr= nlist->colors;
  newIDs=  nlist->cellIDs;
  newVar2= nlist->var2;
  color_loop= 1;
  if(colrtyp < 0) {
    colrskip= 0;
    if(colrtyp < -4) {
      color_loop= 3;
      colrtyp += 16;
    }
    if(colrtyp == -4) comp4= 1;
    else comp4= 0;
  } else {
    if(colrtyp > 4) {
      color_loop= 3;
      colrtyp -= 16;
    }
    colrskip= colrtyp;
    if(colrtyp == 4) comp4= 1;
    else comp4= 0;
  }
  numtri= 0;
  while(list) {
    oldVert= list->xyzverts;
    oldNorm= list->normals;
    oldColr= list->colors;
    oldIDs=  list->cellIDs;
    oldVar2= list->var2;
    count= list->numTri;
    numtri += count;
    for(i=0; i < count; i++) {
      newVert[0]= oldVert[0];
      newVert[1]= oldVert[1];
      newVert[2]= oldVert[2];
      newVert += 3;
      oldVert += 3;
      newNorm[0]= oldNorm[0];
      newNorm[1]= oldNorm[1];
      newNorm[2]= oldNorm[2];
      newNorm += 3;
      oldNorm += 3;
      // handle both color per vertex and color per triangle
      for(il= 0; il < color_loop; il++) {
        newColr[0]= oldColr[0];
        newColr[1]= oldColr[1];
        newColr[2]= oldColr[2];
        if(comp4) {
          newColr[3]= oldColr[3];
          newColr += 4;
        } else {
          newColr += 3;
        }
        oldColr += colrskip;
      }
      if(oldVar2) {
        newVar2[0] = oldVar2[0];
        newVar2[1] = oldVar2[1];
        newVar2[2] = oldVar2[2];
        oldVar2 += 3;
        newVar2 += 3;
      }
      *newIDs++  = *oldIDs++;
    }
    list= list->next;
  }
  nlist->numTri= numtri;
}

void yglDoSortTri3d(long colrtyp, TriArrayGrp *oldtri, TriArrayGrp *newtri)
{
  entry *dist;
  double view[3], viewlen;
  yPoint3D *oldxyz, *newxyz, *oldnorm, *newnorm;
  double *oldVar2, *newVar2;
  float *oldcolr, *newcolr;
  long *oldid, *newid;
  long i, ntri, ndx, base, cbase, xbase;

  /* 
     Perform a depth sort of a triangle
     array. It uses the current viewpoint.
     The resulting triangle array should display correctly
     with translucency when viewed with this viewing transform.
     Rotating the scene may lead to incorrect results
     if the triangles are translucent.
     colrtyp is 3 for RGB colors and 4 for RGBA.
     There MUST be one color for every triangle vertex.
  */
  view[0]= glCurrWin3d->eye[0] - glCurrWin3d->center[0];
  view[1]= glCurrWin3d->eye[1] - glCurrWin3d->center[1];
  view[2]= glCurrWin3d->eye[2] - glCurrWin3d->center[2];
  viewlen= 1.0e-80+sqrt(view[0]*view[0]+view[1]*view[1]+view[2]*view[2]);
  view[0] /= viewlen;
  view[1] /= viewlen;
  view[2] /= viewlen;
  /* compute the average z-depth of the three vertices 
     of each triangle */
  ntri= oldtri->numTri;
  oldxyz= oldtri->xyzverts;
  dist= (entry *) malloc(ntri*sizeof(entry));
  for(i= 0; i < ntri; i++) {
    dist[i].ndx= i;
    base= 3*i;
    dist[i].val= view[0]*(oldxyz[base].x+oldxyz[base+1].x+oldxyz[base+2].x)
		    +view[1]*(oldxyz[base].y+oldxyz[base+1].y+oldxyz[base+2].y)
			+view[2]*(oldxyz[base].z+oldxyz[base+1].z+oldxyz[base+2].z);
  }
  /* sort the polygons from smallest z to largest z */
  myqsort(dist, 0, ntri-1);
  cbase= 0;
  xbase= 0;
  newxyz= newtri->xyzverts;
  oldnorm= oldtri->normals;
  newnorm= newtri->normals;
  oldcolr= oldtri->colors;
  newcolr= newtri->colors;
  oldid= oldtri->cellIDs;
  newid= newtri->cellIDs;
  oldVar2= oldtri->var2;
  newVar2= newtri->var2;
  for(i= 0; i < ntri; i++) {
    ndx= 3*(dist[i].ndx);
	newxyz[xbase]=   oldxyz[ndx];
	newxyz[xbase+1]= oldxyz[ndx+1];
	newxyz[xbase+2]= oldxyz[ndx+2];
	newnorm[xbase]=   oldnorm[ndx];
	newnorm[xbase+1]= oldnorm[ndx+1];
	newnorm[xbase+2]= oldnorm[ndx+2];
	if(oldVar2) {
      newVar2[xbase]=   oldVar2[ndx];
      newVar2[xbase+1]= oldVar2[ndx+1];
      newVar2[xbase+2]= oldVar2[ndx+2];
	}
	xbase += 3;
  }
  if(colrtyp == 4) {
    for(i= 0; i < ntri; i++) {
      ndx= dist[i].ndx;
      newid[i]= oldid[ndx];
	  ndx *= 4;
	  newcolr[cbase]=   oldcolr[ndx];
	  newcolr[cbase+1]= oldcolr[ndx+1];
	  newcolr[cbase+2]= oldcolr[ndx+2];
	  newcolr[cbase+3]= oldcolr[ndx+3];
      cbase += 4;
	}
  } else {
    for(i= 0; i < ntri; i++) {
      ndx= dist[i].ndx;
      newid[i]= oldid[ndx];
	  ndx *= 3;
	  newcolr[cbase]=   oldcolr[ndx];
	  newcolr[cbase+1]= oldcolr[ndx+1];
	  newcolr[cbase+2]= oldcolr[ndx+2];
      cbase += 3;
	}
  }
  free(dist);
}

void myqsort(entry v[], long left, long right)
{
  long i, last, mid;
  entry temp;

  if(left >= right) return;
  mid= (left+right)/2;
  temp= v[left];
  v[left]= v[mid];
  v[mid]= temp;
  last= left;
  for(i= left+1; i <= right; i++) {
    if(v[i].val < v[left].val) {
	  last++;
	  temp= v[i];
	  v[i]= v[last];
	  v[last]= temp;
	}
  }
  temp= v[left];
  v[left]= v[last];
  v[last]= temp;
  myqsort(v, left, last-1);
  myqsort(v, last+1, right);
}

void yglSliceTris3d(long *keep, long *nkeep, double *dp,
                  TriArrayGrp *oldtri, TriArrayGrp *newtri)
{
  long numold, inew, i3new, i, i3, p0, p1, p2, *oldid, *newid;
  yPoint3D *oldxyz, *newxyz, *oldnorm, *newnorm;
  double *oldVar2, *newVar2, frac;
  float *oldcolr, *newcolr;

  /*
    Clip triangles against a plane using nkeep,
    the number of vertices from each input triangle that are "above"
    the slicing plane, keep (above the plane for a specific vertex),
    and dp (the distance to the plane for each vertex). 
  */
  numold=  oldtri->numTri;
  oldxyz=  oldtri->xyzverts;
  oldnorm= oldtri->normals;
  oldcolr= oldtri->colors;
  oldid=   oldtri->cellIDs;
  oldVar2= oldtri->var2;
  newxyz=  newtri->xyzverts;
  newnorm= newtri->normals;
  newcolr= newtri->colors;
  newid=   newtri->cellIDs;
  newVar2= newtri->var2;
  inew= 0;
  for(i= 0; i < numold; i++) {
	i3= 3*i;
	i3new= 3*inew;
    if(nkeep[i] == 3) {
	  /* copy the input triangle to the output triangle list */
	  newxyz[i3new]= oldxyz[i3];
	  newxyz[i3new+1]= oldxyz[i3+1];
	  newxyz[i3new+2]= oldxyz[i3+2];
	  newnorm[i3new]= oldnorm[i3];
	  newnorm[i3new+1]= oldnorm[i3+1];
	  newnorm[i3new+2]= oldnorm[i3+2];
	  newcolr[i3new]= oldcolr[i3];
	  newcolr[i3new+1]= oldcolr[i3+1];
	  newcolr[i3new+2]= oldcolr[i3+2];
	  newid[inew]= oldid[i];
	  if(oldVar2) {
        newVar2[i3new]=   oldVar2[i3];
        newVar2[i3new+1]= oldVar2[i3+1];
        newVar2[i3new+2]= oldVar2[i3+2];
	  }
	  inew++;
	} else if(nkeep[i] == 2) {
	  /* two vertices are "above" the plane, so clipping results
	     in two triangles */
	  if(!keep[i3]) {
		/* the first point must be thrown out */
		p0= 0;
		p1= 1;
		p2= 2;
	  } else if(!keep[i3+1]) {
		/* the second point must be thrown out */
		p0= 1;
		p1= 2;
		p2= 0;
	  } else {
		/* the third point must be thrown out */
		p0= 2;
		p1= 0;
		p2= 1;
	  }
	  frac= -dp[i3+p0]/(dp[i3+p1]-dp[i3+p0]);
	  /* this triangle starts with two points of the existing triangle */
	  newxyz[i3new]= oldxyz[i3+p1];
	  newnorm[i3new]= oldnorm[i3+p1];
	  newxyz[i3new+1]= oldxyz[i3+p2];
	  newnorm[i3new+1]= oldnorm[i3+p2];
	  newxyz[i3new+2].x= oldxyz[i3+p0].x+frac*(oldxyz[i3+p1].x-oldxyz[i3+p0].x);
	  newxyz[i3new+2].y= oldxyz[i3+p0].y+frac*(oldxyz[i3+p1].y-oldxyz[i3+p0].y);
	  newxyz[i3new+2].z= oldxyz[i3+p0].z+frac*(oldxyz[i3+p1].z-oldxyz[i3+p0].z);
	  newnorm[i3new+2].x= oldnorm[i3+p0].x+frac*(oldnorm[i3+p1].x-oldnorm[i3+p0].x);
	  newnorm[i3new+2].y= oldnorm[i3+p0].y+frac*(oldnorm[i3+p1].y-oldnorm[i3+p0].y);
	  newnorm[i3new+2].z= oldnorm[i3+p0].z+frac*(oldnorm[i3+p1].z-oldnorm[i3+p0].z);
	  newcolr[i3new]= oldcolr[i3];
	  newcolr[i3new+1]= oldcolr[i3+1];
	  newcolr[i3new+2]= oldcolr[i3+2];
	  newid[inew]= oldid[i];
	  if(oldVar2) {
        newVar2[i3new]=   oldVar2[i3];
        newVar2[i3new+1]= oldVar2[i3+1];
        newVar2[i3new+2]= oldVar2[i3+2];
	  }
	  inew++;
	  i3new= 3*inew;

	  frac= -dp[i3+p0]/(dp[i3+p2]-dp[i3+p0]);
	  /* this triangle starts with the intersection point from the previous
	     triangle */
	  newxyz[i3new]= newxyz[i3new-1];
	  newnorm[i3new]= newnorm[i3new-1];
	  newxyz[i3new+1]= oldxyz[i3+p2];
	  newnorm[i3new+1]= oldnorm[i3+p2];
      newxyz[i3new+2].x= oldxyz[i3+p0].x+frac*(oldxyz[i3+p2].x-oldxyz[i3+p0].x);
	  newxyz[i3new+2].y= oldxyz[i3+p0].y+frac*(oldxyz[i3+p2].y-oldxyz[i3+p0].y);
	  newxyz[i3new+2].z= oldxyz[i3+p0].z+frac*(oldxyz[i3+p2].z-oldxyz[i3+p0].z);
	  newnorm[i3new+2].x= oldnorm[i3+p0].x+frac*(oldnorm[i3+p2].x-oldnorm[i3+p0].x);
	  newnorm[i3new+2].y= oldnorm[i3+p0].y+frac*(oldnorm[i3+p2].y-oldnorm[i3+p0].y);
	  newnorm[i3new+2].z= oldnorm[i3+p0].z+frac*(oldnorm[i3+p2].z-oldnorm[i3+p0].z);
	  newcolr[i3new]= oldcolr[i3];
	  newcolr[i3new+1]= oldcolr[i3+1];
	  newcolr[i3new+2]= oldcolr[i3+2];
	  newid[inew]= oldid[i];
	  if(oldVar2) {
        newVar2[i3new]=   oldVar2[i3];
        newVar2[i3new+1]= oldVar2[i3+1];
        newVar2[i3new+2]= oldVar2[i3+2];
	  }
	  inew++;
	} else if(nkeep[i] == 1) {
	  /* one vertex is "above" the plane, so clipping results
	     in one triangle */
	  /* one point of the old triangle extends above the clip plane */
	  if(keep[i3]) {
		/* the first point must be kept */
		p0= 0;
		p1= 1;
		p2= 2;
	  } else if(keep[i3+1]) {
		/* the second point must be kept */
		p0= 1;
		p1= 2;
		p2= 0;
	  } else {
		/* the third point must be kept */
		p0= 2;
		p1= 0;
		p2= 1;
	  }
	  newxyz[i3new]= oldxyz[i3+p0];
	  newnorm[i3new]= oldnorm[i3+p0];
	  /* CAUTION - next two points must be in the proper order */
	  frac= -dp[i3+p0]/(dp[i3+p1]-dp[i3+p0]);
	  newxyz[i3new+1].x= oldxyz[i3+p0].x+frac*(oldxyz[i3+p1].x-oldxyz[i3+p0].x);
	  newxyz[i3new+1].y= oldxyz[i3+p0].y+frac*(oldxyz[i3+p1].y-oldxyz[i3+p0].y);
	  newxyz[i3new+1].z= oldxyz[i3+p0].z+frac*(oldxyz[i3+p1].z-oldxyz[i3+p0].z);
	  newnorm[i3new+1].x= oldnorm[i3+p0].x+frac*(oldnorm[i3+p1].x-oldnorm[i3+p0].x);
	  newnorm[i3new+1].y= oldnorm[i3+p0].y+frac*(oldnorm[i3+p1].y-oldnorm[i3+p0].y);
	  newnorm[i3new+1].z= oldnorm[i3+p0].z+frac*(oldnorm[i3+p1].z-oldnorm[i3+p0].z);
	  frac= -dp[i3+p0]/(dp[i3+p2]-dp[i3+p0]);
	  newxyz[i3new+2].x= oldxyz[i3+p0].x+frac*(oldxyz[i3+p2].x-oldxyz[i3+p0].x);
	  newxyz[i3new+2].y= oldxyz[i3+p0].y+frac*(oldxyz[i3+p2].y-oldxyz[i3+p0].y);
	  newxyz[i3new+2].z= oldxyz[i3+p0].z+frac*(oldxyz[i3+p2].z-oldxyz[i3+p0].z);
	  newnorm[i3new+2].x= oldnorm[i3+p0].x+frac*(oldnorm[i3+p2].x-oldnorm[i3+p0].x);
	  newnorm[i3new+2].y= oldnorm[i3+p0].y+frac*(oldnorm[i3+p2].y-oldnorm[i3+p0].y);
	  newnorm[i3new+2].z= oldnorm[i3+p0].z+frac*(oldnorm[i3+p2].z-oldnorm[i3+p0].z);
	  newcolr[i3new]= oldcolr[i3];
	  newcolr[i3new+1]= oldcolr[i3+1];
	  newcolr[i3new+2]= oldcolr[i3+2];
	  newid[inew]= oldid[i];
	  if(oldVar2) {
        newVar2[i3new]=   oldVar2[i3];
        newVar2[i3new+1]= oldVar2[i3+1];
        newVar2[i3new+2]= oldVar2[i3+2];
	  }
	  inew++;
	}
  }
}

void yglMap2ColorRaw3d(long ncolr, unsigned char *red, unsigned char *green, 
					 unsigned char *blue, double vmin, double vmax, double *var, 
					 long ntri, long *cellids, float *colr)
{
  long i, ndx;
  double x;

  /*
	 This function takes as input a vector of cell numbers (cellids).
	 The number of such values is ntri.
	 red, green, and blue are a color table with ncolr elements.
	 For a cell with index ii, this function maps the variable
	 value, var[cellids[ii]], into an RGB color. 
	 The mapping takes the range from vmin to vmax and generates
	 indices from zero to ncolr.
  */
  for(i= 0; i < ntri; i++) {
    x= var[cellids[i]];
    if(x < vmin) x= vmin;
	if(x > vmax) x= vmax;
    ndx= (long) (ncolr*(x-vmin)/(vmax-vmin));
	if(ndx > ncolr-1) ndx= ncolr-1;
	colr[3*i  ]= (float) (red[ndx]  /256.0);
	colr[3*i+1]= (float) (green[ndx]/256.0);
	colr[3*i+2]= (float) (blue[ndx] /256.0);
  }
}
