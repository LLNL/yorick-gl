/*
 * $Id: slicetree.c,v 1.1.1.1 2005-09-18 22:08:03 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "Contour3D.h"
#include "slicetree.h"
#include <stdlib.h>

/* #include "ydata.h" */

static int t_octant;
static long *t_sizes;
static long *t_start;
static long *t_chunk;
static double *t_deltas;
static double *t_origin; 
static double *t_point;
static double *t_normal;
static double *t_xyz;
static TriArrayGrp *t_triangles;
#ifdef OLD_ISO_STUFF
static OctSTree *t_tree;
static long t_maxdepth;
#endif
static long *t_trsiz;
static long *t_offsets;
static OctSpan *t_ranges;
static double *t_var2;

#define DO_STATS
#ifdef DO_STATS
extern long numscan, numcross;
#endif

#ifdef OLD_ISO_STUFF
static long edges[12][2] = { {0,1}, {1,2}, {2,3}, {3,0},
                            {4,5}, {5,6}, {6,7}, {7,4},
                            {0,4}, {1,5}, {3,7}, {2,6}};
static int CASE_MASK[8] = {1,2,4,8,16,32,64,128};
#endif

int ycMakeSliceTreeCrv(double *xyz, OctSTree *tree)
{
  long maxdepth= tree->maxdepth;
  long *trsiz= tree->trsiz;
  long *offsets= tree->offsets;
  OctSpan *ranges= tree->ranges;
  long i;
  OctSpan *oldr, *newr;
  
  /* This function builds an octree from the input data values.
     The caller has pre-allocated enough space for the tree.
     The caller must free the storage.
     tree->size is the number of points in each direction in
     the whole array. 
     tree->chunk is the number of points in this chunk (one 
     more than the number of cells).
     tree->start is the point in the larger grid where this 
     chunk starts.
     tree->trsiz contains the number of "cells" in the various levels
     of the octree. The first three elements are the same as
     chunk.
  */

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();
  /* build an octree containing the data range for sub-blocks of
     the input array */
  firstSblk(tree->start, tree->size, trsiz, xyz, ranges);
  for(i= 1; i < maxdepth; i++) {
    /* NOTE: nextblk needs both the size of the tree level it is 
       filling in and the size of the parent level. 
       It also needs both the parent ranges and space to fill in 
       the ranges for this level. */
    oldr= ranges+offsets[i-1];
    newr= ranges+offsets[i];
    nextSblk(trsiz+3*(i-1), oldr, newr);
  }
  return 1;
}

void firstSblk(long *start, long *sizes, long *trsiz, double *xyz, 
			   OctSpan *rng)
{
  /* number of vertices in the full array (including guard cells) */
  long ngx= sizes[0];
  long ngy= sizes[1];
  /* starting point of this chunk in the full array */
  long stx= start[0];
  long sty= start[1];
  long stz= start[2];
  /* number of "cells" in this chunk (all are real cells) */
  long nx= trsiz[0];
  long ny= trsiz[1];
  long nz= trsiz[2];
  long rowsiz= 3*ngx;
  long plnsiz= 3*ngx*ngy;
  long i, j, k, base, nbase;
  double xmin, xmax, ymin, ymax, zmin, zmax, val;

  /* NOTE: nx, ny, and nz are numbers of cells.
     The number of vertices is one more.
     start gives the starting point of this chunk in the 
     larger array. 
     sizes is the number of vertices of the larger array.
     The range array is the same size as an array of cells 
     (for this chunk).
  */

  /* the following loop is over all cells */
  for(k= 0; k < nz; k++) {
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        /* base is a point index into the full input array */
        base= (i+stx)*3+(j+sty)*rowsiz+(k+stz)*plnsiz;
        xmin= xmax= xyz[base];
        ymin= ymax= xyz[base+1];
        zmin= zmax= xyz[base+2];
        val= xyz[base+3];
        if(val < xmin) xmin= val;
        if(val > xmax) xmax= val;
        val= xyz[base+4];
        if(val < ymin) ymin= val;
        if(val > ymax) ymax= val;
        val= xyz[base+5];
        if(val < zmin) zmin= val;
        if(val > zmax) zmax= val;
        val= xyz[base+  rowsiz];
        if(val < xmin) xmin= val;
        if(val > xmax) xmax= val;
        val= xyz[base+1+rowsiz];
        if(val < ymin) ymin= val;
        if(val > ymax) ymax= val;
        val= xyz[base+2+rowsiz];
        if(val < zmin) zmin= val;
        if(val > zmax) zmax= val;
        val= xyz[base+3+rowsiz];
        if(val < xmin) xmin= val;
        if(val > xmax) xmax= val;
        val= xyz[base+4+rowsiz];
        if(val < ymin) ymin= val;
        if(val > ymax) ymax= val;
        val= xyz[base+5+rowsiz];
        if(val < zmin) zmin= val;
        if(val > zmax) zmax= val;
        val= xyz[base         +plnsiz];
        if(val < xmin) xmin= val;
        if(val > xmax) xmax= val;
        val= xyz[base+1       +plnsiz];
        if(val < ymin) ymin= val;
        if(val > ymax) ymax= val;
        val= xyz[base+2       +plnsiz];
        if(val < zmin) zmin= val;
        if(val > zmax) zmax= val;
        val= xyz[base+3       +plnsiz];
        if(val < xmin) xmin= val;
        if(val > xmax) xmax= val;
        val= xyz[base+4       +plnsiz];
        if(val < ymin) ymin= val;
        if(val > ymax) ymax= val;
        val= xyz[base+5       +plnsiz];
        if(val < zmin) zmin= val;
        if(val > zmax) zmax= val;
        val= xyz[base  +rowsiz+plnsiz];
        if(val < xmin) xmin= val;
        if(val > xmax) xmax= val;
        val= xyz[base+1+rowsiz+plnsiz];
        if(val < ymin) ymin= val;
        if(val > ymax) ymax= val;
        val= xyz[base+2+rowsiz+plnsiz];
        if(val < zmin) zmin= val;
        if(val > zmax) zmax= val;
        val= xyz[base+3+rowsiz+plnsiz];
        if(val < xmin) xmin= val;
        if(val > xmax) xmax= val;
        val= xyz[base+4+rowsiz+plnsiz];
        if(val < ymin) ymin= val;
        if(val > ymax) ymax= val;
        val= xyz[base+5+rowsiz+plnsiz];
        if(val < zmin) zmin= val;
        if(val > zmax) zmax= val;
        /* nbase is an index into this chunk of cells */
        nbase= i+j*nx+k*nx*ny;
        rng[nbase].xmin= xmin;
        rng[nbase].xmax= xmax;
        rng[nbase].ymin= ymin;
        rng[nbase].ymax= ymax;
        rng[nbase].zmin= zmin;
        rng[nbase].zmax= zmax;
      }
    }
  }
}

void nextSblk(long *trsiz, OctSpan *oldr, OctSpan *nrng)
{
  long oldnx= trsiz[0];
  long oldny= trsiz[1];
  long oldnz= trsiz[2];
  long newnx= trsiz[3];
  long newny= trsiz[4];
  long i, j, k, base, nbase;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  OctSpan span;

  /* NOTE: nx, ny, and nz are numbers of cells.
     The number of vertices is one more.
	 The range array is the same size as an array of cells.
  */

  /* Run through all 2x2x2 groups of cells.
     NOTE: The loops must check for cases where the block
     being collapsed isn't a full 2x2x2.
  */
  for(k= 0; k < oldnz-1; k+= 2) {
    for(j= 0; j < oldny-1; j += 2) {
      for(i= 0; i < oldnx-1; i += 2) {
	/* "cell cluster" index for parent level */
       	base= i+j*oldnx+k*oldnx*oldny;
       	xmin= oldr[base].xmin;
       	xmax= oldr[base].xmax;
       	ymin= oldr[base].ymin;
       	ymax= oldr[base].ymax;
       	zmin= oldr[base].zmin;
       	zmax= oldr[base].zmax;
       	span= oldr[base+1];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+oldnx];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+1+oldnx];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+1+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+oldnx+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+1+oldnx+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	/* "cell cluster" index for child level */
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].xmin= xmin;
       	nrng[nbase].xmax= xmax;
       	nrng[nbase].ymin= ymin;
       	nrng[nbase].ymax= ymax;
       	nrng[nbase].zmin= zmin;
       	nrng[nbase].zmax= zmax;
      }
      if(oldnx & 1) {
       	i= oldnx-1;
       	base= i+j*oldnx+k*oldnx*oldny;
       	xmin= oldr[base].xmin;
       	xmax= oldr[base].xmax;
       	ymin= oldr[base].ymin;
       	ymax= oldr[base].ymax;
       	zmin= oldr[base].zmin;
       	zmax= oldr[base].zmax;
       	span= oldr[base+oldnx];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+oldnx+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].xmin= xmin;
       	nrng[nbase].xmax= xmax;
       	nrng[nbase].ymin= ymin;
       	nrng[nbase].ymax= ymax;
       	nrng[nbase].zmin= zmin;
       	nrng[nbase].zmax= zmax;
      }
    }
    if(oldny & 1) {
      j= oldny-1;
      for(i= 0; i < oldnx-1; i += 2) {
       	base= i+j*oldnx+k*oldnx*oldny;
       	xmin= oldr[base].xmin;
       	xmax= oldr[base].xmax;
       	ymin= oldr[base].ymin;
       	ymax= oldr[base].ymax;
       	zmin= oldr[base].zmin;
       	zmax= oldr[base].zmax;
       	span= oldr[base+1];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+1+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].xmin= xmin;
       	nrng[nbase].xmax= xmax;
       	nrng[nbase].ymin= ymin;
       	nrng[nbase].ymax= ymax;
       	nrng[nbase].zmin= zmin;
       	nrng[nbase].zmax= zmax;
      }
      if(oldnx & 1) {
       	i= oldnx-1;
       	base= i+j*oldnx+k*oldnx*oldny;
       	xmin= oldr[base].xmin;
       	xmax= oldr[base].xmax;
       	ymin= oldr[base].ymin;
       	ymax= oldr[base].ymax;
       	zmin= oldr[base].zmin;
       	zmax= oldr[base].zmax;
       	span= oldr[base+oldnx*oldny];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].xmin= xmin;
       	nrng[nbase].xmax= xmax;
       	nrng[nbase].ymin= ymin;
       	nrng[nbase].ymax= ymax;
       	nrng[nbase].zmin= zmin;
       	nrng[nbase].zmax= zmax;
      }
    }
  }
  if(oldnz & 1) {
    k= oldnz-1;
    for(j= 0; j < oldny-1; j += 2) {
      for(i= 0; i < oldnx-1; i += 2) {
       	base= i+j*oldnx+k*oldnx*oldny;
       	xmin= oldr[base].xmin;
       	xmax= oldr[base].xmax;
       	ymin= oldr[base].ymin;
       	ymax= oldr[base].ymax;
       	zmin= oldr[base].zmin;
       	zmax= oldr[base].zmax;
       	span= oldr[base+1];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+oldnx];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	span= oldr[base+1+oldnx];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].xmin= xmin;
       	nrng[nbase].xmax= xmax;
       	nrng[nbase].ymin= ymin;
       	nrng[nbase].ymax= ymax;
       	nrng[nbase].zmin= zmin;
       	nrng[nbase].zmax= zmax;
      }
      if(oldnx & 1) {
       	i= oldnx-1;
       	base= i+j*oldnx+k*oldnx*oldny;
       	xmin= oldr[base].xmin;
       	xmax= oldr[base].xmax;
       	ymin= oldr[base].ymin;
       	ymax= oldr[base].ymax;
       	zmin= oldr[base].zmin;
       	zmax= oldr[base].zmax;
       	span= oldr[base+oldnx];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].xmin= xmin;
       	nrng[nbase].xmax= xmax;
       	nrng[nbase].ymin= ymin;
       	nrng[nbase].ymax= ymax;
       	nrng[nbase].zmin= zmin;
       	nrng[nbase].zmax= zmax;
      }
    }
    if(oldny & 1) {
      j= oldny-1;
      for(i= 0; i < oldnx-1; i += 2) {
       	base= i+j*oldnx+k*oldnx*oldny;
       	xmin= oldr[base].xmin;
       	xmax= oldr[base].xmax;
       	ymin= oldr[base].ymin;
       	ymax= oldr[base].ymax;
       	zmin= oldr[base].zmin;
       	zmax= oldr[base].zmax;
       	span= oldr[base+1];
       	if(xmin > span.xmin) xmin= span.xmin;
       	if(xmax < span.xmax) xmax= span.xmax;
       	if(ymin > span.ymin) ymin= span.ymin;
       	if(ymax < span.ymax) ymax= span.ymax;
       	if(zmin > span.zmin) zmin= span.zmin;
       	if(zmax < span.zmax) zmax= span.zmax;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].xmin= xmin;
       	nrng[nbase].xmax= xmax;
       	nrng[nbase].ymin= ymin;
       	nrng[nbase].ymax= ymax;
       	nrng[nbase].zmin= zmin;
       	nrng[nbase].zmax= zmax;
      }
      if(oldnx & 1) {
       	i= oldnx-1;
       	base= i+j*oldnx+k*oldnx*oldny;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].xmin= oldr[base].xmin;
       	nrng[nbase].xmax= oldr[base].xmax;
       	nrng[nbase].ymin= oldr[base].ymin;
       	nrng[nbase].ymax= oldr[base].ymax;
       	nrng[nbase].zmin= oldr[base].zmin;
       	nrng[nbase].zmax= oldr[base].zmax;
      }
    }
  }
}

int ycSliceTreeCrv(double point[3], double normal[3], double *xyz, 
		   double *var, TriArrayGrp *triangles, OctSTree *tree)
{
  long depth;

  /* This function finds the triangles making up a slicing
     plane using a pre-computed octree.
	 If var2 is non-zero, it will be interpolated to the
	 vertices on the slicing plane.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_xyz= xyz;
  t_var2= var;
  t_point= point;
  t_normal= normal;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  if(t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  t_octant= 0;
  if(t_normal[0] >= 0) t_octant |= 1;
  if(t_normal[1] >= 0) t_octant |= 2;
  if(t_normal[2] >= 0) t_octant |= 4;
  t_triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_SblkCrv(0, 0, 0, depth);

  if(t_triangles->numTri) return 1;
  else return 0;
}

long do_SblkCrv(long i, long j, long k, long depth)
{
  long child_depth, nx, ny, nz, ndx, newi, newj, newk;
  long ilo, ihi, jlo, jhi, klo, khi;
  double s0, s1;
  OctSpan *the_rng;

  /* Check to see if the slicing plane passes through the
     chunk of grid between xmin, xmax, ymin, ymax, zmin, zmax.
	 This can be done by testing the two corners that are 
	 most nearly along the direction of the plane's normal.
     Immediately reject sub-blocks that are not cut.
	 Descend sub-blocks that are cut until reaching a leaf.
  */
#ifdef DO_STATS
  numscan++;
#endif
  nx= t_trsiz[3*depth];
  ny= t_trsiz[3*depth+1];
  nz= t_trsiz[3*depth+2];
  ndx= i+j*nx+k*nx*ny;
  the_rng= t_ranges+t_offsets[depth]+ndx;
  
  switch(t_octant) {
  case 0:
    /* the normal is in the neg-x, neg-y, neg-z quadrant */
    s0=  (the_rng->xmax-t_point[0])*t_normal[0]
        +(the_rng->ymax-t_point[1])*t_normal[1]
        +(the_rng->zmax-t_point[2])*t_normal[2];
    s1=  (the_rng->xmin-t_point[0])*t_normal[0]
        +(the_rng->ymin-t_point[1])*t_normal[1]
        +(the_rng->zmin-t_point[2])*t_normal[2];
    break;
  case 1:
    /* the normal is in the plus-x, neg-y, neg-z quadrant */
    s0=  (the_rng->xmin-t_point[0])*t_normal[0]
        +(the_rng->ymax-t_point[1])*t_normal[1]
        +(the_rng->zmax-t_point[2])*t_normal[2];
    s1=  (the_rng->xmax-t_point[0])*t_normal[0]
        +(the_rng->ymin-t_point[1])*t_normal[1]
        +(the_rng->zmin-t_point[2])*t_normal[2];
    break;
  case 2:
    /* the normal is in the neg-x, plus-y, neg-z quadrant */
    s0=  (the_rng->xmax-t_point[0])*t_normal[0]
        +(the_rng->ymin-t_point[1])*t_normal[1]
        +(the_rng->zmax-t_point[2])*t_normal[2];
    s1=  (the_rng->xmin-t_point[0])*t_normal[0]
        +(the_rng->ymax-t_point[1])*t_normal[1]
        +(the_rng->zmin-t_point[2])*t_normal[2];
    break;
  case 3:
    /* the normal is in the plus-x, plus-y, neg-z quadrant */
    s0=  (the_rng->xmin-t_point[0])*t_normal[0]
        +(the_rng->ymin-t_point[1])*t_normal[1]
        +(the_rng->zmax-t_point[2])*t_normal[2];
    s1=  (the_rng->xmax-t_point[0])*t_normal[0]
        +(the_rng->ymax-t_point[1])*t_normal[1]
        +(the_rng->zmin-t_point[2])*t_normal[2];
    break;
  case 4:
    /* the normal is in the neg-x, neg-y, plus-z quadrant */
    s0=  (the_rng->xmax-t_point[0])*t_normal[0]
        +(the_rng->ymax-t_point[1])*t_normal[1]
        +(the_rng->zmin-t_point[2])*t_normal[2];
    s1=  (the_rng->xmin-t_point[0])*t_normal[0]
        +(the_rng->ymin-t_point[1])*t_normal[1]
        +(the_rng->zmax-t_point[2])*t_normal[2];
    break;
  case 5:
    /* the normal is in the plus-x, neg-y, plus-z quadrant */
    s0=  (the_rng->xmin-t_point[0])*t_normal[0]
        +(the_rng->ymax-t_point[1])*t_normal[1]
        +(the_rng->zmin-t_point[2])*t_normal[2];
    s1=  (the_rng->xmax-t_point[0])*t_normal[0]
        +(the_rng->ymin-t_point[1])*t_normal[1]
        +(the_rng->zmax-t_point[2])*t_normal[2];
    break;
  case 6:
    /* the normal is in the neg-x, plus-y, plus-z quadrant */
    s0=  (the_rng->xmax-t_point[0])*t_normal[0]
        +(the_rng->ymin-t_point[1])*t_normal[1]
        +(the_rng->zmin-t_point[2])*t_normal[2];
    s1=  (the_rng->xmin-t_point[0])*t_normal[0]
        +(the_rng->ymax-t_point[1])*t_normal[1]
	+(the_rng->zmax-t_point[2])*t_normal[2];
    break;
  case 7:
    /* the normal is in the plus-x, plus-y, plus-z quadrant */
    s0=  (the_rng->xmin-t_point[0])*t_normal[0]
	+(the_rng->ymin-t_point[1])*t_normal[1]
	+(the_rng->zmin-t_point[2])*t_normal[2];
    s1=  (the_rng->xmax-t_point[0])*t_normal[0]
        +(the_rng->ymax-t_point[1])*t_normal[1]
	+(the_rng->zmax-t_point[2])*t_normal[2];
    break;
  }
  /* If the corners are on the same side of the plane,
     the block is not cut. No need to sub-divide. 
     Treat the case with s0 == 0 a not cut and the 
     case with s1 == 0 as cut. */
  if(s0 >= 0.0 || s1 < 0.0) return 0;
  if(depth == 0) {
    /* Have reached a leaf of the tree that is cut by the
       slice plane, so generate triangles.  */
    grab_StrisCrv(i, j, k);
    return 1;
  }
  
  /* The slicing plane cuts this block. Recursively call this
     function for all sub-blocks. */
  child_depth= depth-1;  /* yorick style depth of the child blocks */
  nx= t_trsiz[3*child_depth];
  ny= t_trsiz[3*child_depth+1];
  nz= t_trsiz[3*child_depth+2];
  klo= 2*k;  khi= klo+1; if(khi > nz-1) khi= nz-1;
  jlo= 2*j;  jhi= jlo+1; if(jhi > ny-1) jhi= ny-1;
  ilo= 2*i;  ihi= ilo+1; if(ihi > nx-1) ihi= nx-1;
  for(newk= klo; newk <= khi; newk++) {
    for(newj= jlo; newj <= jhi; newj++) {
      for(newi= ilo; newi <= ihi; newi++) {
        do_SblkCrv(newi, newj, newk, child_depth);
      }
    }
  }
  return 1;
}

long grab_StrisCrv(long i, long j, long k)
{
  long nx= t_sizes[0];
  long ny= t_sizes[1];
  double s[8], vars[8], *var2;
  int rowsiz, plnsiz, sliceSize;
  int idxvg, idxcg, idv, ii, index, mask;
  long numTri, *cellIDs;
  yPoint3D *xyzverts, pts[8];

  cellIDs = t_triangles->cellIDs;
  xyzverts = t_triangles->xyzverts;
  var2= t_triangles->var2;
  numTri= t_triangles->numTri;
#ifdef DO_STATS
  numcross++;
#endif

  /*
     This cell is (probably) cut by the slicing plane.
     Generate the proper triangles using MarchingCubes style
     case tables.
  */
  /* NOTE: on input i,j,k are indices into the chunk 
     that is currently being processed. 
     Convert into indices into the full grid. */
  i += t_start[0];
  j += t_start[1];
  k += t_start[2];
  sliceSize = nx*ny;
  rowsiz= 3*nx;
  plnsiz= 3*nx*ny;
  /* idxvg is a vertex index into the global array */
  idxvg = 3*(i + j*nx + k*sliceSize);
  /* compute the distance to the plane at the 8 corners */
  s[0] =  (t_xyz[idxvg  ]-t_point[0])*t_normal[0]
	 +(t_xyz[idxvg+1]-t_point[1])*t_normal[1]
	 +(t_xyz[idxvg+2]-t_point[2])*t_normal[2];
  s[1] =  (t_xyz[idxvg+3]-t_point[0])*t_normal[0]
         +(t_xyz[idxvg+4]-t_point[1])*t_normal[1]
         +(t_xyz[idxvg+5]-t_point[2])*t_normal[2];
  s[2] =  (t_xyz[idxvg+3+rowsiz]-t_point[0])*t_normal[0]
	 +(t_xyz[idxvg+4+rowsiz]-t_point[1])*t_normal[1]
	 +(t_xyz[idxvg+5+rowsiz]-t_point[2])*t_normal[2];
  s[3] =  (t_xyz[idxvg  +rowsiz]-t_point[0])*t_normal[0]
	 +(t_xyz[idxvg+1+rowsiz]-t_point[1])*t_normal[1]
	 +(t_xyz[idxvg+2+rowsiz]-t_point[2])*t_normal[2];
  s[4] =  (t_xyz[idxvg       +plnsiz]-t_point[0])*t_normal[0]
	 +(t_xyz[idxvg+1     +plnsiz]-t_point[1])*t_normal[1]
	 +(t_xyz[idxvg+2     +plnsiz]-t_point[2])*t_normal[2];
  s[5] =  (t_xyz[idxvg+3     +plnsiz]-t_point[0])*t_normal[0]
	 +(t_xyz[idxvg+4     +plnsiz]-t_point[1])*t_normal[1]
	 +(t_xyz[idxvg+5     +plnsiz]-t_point[2])*t_normal[2];
  s[6] =  (t_xyz[idxvg+3+rowsiz+plnsiz]-t_point[0])*t_normal[0]
	 +(t_xyz[idxvg+4+rowsiz+plnsiz]-t_point[1])*t_normal[1]
         +(t_xyz[idxvg+5+rowsiz+plnsiz]-t_point[2])*t_normal[2];
  s[7] =  (t_xyz[idxvg+  rowsiz+plnsiz]-t_point[0])*t_normal[0]
         +(t_xyz[idxvg+1+rowsiz+plnsiz]-t_point[1])*t_normal[1]
         +(t_xyz[idxvg+2+rowsiz+plnsiz]-t_point[2])*t_normal[2];

  /* Build the case table */
  for ( ii=0, index = 0, mask= 1; ii < 8; ii++, mask += mask) {
    if ( s[ii] >= 0.0 ) 
      index |= mask;
  }
#define DEBUG_ISO
#ifdef DEBUG_ISO
  if(index == 0 || index == 255) {
    /* because the test in the outer loop was against
       the bounding box, this cell may not be cut. */
    return 0;
  }
#endif

  /* set coordinates for the corners of the cell */
  pts[0].x= t_xyz[idxvg  ];
  pts[0].y= t_xyz[idxvg+1];
  pts[0].z= t_xyz[idxvg+2];
  pts[1].x= t_xyz[idxvg+3];
  pts[1].y= t_xyz[idxvg+4];
  pts[1].z= t_xyz[idxvg+5];
  pts[2].x= t_xyz[idxvg+3+rowsiz];
  pts[2].y= t_xyz[idxvg+4+rowsiz];
  pts[2].z= t_xyz[idxvg+5+rowsiz];
  pts[3].x= t_xyz[idxvg  +rowsiz];
  pts[3].y= t_xyz[idxvg+1+rowsiz];
  pts[3].z= t_xyz[idxvg+2+rowsiz];
  pts[4].x= t_xyz[idxvg       +plnsiz];
  pts[4].y= t_xyz[idxvg+1     +plnsiz];
  pts[4].z= t_xyz[idxvg+2     +plnsiz];
  pts[5].x= t_xyz[idxvg+3     +plnsiz];
  pts[5].y= t_xyz[idxvg+4     +plnsiz];
  pts[5].z= t_xyz[idxvg+5     +plnsiz];
  pts[6].x= t_xyz[idxvg+3+rowsiz+plnsiz];
  pts[6].y= t_xyz[idxvg+4+rowsiz+plnsiz];
  pts[6].z= t_xyz[idxvg+5+rowsiz+plnsiz];
  pts[7].x= t_xyz[idxvg+  rowsiz+plnsiz];
  pts[7].y= t_xyz[idxvg+1+rowsiz+plnsiz];
  pts[7].z= t_xyz[idxvg+2+rowsiz+plnsiz];
  /* set variable at corners of cell if provided */
  idv = i + j*nx + k*sliceSize;
  if(t_var2) {
    vars[0]= t_var2[idv];
    vars[1]= t_var2[idv+1];
    vars[2]= t_var2[idv+1+nx];
    vars[3]= t_var2[idv  +nx];
    vars[4]= t_var2[idv     +sliceSize];
    vars[5]= t_var2[idv+1   +sliceSize];
    vars[6]= t_var2[idv+1+nx+sliceSize];
    vars[7]= t_var2[idv  +nx+sliceSize];
  }

  idxcg= i+j*(nx-1)+k*(nx-1)*(ny-1);  
  extract_slicetris_tet(index, idxcg, t_var2, &numTri, s, pts, 
               vars, cellIDs, xyzverts, var2);
  
  /* WARNING: caller must be sure triangle array doesn't overflow. */
  t_triangles->numTri= numTri;
  return 0;
}

int ycSliceTree(long maxdepth, long sizes[3], long chunk[3], long start[3],
                double deltas[3], double origin[3], double point[3], 
                double normal[3], double *var, TriArrayGrp *triangles)
{
  long depth;
  long ilo, ihi, jlo, jhi, klo, khi;

  /* This function draws a slicing plane through a regular
     orthogonal grid. It uses a binary search for grid 
     chunks that are cut by the plane. No Octree is needed
     because the grid is completely regular.

     The first step is to find the octant in which the normal
     lies. This determines the two corners that extends
     furthest along the normal, and only those two must be
     tested to see if the plane cuts a chunk.
  */
  t_sizes= sizes;
  t_chunk= chunk;
  t_start= start;
  t_deltas= deltas;
  t_origin= origin; 
  t_point= point;
  t_normal= normal;
  t_var2= var;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_maxdepth= maxdepth;
#endif
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  if(t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  t_octant= 0;
  if(t_normal[0] >= 0) t_octant |= 1;
  if(t_normal[1] >= 0) t_octant |= 2;
  if(t_normal[2] >= 0) t_octant |= 4;
  t_triangles->numTri= 0;
  depth= maxdepth-1;  /* depth as a C-style index */
  /* The initial block includes the whole chunk */
  ilo= t_start[0];
  ihi= t_start[0]+t_chunk[0]-1;
  jlo= t_start[1];
  jhi= t_start[1]+t_chunk[1]-1;
  klo= t_start[2];
  khi= t_start[2]+t_chunk[2]-1;
  do_Sblk(ilo, ihi, jlo, jhi, klo, khi, depth);

  if(t_triangles->numTri) return 1;
  else return 0;
}

long do_Sblk(long ilo, long ihi, long jlo, long jhi, long klo,
             long khi, long depth)
{
  long child_depth, imid, jmid, kmid;
  double s0, s1, s[8];
  double delx= t_deltas[0];
  double dely= t_deltas[1];
  double delz= t_deltas[2];
  double offx= t_origin[0]-t_point[0];
  double offy= t_origin[1]-t_point[1];
  double offz= t_origin[2]-t_point[2];

  /* Check to see if the slicing plane passes through the
     chunk of grid between ilo and ihi, jlo and jhi, 
     and klo and khi.
     This can be done by testing at the two corners that are 
     most nearly along the direction of the plane's normal.
     Immediately reject sub-blocks that are not cut.
     Descend sub-blocks that are cut until reaching a leaf.
  */
  if(ilo >= ihi || jlo >= jhi || klo >= khi) return 0;
#ifdef DO_STATS
  numscan++;
#endif
  /*
     NOTE: with a Cartesian grid, intersection of the slice
     plane and the zone can be checked for by evaluating the
     distance to the plane at only 2 vertices.
     s0 is the distance to the "entrance" vertex and
     s1 is the distance to the "exit" vertex. The plane
     cuts the "cell cluster" if s0 < 0 and s1 > 0.
     A simpler test is that s0 and s1 have different signs.
  */
  switch(t_octant) {
  case 0:
    /* the normal is in the neg-x, neg-y, neg-z quadrant */
    s0=  (ihi*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
        +(khi*delz+offz)*t_normal[2];
    s1=  (ilo*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
        +(klo*delz+offz)*t_normal[2];
    break;
  case 1:
    /* the normal is in the plus-x, neg-y, neg-z quadrant */
    s0=  (ilo*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
        +(khi*delz+offz)*t_normal[2];
    s1=  (ihi*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
        +(klo*delz+offz)*t_normal[2];
    break;
  case 2:
    /* the normal is in the neg-x, plus-y, neg-z quadrant */
    s0=  (ihi*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
        +(khi*delz+offz)*t_normal[2];
    s1=  (ilo*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
        +(klo*delz+offz)*t_normal[2];
    break;
  case 3:
    /* the normal is in the plus-x, plus-y, neg-z quadrant */
    s0=  (ilo*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
        +(khi*delz+offz)*t_normal[2];
    s1=  (ihi*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
        +(klo*delz+offz)*t_normal[2];
    break;
  case 4:
    /* the normal is in the neg-x, neg-y, plus-z quadrant */
    s0=  (ihi*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
        +(klo*delz+offz)*t_normal[2];
    s1=  (ilo*delx+offx)*t_normal[0]
        +(jlo*dely+offy)*t_normal[1]+(khi*delz+offz)*t_normal[2];
    break;
  case 5:
    /* the normal is in the plus-x, neg-y, plus-z quadrant */
    s0=  (ilo*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
        +(klo*delz+offz)*t_normal[2];
    s1=  (ihi*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
        +(khi*delz+offz)*t_normal[2];
    break;
  case 6:
    /* the normal is in the neg-x, plus-y, plus-z quadrant */
    s0=  (ihi*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
        +(klo*delz+offz)*t_normal[2];
    s1=  (ilo*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
        +(khi*delz+offz)*t_normal[2];
    break;
  case 7:
    /* the normal is in the plus-x, plus-y, plus-z quadrant */
    s0=  (ilo*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
        +(klo*delz+offz)*t_normal[2];
    s1=  (ihi*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
        +(khi*delz+offz)*t_normal[2];
    break;
  }
#define DBG_S0_S1
#ifdef DBG_S0_S1
  if(s0 > s1) {
    /* impossible error */
    return 0;
  }
#endif
  /* If the corners are on the same side of the plane,
     the block is not cut. No need to sub-divide. 
     Treat the case with s0 == 0 as not cut and the 
     case with s1 == 0 as cut. */
  if(s0 >= 0.0 || s1 < 0.0) return 0;
  if(ihi-ilo == 1 && jhi-jlo == 1 && khi-klo == 1) {
    /* Have reached a leaf of the tree that is cut by the
       slice plane, so generate triangles. Be sure to
       do this in exactly the same way as the computation
       of s0 and s1 above. */
    s[0]= (ilo*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
         +(klo*delz+offz)*t_normal[2];
    s[1]= (ihi*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
         +(klo*delz+offz)*t_normal[2];
    s[2]= (ihi*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
         +(klo*delz+offz)*t_normal[2];
    s[3]= (ilo*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
         +(klo*delz+offz)*t_normal[2];
    s[4]= (ilo*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
         +(khi*delz+offz)*t_normal[2];
    s[5]= (ihi*delx+offx)*t_normal[0]+(jlo*dely+offy)*t_normal[1]
         +(khi*delz+offz)*t_normal[2];
    s[6]= (ihi*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
         +(khi*delz+offz)*t_normal[2];
    s[7]= (ilo*delx+offx)*t_normal[0]+(jhi*dely+offy)*t_normal[1]
         +(khi*delz+offz)*t_normal[2];
    grab_Stris(ilo, jlo, klo, s);
    return 1;
  }
  
  /* The slicing plane cuts this block. Recursively call this
     function for all sub-blocks. */
  child_depth= depth-1;  /* yorick style depth of the child blocks */
  imid= (ilo+ihi+1)/2;
  jmid= (jlo+jhi+1)/2;
  kmid= (klo+khi+1)/2;
  /* Call do_Sblk with a range of cells to handle. Empty ranges
     are rejected at the top of do_Sblk. */
  do_Sblk(ilo,  imid, jlo,  jmid, klo,  kmid, child_depth);
  do_Sblk(imid, ihi,  jlo,  jmid, klo,  kmid, child_depth);
  do_Sblk(ilo,  imid, jmid, jhi,  klo,  kmid, child_depth);
  do_Sblk(imid, ihi,  jmid, jhi,  klo,  kmid, child_depth);
  do_Sblk(ilo,  imid, jlo,  jmid, kmid, khi, child_depth);
  do_Sblk(imid, ihi,  jlo,  jmid, kmid, khi, child_depth);
  do_Sblk(ilo,  imid, jmid, jhi,  kmid, khi, child_depth);
  do_Sblk(imid, ihi,  jmid, jhi,  kmid, khi, child_depth);
  return 1;
}

long grab_Stris(long i, long j, long k, double s[8])
{
  /* nx, ny, nz are the number of vertices in the global array */
  long nx= t_sizes[0];
  long ny= t_sizes[1];
  double dx= t_deltas[0];
  double dy= t_deltas[1];
  double dz= t_deltas[2];
  double x0= t_origin[0];
  double y0= t_origin[1];
  double z0= t_origin[2];
  double vars[8];
  int idxcg, idxvg, ii, index, mask;
  double *var2;
  long numTri, sliceSize, *cellIDs;
  yPoint3D *xyzverts, pts[8];

  cellIDs = t_triangles->cellIDs;
  xyzverts = t_triangles->xyzverts;
  var2= t_triangles->var2;
  numTri= t_triangles->numTri;
#ifdef DO_STATS
  numcross++;
#endif

  /*
     This cell is cut by the slicing plane.
     Generate the proper triangles using MarchingCubes style case tables.
  */

  /* Build the case table */
  for ( ii=0, index = 0, mask= 1; ii < 8; ii++, mask += mask) {
    if ( s[ii] >= 0.0 ) 
      index |= mask;
  }
#define DEBUG_ISO
#ifdef DEBUG_ISO
  if(index == 0 || index == 255) {
    YError("contouring a zone that should have been rejected");
  }
#endif

  /* set coordinates for the corners of the cell */
  pts[0].z = z0 + k*dz;
  pts[0].y = y0 + j*dy;
  pts[0].x = x0 + i*dx;

  pts[1].x = pts[0].x + dx;  
  pts[1].y = pts[0].y;
  pts[1].z = pts[0].z;

  pts[2].x = pts[0].x + dx;  
  pts[2].y = pts[0].y + dy;
  pts[2].z = pts[0].z;

  pts[3].x = pts[0].x;
  pts[3].y = pts[0].y + dy;
  pts[3].z = pts[0].z;

  pts[4].x = pts[0].x;
  pts[4].y = pts[0].y;
  pts[4].z = pts[0].z + dz;

  pts[5].x = pts[0].x + dx;  
  pts[5].y = pts[0].y;
  pts[5].z = pts[0].z + dz;

  pts[6].x = pts[0].x + dx;  
  pts[6].y = pts[0].y + dy;
  pts[6].z = pts[0].z + dz;

  pts[7].x = pts[0].x;
  pts[7].y = pts[0].y + dy;
  pts[7].z = pts[0].z + dz;
  /* set variable at corners of cell if provided */
  sliceSize = nx*ny;
  idxvg = i + j*nx + k*sliceSize;
  if(t_var2) {
    vars[0]= t_var2[idxvg];
    vars[1]= t_var2[idxvg+1];
    vars[2]= t_var2[idxvg+1+nx];
    vars[3]= t_var2[idxvg  +nx];
    vars[4]= t_var2[idxvg     +sliceSize];
    vars[5]= t_var2[idxvg+1   +sliceSize];
    vars[6]= t_var2[idxvg+1+nx+sliceSize];
    vars[7]= t_var2[idxvg  +nx+sliceSize];
  }

  idxcg= i+j*(nx-1)+k*(nx-1)*(ny-1);  
  extract_slicetris_tet(index, idxcg, t_var2, &numTri, s, pts, 
               vars, cellIDs, xyzverts, var2);
  
  /* WARNING: caller must be sure triangle array doesn't overflow. */
  t_triangles->numTri= numTri;
  return 0;
}
