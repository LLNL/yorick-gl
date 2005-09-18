/*
 * $Id: isotree.c,v 1.1 2005-09-18 22:08:02 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "Contour3D.h"
#include "isotree.h"
#include <stdlib.h>

static long *t_sizes;
static long *t_start;
static long *t_chunk;
static double *t_deltas;
static double *t_origin; 
static double t_level;
static double *t_var;
static double *t_vcen;
static double *t_var2;
static yPoint3D *t_xyz;
static TriArrayGrp *t_triangles;
#ifdef OLD_ISO_STUFF
static OctTree *t_tree;
static long t_maxdepth;
#endif
static long *t_trsiz;
static long *t_offsets;
static OctRange *t_ranges;
static long make_strip= 0;
static long *t_ptndx;
static long v_edg_offset[12];
 
#define DO_STATS
#ifdef DO_STATS
long numscan, numcross;
#endif

#ifdef OLD_ISO_STUFF
static long edges[12][2] = { {0,1}, {1,2}, {2,3}, {3,0},
                            {4,5}, {5,6}, {6,7}, {7,4},
                            {0,4}, {1,5}, {3,7}, {2,6}};
static int CASE_MASK[8] = {1,2,4,8,16,32,64,128};
#endif

static void extract_tris(int index, long idxvg, long *numTri,
                  double lev, double s[8], yPoint3D pts[8],
                  yPoint3D gradients[8], double vars[8], long *cellIDs,
                  yPoint3D *xyzverts, yPoint3D *normals, double *var2);

int ycMakeContourTree(double *var, OctTree *tree)
{
  long maxdepth= tree->maxdepth;
  long *chk= tree->chunk;
  long *trsiz= tree->trsiz;
  long *offsets= tree->offsets;
  OctRange *ranges= tree->ranges;
  long nx= chk[0];
  long ny= chk[1];
  long nz= chk[2];
  long i;
  OctRange *oldr, *newr;
  
  /* This function builds an octree from the input data values.
     The caller has pre-allocated enough space for the tree.
     The caller must free the storage.
     tree->size is the number of points in each direction in
     the global (whole) array. The global array has a 
     layer of guard cells around it. The guard cells
     don't have a large effect in this routine because the
     algorithm must already deal with the difference
     between chunk indices and global indices. 
     NOTE: the chunk specifies points surrounding cells to
     be contoured. The chunk will never include the first or
     last points in the global array.
     tree->chunk is the number of points in this chunk (one 
     more than the number of cells).
     tree->start is the point in the larger grid where this 
     chunk starts.
     tree->trsiz contains the number of "cells" in the various levels
     of the octree. The first three elements are the same as
     chunk.
     var contains the point-centered data values.
  */
  if (!var || nx < 4 || ny < 4 || nz < 4) {
    return 0;
  }

  /* build an octree containing the data range for sub-blocks of
     the input array */
  firstblk(var, tree->start, tree->size, trsiz, ranges);
  for(i= 1; i < maxdepth; i++) {
    /* NOTE: nextblk needs both the size of the tree level it is 
       filling in and the size of the parent level. 
       It also needs both the parent ranges and space to fill in 
       the ranges for this level. */
    oldr= ranges+offsets[i-1];
    newr= ranges+offsets[i];
    nextblk(trsiz+3*(i-1), oldr, newr);
  }
  return 1;
}

void firstblk(double *data, long *start, long *sizes, long *trsiz, OctRange *rng)
{
  /* number of "cells" in this chunk */
  long nx= trsiz[0];
  long ny= trsiz[1];
  long nz= trsiz[2];
  /* number of vertices in the full array */
  long ngx= sizes[0];
  long ngy= sizes[1];
  /* starting point of this chunk in the full array */
  long stx= start[0];
  long sty= start[1];
  long stz= start[2];
  long i, j, k, base, nbase;
  double lo, hi, val;

  /* NOTE: nx, ny, and nz are numbers of cells.
     The number of vertices is one more.
     start gives the starting point of this chunk in the 
     larger array. 
     sizes is the number of vertices of the larger array.
     The range array is the same size as an array of vertices 
     (for this chunk).
     The uppermost element of the range is never used.
  */

  for(k= 0; k < nz; k++) {
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        base= i+stx+(j+sty)*ngx+(k+stz)*ngx*ngy;
	lo= data[base];
	hi= lo;
	val= data[base+1];
	if(val < lo) lo= val;
	if(val > hi) hi= val;
	val= data[base+ngx];
        if(val < lo) lo= val;
        if(val > hi) hi= val;
       	val= data[base+1+ngx];
       	if(val < lo) lo= val;
       	if(val > hi) hi= val;
       	val= data[base+ngx*ngy];
       	if(val < lo) lo= val;
       	if(val > hi) hi= val;
       	val= data[base+1+ngx*ngy];
       	if(val < lo) lo= val;
       	if(val > hi) hi= val;
       	val= data[base+ngx+ngx*ngy];
       	if(val < lo) lo= val;
       	if(val > hi) hi= val;
       	val= data[base+1+ngx+ngx*ngy];
       	if(val < lo) lo= val;
       	if(val > hi) hi= val;
       	nbase= i+j*nx+k*nx*ny;
       	rng[nbase].lo= lo;
       	rng[nbase].hi= hi;
      }
    }
  }
}

void nextblk(long *trsiz, OctRange *oldr, OctRange *nrng)
{
  long oldnx= trsiz[0];
  long oldny= trsiz[1];
  long oldnz= trsiz[2];
  long newnx= trsiz[3];
  long newny= trsiz[4];
  long i, j, k, base, nbase;
  double lo, hi, loval, hival;

  /* NOTE: nx, ny, and nz are numbers of cells.
     The number of vertices is one more.
     The range array is the same size as an array of cells.
     The uppermost element of the range is never used.
  */

  /* Run through all 2x2x2 groups of cells.
     NOTE: The loops must check for cases where the block
     being collapsed isn't a full 2x2x2.
  */
  for(k= 0; k < oldnz-1; k+= 2) {
    for(j= 0; j < oldny-1; j += 2) {
      for(i= 0; i < oldnx-1; i += 2) {
       	base= i+j*oldnx+k*oldnx*oldny;
       	lo= oldr[base].lo;
       	hi= oldr[base].hi;
       	loval= oldr[base+1].lo;
       	hival= oldr[base+1].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+oldnx].lo;
       	hival= oldr[base+oldnx].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+1+oldnx].lo;
       	hival= oldr[base+1+oldnx].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+oldnx*oldny].lo;
       	hival= oldr[base+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+1+oldnx*oldny].lo;
       	hival= oldr[base+1+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+oldnx+oldnx*oldny].lo;
       	hival= oldr[base+oldnx+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+1+oldnx+oldnx*oldny].lo;
       	hival= oldr[base+1+oldnx+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].lo= lo;
       	nrng[nbase].hi= hi;
      }
      if(oldnx & 1) {
       	i= oldnx-1;
       	base= i+j*oldnx+k*oldnx*oldny;
       	lo= oldr[base].lo;
       	hi= oldr[base].hi;
       	loval= oldr[base+oldnx].lo;
       	hival= oldr[base+oldnx].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+oldnx*oldny].lo;
       	hival= oldr[base+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+oldnx+oldnx*oldny].lo;
       	hival= oldr[base+oldnx+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].lo= lo;
       	nrng[nbase].hi= hi;
      }
    }
    if(oldny & 1) {
      j= oldny-1;
      for(i= 0; i < oldnx-1; i += 2) {
       	base= i+j*oldnx+k*oldnx*oldny;
       	lo= oldr[base].lo;
       	hi= oldr[base].hi;
       	loval= oldr[base+1].lo;
       	hival= oldr[base+1].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+oldnx*oldny].lo;
       	hival= oldr[base+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+1+oldnx*oldny].lo;
       	hival= oldr[base+1+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].lo= lo;
       	nrng[nbase].hi= hi;
      }
      if(oldnx & 1) {
       	i= oldnx-1;
       	base= i+j*oldnx+k*oldnx*oldny;
       	lo= oldr[base].lo;
       	hi= oldr[base].hi;
       	loval= oldr[base+oldnx*oldny].lo;
       	hival= oldr[base+oldnx*oldny].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].lo= lo;
       	nrng[nbase].hi= hi;
      }
    }
  }
  if(oldnz & 1) {
    k= oldnz-1;
    for(j= 0; j < oldny-1; j += 2) {
      for(i= 0; i < oldnx-1; i += 2) {
       	base= i+j*oldnx+k*oldnx*oldny;
       	lo= oldr[base].lo;
       	hi= oldr[base].hi;
       	loval= oldr[base+1].lo;
       	hival= oldr[base+1].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+oldnx].lo;
       	hival= oldr[base+oldnx].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	loval= oldr[base+1+oldnx].lo;
       	hival= oldr[base+1+oldnx].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].lo= lo;
       	nrng[nbase].hi= hi;
      }
      if(oldnx & 1) {
       	i= oldnx-1;
       	base= i+j*oldnx+k*oldnx*oldny;
       	lo= oldr[base].lo;
       	hi= oldr[base].hi;
       	loval= oldr[base+oldnx].lo;
       	hival= oldr[base+oldnx].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].lo= lo;
       	nrng[nbase].hi= hi;
      }
    }
    if(oldny & 1) {
      j= oldny-1;
      for(i= 0; i < oldnx-1; i += 2) {
       	base= i+j*oldnx+k*oldnx*oldny;
       	lo= oldr[base].lo;
       	hi= oldr[base].hi;
       	loval= oldr[base+1].lo;
       	hival= oldr[base+1].hi;
       	if(loval < lo) lo= loval;
       	if(hival > hi) hi= hival;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].lo= lo;
       	nrng[nbase].hi= hi;
      }
      if(oldnx & 1) {
       	i= oldnx-1;
       	base= i+j*oldnx+k*oldnx*oldny;
       	lo= oldr[base].lo;
       	hi= oldr[base].hi;
       	nbase= i/2+(j/2)*newnx+(k/2)*newnx*newny;
       	nrng[nbase].lo= lo;
       	nrng[nbase].hi= hi;
      }
    }
  }
}

int ycContourTree(double deltas[3], double origin[3], double level, 
		  double *var, TriArrayGrp *triangles, OctTree *tree)
{
  long depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= deltas;
  t_origin= origin; 
  t_level= level;
  t_var= var;
  t_vcen= 0;
  t_var2= 0;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= 0;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTree2(double deltas[3], double origin[3], double level, 
		   double *var, double *var2, TriArrayGrp *triangles, 
		   OctTree *tree)
{
  long depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= deltas;
  t_origin= origin; 
  t_level= level;
  t_var= var;
  t_vcen= 0;
  t_var2= var2;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= 0;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeCrv(double level, yPoint3D *xyz, double *var, 
		     TriArrayGrp *triangles, OctTree *tree)
{
  long depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= 0;
  t_origin= 0; 
  t_level= level;
  t_var= var;
  t_vcen= 0;
  t_var2= 0;
  t_xyz= xyz;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= 0;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeCrv2(double level, yPoint3D *xyz, double *var, 
		      double *var2, TriArrayGrp *triangles, 
		      OctTree *tree)
{
  long depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= 0;
  t_origin= 0; 
  t_level= level;
  t_var= var;
  t_vcen= 0;
  t_var2= var2;
  t_xyz= xyz;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= 0;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeZcen(double deltas[3], double origin[3], double level, 
		      double *var, double *vcen, TriArrayGrp *triangles, 
		      OctTree *tree)
{
  long depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= deltas;
  t_origin= origin; 
  t_level= level;
  t_var= var;
  t_vcen= vcen;
  t_var2= 0;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= 0;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeZcen2(double deltas[3], double origin[3], double level, 
                       double *var, double *vcen, double *var2,
                       TriArrayGrp *triangles, OctTree *tree)
{
  long depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= deltas;
  t_origin= origin; 
  t_level= level;
  t_var= var;
  t_vcen= vcen;
  t_var2= var2;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= 0;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeCrvZcen(double level, yPoint3D *xyz, 
                         double *var, double *vcen, TriArrayGrp *triangles, 
                         OctTree *tree)
{
  long depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= 0;
  t_origin= 0; 
  t_level= level;
  t_var= var;
  t_vcen= vcen;
  t_var2= 0;
  t_xyz= xyz;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= 0;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeCrvZcen2(double level, yPoint3D *xyz, double *var, 
                          double *vcen, double *var2, TriArrayGrp *triangles, 
                          OctTree *tree)
{
  long depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= 0;
  t_origin= 0; 
  t_level= level;
  t_var= var;
  t_vcen= vcen;
  t_var2= var2;
  t_xyz= xyz;
  t_triangles= triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= 0;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeVarr(double deltas[3], double origin[3], double level, 
                      double *var, TriVertexGrp *triangles, 
                      OctTree *tree, long *edgndx)
{
  long i, nedg, depth, nx, ny, xy_siz;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();
  
  /* This function draws an iso-surface using a pre-computed octree. */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= deltas;
  t_origin= origin; 
  t_level= level;
  t_var= var;
  t_vcen= 0;
  t_var2= 0;
  t_triangles= (TriArrayGrp *) triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= edgndx;
  nx= t_chunk[0];
  ny= t_chunk[1];
  xy_siz= nx*ny;
  v_edg_offset[0]= 0;
  v_edg_offset[1]= 4;
  v_edg_offset[2]= 3*nx;
  v_edg_offset[3]= 1;
  v_edg_offset[4]= 3*xy_siz;
  v_edg_offset[5]= 3*xy_siz+4;
  v_edg_offset[6]= 3*xy_siz+3*nx;
  v_edg_offset[7]= 3*xy_siz+1;
  v_edg_offset[8]= 2;
  v_edg_offset[9]= 5;
  v_edg_offset[10]= 3*nx+5;
  v_edg_offset[11]= 3*nx+2;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  /* MUST clear both the number of triangles and the number of edges cut here */
  triangles->numTri= 0;
  triangles->numEdg= 0;
  /* adjust so that the array indicates no vertices have yet been found */
  nedg= 3*t_chunk[0]*t_chunk[1]*t_chunk[2];
  for(i= 0; i <= nedg; i++) edgndx[i]= -1;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeVarr2(double deltas[3], double origin[3], double level, 
                       double *var, double *var2, TriVertexGrp *triangles, 
                       OctTree *tree, long *edgndx)
{
  long i, nedg, depth, nx, ny, xy_siz;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= deltas;
  t_origin= origin; 
  t_level= level;
  t_var= var;
  t_vcen= 0;
  t_var2= var2;
  t_triangles= (TriArrayGrp *) triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= edgndx;
  nx= t_chunk[0];
  ny= t_chunk[1];
  xy_siz= nx*ny;
  v_edg_offset[0]= 0;
  v_edg_offset[1]= 4;
  v_edg_offset[2]= 3*nx;
  v_edg_offset[3]= 1;
  v_edg_offset[4]= 3*xy_siz;
  v_edg_offset[5]= 3*xy_siz+4;
  v_edg_offset[6]= 3*xy_siz+3*nx;
  v_edg_offset[7]= 3*xy_siz+1;
  v_edg_offset[8]= 2;
  v_edg_offset[9]= 5;
  v_edg_offset[10]= 3*nx+5;
  v_edg_offset[11]= 3*nx+2;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  /* MUST clear both the number of triangles and the number of edges cut here */
  triangles->numTri= 0;
  triangles->numEdg= 0;
  /* adjust so that the array indicates no vertices have yet been found */
  nedg= 3*t_chunk[0]*t_chunk[1]*t_chunk[2];
  for(i= 0; i <= nedg; i++) edgndx[i]= -1;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

int ycContourTreeVarrStr(double deltas[3], double origin[3], double level, 
                         double *var, TriVertexGrp *triangles, 
                         OctTree *tree, long *edgndx)
{
  long i, nedg, depth;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /* This function draws an iso-surface using a pre-computed
     octree.
  */
  t_sizes= tree->size;
  t_chunk= tree->chunk;
  t_start= tree->start;
  t_deltas= deltas;
  t_origin= origin; 
  t_level= level;
  t_var= var;
  t_vcen= 0;
  t_var2= 0;
  t_triangles= (TriArrayGrp *) triangles;
#ifdef OLD_ISO_STUFF
  t_tree= tree;
  t_maxdepth= tree->maxdepth;
#endif
  t_trsiz= tree->trsiz;
  t_offsets= tree->offsets;
  t_ranges= tree->ranges;
  t_ptndx= edgndx;
#ifdef DO_STATS
  numscan= 0;
  numcross= 0;
#endif
    
  /* NOTE: a chunk does not include any guard cells */
  if (!var || t_chunk[0] < 2 || t_chunk[1] < 2 || t_chunk[2] < 2) {
    return 0;
  }
  triangles->numTri= 0;
  triangles->numEdg= 0;
  /* adjust so that the array indicates no vertices have yet been found */
  nedg= 3*t_chunk[0]*t_chunk[1]*t_chunk[2];
  for(i= 0; i <= nedg; i++) edgndx[i]= -1;
  depth= tree->maxdepth-1;  /* depth as a C-style index */
  do_blk(0, 0, 0, depth);

  if(triangles->numTri) return 1;
  else return 0;
}

long do_blk(long i, long j, long k, long depth)
{
  long child_depth, nx, ny, nz, ndx, newi, newj, newk;
  long ilo, ihi, jlo, jhi, klo, khi;
  OctRange *the_rng;

  /* Search the octree for sub-blocks spanning the iso-level.
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
  if(t_level <= the_rng->lo || t_level >= the_rng->hi) {
    return 0;
  }
  if(depth == 0) {
    /* Have reached a leaf of the tree that is cut by the
       iso-surface, so generate triangles.  */
    grab_tris(i, j, k);
    return 1;
  }
  
  /* The contour level cuts this block. Recursively call this
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
        do_blk(newi, newj, newk, child_depth);
      }
    }
  }
  return 1;
}

long grab_tris(long i, long j, long k)
{
  long iSize= t_sizes[0];
  long jSize= t_sizes[1];
  long kSize= t_sizes[2];
  int res;

  /* NOTE: on input i,j,k are indices into the chunk 
     that is currently being processed. 
     Convert into indices into the full grid. */
  i += t_start[0];
  j += t_start[1];
  k += t_start[2];
  if(k <= 0 || k >= kSize-2) {
    /* k is bad - it must be in the interior of the array
       because there are guard cells */
  }
  if(j <= 0 || j >= jSize-2) {
    /* j is bad - it must be in the interior of the array
       because there are guard cells */
  }
  if(i <= 0 || i >= iSize-2) {
    /* i is bad - it must be in the interior of the array
       because there are guard cells */
  }
  /* call the appropriate routine for this grid etc. */
  if(t_ptndx) {
    /* call the routine for point centered data and vertex arrays */
    res= grab_tris_varr(i, j, k);
    return res;
  }
  if(t_vcen) {
    /* call the routine for cell centered data */
    res= grab_tris_zcen(i, j, k);
    return res;
  }
  if(t_xyz) {
    /* call the routine for curvilinear grids */
    res= grab_tris_crv(i, j, k);
    return res;
  }
  /* the default is point centered data on a regular Cartesian grid */
  res= grab_tris_ijk(i, j, k);
  return res;
}

long grab_tris_ijk(long i, long j, long k)
{
  long iSize= t_sizes[0];
  long jSize= t_sizes[1];
  long kSize= t_sizes[2];
  double dx, dy, dz, x0, y0, z0;
  double s[8], vars[8], *var2;
  long ijSize, numTri, *cellIDs;
  int idxvg, idxv, ii, index, mask;
  yPoint3D *xyzverts, *normals, pts[8], gradients[8];
  
  /*
     This cell spans the iso-level.
     Generate the proper triangles using MarchingCubes style
     case tables.
  */
  /* NOTE: on input i,j,k are indices into the full grid. */
  ijSize = iSize * jSize;
  dx= t_deltas[0];
  dy= t_deltas[1];
  dz= t_deltas[2];
  x0= t_origin[0];
  y0= t_origin[1];
  z0= t_origin[2];
  cellIDs = t_triangles->cellIDs;
  xyzverts = t_triangles->xyzverts;
  var2 = t_triangles->var2;
  normals = t_triangles->normals;
  numTri= t_triangles->numTri;
#ifdef DO_STATS
  numcross++;
#endif

  /* get data values at the corners */
  idxvg = i + j*iSize + k*ijSize;
  idxv  = i-t_start[0] + (j-t_start[1])*t_chunk[0] + (k-t_start[2])*t_chunk[0]*t_chunk[1];
  s[0] = t_var[idxvg];
  s[1] = t_var[idxvg+1];
  s[2] = t_var[idxvg+1 + iSize];
  s[3] = t_var[idxvg   + iSize];
  s[4] = t_var[idxvg           + ijSize];
  s[5] = t_var[idxvg+1         + ijSize];
  s[6] = t_var[idxvg+1 + iSize + ijSize];
  s[7] = t_var[idxvg   + iSize + ijSize];

  /* Build the case table */
  for ( ii=0, index = 0, mask= 1; ii < 8; ii++, mask += mask) {
    if ( s[ii] >= t_level ) 
      index |= mask;
  }
  if(t_var2) {
    vars[0] = t_var2[idxvg];
    vars[1] = t_var2[idxvg+1];
    vars[2] = t_var2[idxvg+1 + iSize];
    vars[3] = t_var2[idxvg   + iSize];
    vars[4] = t_var2[idxvg           + ijSize];
    vars[5] = t_var2[idxvg+1         + ijSize];
    vars[6] = t_var2[idxvg+1 + iSize + ijSize];
    vars[7] = t_var2[idxvg   + iSize + ijSize];
  }

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
  
  /* create gradients (i,j,k are global indices) */
  ycPointGradientAll(i, j, k, iSize, jSize, kSize, t_var, dx, dy, dz, gradients);

  extract_tris_tet(index, idxvg, t_var2, &numTri, t_level, s, pts, gradients,
               vars, cellIDs, xyzverts, normals, var2);
  
  /* WARNING: caller must be sure triangle array doesn't overflow. */
  t_triangles->numTri= numTri;
  return 0;
}

long grab_tris_crv(long i, long j, long k)
{
  long iSize= t_sizes[0];
  long jSize= t_sizes[1];
  long kSize= t_sizes[2];
  double s[8], vars[8], *var2;
  long ijSize, numTri, *cellIDs;
  int idxv, idxvg, ii, index, mask;
  yPoint3D *xyzverts, *normals, pts[8], gradients[8];

  cellIDs = t_triangles->cellIDs;
  xyzverts = t_triangles->xyzverts;
  var2 = t_triangles->var2;
  normals = t_triangles->normals;
  numTri= t_triangles->numTri;
#ifdef DO_STATS
  numcross++;
#endif

  /*
     This cell spans the iso-level.
     Generate the proper triangles using MarchingCubes style
     case tables.
  */
  /* NOTE: on input i,j,k are indices into the full grid. */
  ijSize = iSize * jSize;
  /* get data values at the corners */
  idxvg = i + j*iSize + k*ijSize;
  s[0] = t_var[idxvg];
  s[1] = t_var[idxvg+1];
  s[2] = t_var[idxvg+1 + iSize];
  s[3] = t_var[idxvg   + iSize];
  s[4] = t_var[idxvg           + ijSize];
  s[5] = t_var[idxvg+1         + ijSize];
  s[6] = t_var[idxvg+1 + iSize + ijSize];
  s[7] = t_var[idxvg   + iSize + ijSize];

  /* Build the case table */
  for ( ii=0, index = 0, mask= 1; ii < 8; ii++, mask += mask) {
    if ( s[ii] >= t_level ) 
      index |= mask;
  }
  if(t_var2) {
    vars[0] = t_var2[idxvg];
    vars[1] = t_var2[idxvg+1];
    vars[2] = t_var2[idxvg+1 + iSize];
    vars[3] = t_var2[idxvg   + iSize];
    vars[4] = t_var2[idxvg           + ijSize];
    vars[5] = t_var2[idxvg+1         + ijSize];
    vars[6] = t_var2[idxvg+1 + iSize + ijSize];
    vars[7] = t_var2[idxvg   + iSize + ijSize];
  }

  /* set coordinates for the corners of the cell */
  pts[0] = t_xyz[idxvg];
  pts[1] = t_xyz[idxvg+1];
  pts[2] = t_xyz[idxvg+1+iSize];
  pts[3] = t_xyz[idxvg  +iSize];
  pts[4] = t_xyz[idxvg        +ijSize];
  pts[5] = t_xyz[idxvg+1      +ijSize];
  pts[6] = t_xyz[idxvg+1+iSize+ijSize];
  pts[7] = t_xyz[idxvg  +iSize+ijSize];
  
  /* create gradients */
  ycPointGradientCrv(i,   j,   k,   iSize, jSize, kSize, t_xyz, t_var, gradients  );
  ycPointGradientCrv(i+1, j,   k,   iSize, jSize, kSize, t_xyz, t_var, gradients+1);
  ycPointGradientCrv(i+1, j+1, k,   iSize, jSize, kSize, t_xyz, t_var, gradients+2);
  ycPointGradientCrv(i,   j+1, k,   iSize, jSize, kSize, t_xyz, t_var, gradients+3);
  ycPointGradientCrv(i,   j,   k+1, iSize, jSize, kSize, t_xyz, t_var, gradients+4);
  ycPointGradientCrv(i+1, j,   k+1, iSize, jSize, kSize, t_xyz, t_var, gradients+5);
  ycPointGradientCrv(i+1, j+1, k+1, iSize, jSize, kSize, t_xyz, t_var, gradients+6);
  ycPointGradientCrv(i,   j+1, k+1, iSize, jSize, kSize, t_xyz, t_var, gradients+7);

  idxv  = i-t_start[0] + (j-t_start[1])*t_chunk[0] + (k-t_start[2])*t_chunk[0]*t_chunk[1];
  extract_tris_tet(index, idxvg, t_var2, &numTri, t_level, s, pts, gradients,
               vars, cellIDs, xyzverts, normals, var2);
  
  /* WARNING: caller must be sure triangle array doesn't overflow. */
  t_triangles->numTri= numTri;
  return 0;
}

long grab_tris_zcen(long i, long j, long k)
{
  long iSize= t_sizes[0];
  long jSize= t_sizes[1];
  double dx, dy, dz, x0, y0, z0;
  double s[8], vars[8], *var2;
  long ijSize, numTri, *cellIDs;
  int res;
  int idxv, idxvg, ii, index, mask;
  yPoint3D *xyzverts, *normals, pts[8], gradients[8];

  if(t_xyz) {
    /* call the routine for curvilinear grids with cell centered data */
    res= grab_tris_zcen_crv(i, j, k);
	return res;
  }
  dx= t_deltas[0];
  dy= t_deltas[1];
  dz= t_deltas[2];
  x0= t_origin[0];
  y0= t_origin[1];
  z0= t_origin[2];
  cellIDs = t_triangles->cellIDs;
  xyzverts = t_triangles->xyzverts;
  var2 = t_triangles->var2;
  normals = t_triangles->normals;
  numTri= t_triangles->numTri;
#ifdef DO_STATS
  numcross++;
#endif

  /*
     This cell spans the iso-level.
     Generate the proper triangles using MarchingCubes style
     case tables.
  */
  /* NOTE: on input i,j,k are indices into the full grid. */
  ijSize = iSize * jSize;
  /* get data values at the corners */
  idxvg = i + j*iSize + k*ijSize;
  s[0] = t_vcen[idxvg];
  s[1] = t_vcen[idxvg+1];
  s[2] = t_vcen[idxvg+1 + iSize];
  s[3] = t_vcen[idxvg   + iSize];
  s[4] = t_vcen[idxvg           + ijSize];
  s[5] = t_vcen[idxvg+1         + ijSize];
  s[6] = t_vcen[idxvg+1 + iSize + ijSize];
  s[7] = t_vcen[idxvg   + iSize + ijSize];

  /* Build the case table */
  for ( ii=0, index = 0, mask= 1; ii < 8; ii++, mask += mask) {
    if ( s[ii] >= t_level ) 
      index |= mask;
  }
  if(t_var2) {
    /* NOTE: var2 must be point centered by the caller */
    vars[0] = t_var2[idxvg];
    vars[1] = t_var2[idxvg+1];
    vars[2] = t_var2[idxvg+1 + iSize];
    vars[3] = t_var2[idxvg   + iSize];
    vars[4] = t_var2[idxvg           + ijSize];
    vars[5] = t_var2[idxvg+1         + ijSize];
    vars[6] = t_var2[idxvg+1 + iSize + ijSize];
    vars[7] = t_var2[idxvg   + iSize + ijSize];
  }

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
  
  /* create gradients */
  ycPointGradientIntGrdAllZcen(i, j, k, iSize, jSize, dx, dy, dz, t_vcen, gradients);

  idxv  = i-t_start[0] + (j-t_start[1])*t_chunk[0] + (k-t_start[2])*t_chunk[0]*t_chunk[1];
  extract_tris_tet(index, idxvg, t_var2, &numTri, t_level, s, pts, gradients,
               vars, cellIDs, xyzverts, normals, var2);
  
  /* WARNING: caller must be sure triangle array doesn't overflow. */
  t_triangles->numTri= numTri;
  return 0;
}

long grab_tris_zcen_crv(long i, long j, long k)
{
  long iSize= t_sizes[0];
  long jSize= t_sizes[1];
  double s[8], vars[8];
  long ijSize, numTri, *cellIDs;
  int idxv, idxvg, ii, index, mask;
  double *var2;
  yPoint3D *xyzverts, *normals, pts[8], gradients[8];

  cellIDs = t_triangles->cellIDs;
  xyzverts = t_triangles->xyzverts;
  var2 = t_triangles->var2;
  normals = t_triangles->normals;
  numTri= t_triangles->numTri;
#ifdef DO_STATS
  numcross++;
#endif

  /*
     This cell spans the iso-level.
     Generate the proper triangles using MarchingCubes style
     case tables.
  */
  /* NOTE: on input i,j,k are indices into the full grid. */
  ijSize = iSize * jSize;
  /* get data values at the corners */
  idxvg = i + j*iSize + k*ijSize;
  s[0] = t_vcen[idxvg];
  s[1] = t_vcen[idxvg+1];
  s[2] = t_vcen[idxvg+1 + iSize];
  s[3] = t_vcen[idxvg   + iSize];
  s[4] = t_vcen[idxvg           + ijSize];
  s[5] = t_vcen[idxvg+1         + ijSize];
  s[6] = t_vcen[idxvg+1 + iSize + ijSize];
  s[7] = t_vcen[idxvg   + iSize + ijSize];

  /* Build the case table */
  for ( ii=0, index = 0, mask= 1; ii < 8; ii++, mask += mask) {
    if ( s[ii] >= t_level ) 
      index |= mask;
  }
  if(t_var2) {
    /* NOTE: var2 must be point centered by the caller */
    vars[0] = t_var2[idxvg];
    vars[1] = t_var2[idxvg+1];
    vars[2] = t_var2[idxvg+1 + iSize];
    vars[3] = t_var2[idxvg   + iSize];
    vars[4] = t_var2[idxvg           + ijSize];
    vars[5] = t_var2[idxvg+1         + ijSize];
    vars[6] = t_var2[idxvg+1 + iSize + ijSize];
    vars[7] = t_var2[idxvg   + iSize + ijSize];
  }

  /* set coordinates for the corners of the cell */
  pts[0] = t_xyz[idxvg];
  pts[1] = t_xyz[idxvg+1];
  pts[2] = t_xyz[idxvg+1+iSize];
  pts[3] = t_xyz[idxvg  +iSize];
  pts[4] = t_xyz[idxvg        +ijSize];
  pts[5] = t_xyz[idxvg+1      +ijSize];
  pts[6] = t_xyz[idxvg+1+iSize+ijSize];
  pts[7] = t_xyz[idxvg  +iSize+ijSize];
  
  /* create gradients */
  ycPointGradientCrvgAllZcen(idxvg, iSize, jSize, t_xyz+idxvg, t_vcen, gradients);

  idxv  = i-t_start[0] + (j-t_start[1])*t_chunk[0] + (k-t_start[2])*t_chunk[0]*t_chunk[1];
  extract_tris_tet(index, idxvg, t_var2, &numTri, t_level, s, pts, gradients,
               vars, cellIDs, xyzverts, normals, var2);
  
  /* WARNING: caller must be sure triangle array doesn't overflow. */
  t_triangles->numTri= numTri;
  return 0;
}

void ycPointGradientIntGrdAllZcen(long i, long j, long k, long iSize, 
			long jSize, double dx, double dy, double dz, 
			double *var, yPoint3D gradient[8])
{
  long ijSize= iSize*jSize;
  long ibase;
  long ii, jj, kk, pp;
  long ioff[]= {0,1,1,0,0,1,1,0};
  long joff[]= {0,0,1,1,0,0,1,1};
  long koff[]= {0,0,0,0,1,1,1,1};
  double v15, v35, v55, v75;

  /* Approximate the gradient at the vertex from the data at the 
     surrounding zone centers. Guaranteed that the vertex is in 
	 the interior of the mesh. 
     On entry, (i,j,k) is the index of the point. Convert to the 
	 index of the cell to the lower left bottom from this point
	 before accessing var.
  */
  for(pp= 0; pp < 8; pp++) {
    ii= i+ioff[pp];
    jj= j+joff[pp];
    kk= k+koff[pp];
    ibase= ii + jj*iSize + kk*ijSize;
    v15= var[ibase]+var[ibase-1];
    v35= var[ibase-iSize]+var[ibase-1-iSize];
    v55= var[ibase-ijSize]+var[ibase-1-ijSize];
    v75= var[ibase-iSize-ijSize]+var[ibase-1-iSize-ijSize];
    /* x-direction */
    gradient[pp].x = 0.25/dx*(var[ibase]-var[ibase-1]+var[ibase-iSize]-var[ibase-1-iSize]
                  +var[ibase-ijSize]-var[ibase-1-ijSize]
                  +var[ibase-iSize-ijSize]-var[ibase-1-iSize-ijSize]) ;
    /* y-direction */
    gradient[pp].y = 0.25/dy*(v55-v75+v15-v35);
    /* z-direction */
    gradient[pp].z = 0.25/dz*(v15-v55+v35-v75);
  }
}

void ycPointGradientCrvgAllZcen(long idxvg, long iSize, long jSize, 
     yPoint3D *x, double *var, yPoint3D gradient[8])
{
  double dvar, del2;
  long ijSize= iSize*jSize;
  long gbase, pp;
  yPoint3D delta;
  long voff[8];
  voff[0]= 0; voff[1]= 1; voff[2]= 1+iSize; voff[3]= iSize;
  voff[4]= ijSize; voff[5]= 1+ijSize; voff[6]= 1+iSize+ijSize;
  voff[7]= iSize+ijSize;

  /* idxvg is the index of the point at the "lower" corner
     of the cell. */
  for(pp= 0; pp < 8; pp++) {
    gbase= idxvg+voff[pp];
    /* i-direction */
    dvar = var[gbase+1] - var[gbase-1];
    delta.x= x[gbase+1].x-x[gbase-1].x;
    delta.y= x[gbase+1].y-x[gbase-1].y;
    delta.z= x[gbase+1].z-x[gbase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    gradient[pp].x = dvar*delta.x / del2;
    gradient[pp].y = dvar*delta.y / del2;
    gradient[pp].z = dvar*delta.z / del2;

    /* j-direction */
    dvar = var[gbase+iSize] - var[gbase-iSize];
    delta.x= x[gbase+iSize].x-x[gbase-iSize].x;
    delta.y= x[gbase+iSize].y-x[gbase-iSize].y;
    delta.z= x[gbase+iSize].z-x[gbase-iSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    gradient[pp].x += dvar*delta.x / del2;
    gradient[pp].y += dvar*delta.y / del2;
    gradient[pp].z += dvar*delta.z / del2;

    /* k-direction */
    dvar = var[gbase+ijSize] - var[gbase-ijSize];
    delta.x= x[gbase+ijSize].x-x[gbase-ijSize].x;
    delta.y= x[gbase+ijSize].y-x[gbase-ijSize].y;
    delta.z= x[gbase+ijSize].z-x[gbase-ijSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    gradient[pp].x += dvar*delta.x / del2;
    gradient[pp].y += dvar*delta.y / del2;
    gradient[pp].z += dvar*delta.z / del2;
  }
}


/*
   Extract an iso-surface from data on a regular 3D grid.

   ndx is an edge centered scratch array. All elements are initially -1.
   Each time an edge is cut by the iso-surface, ndx is checked for 
   that edge. If it is -1, the vertex and gradient are appended to the
   output lists and the position in the list is stored in ndx. This means
   that the vertices and gradients are only stored once, even though they
   are ordinarily shared by at least 4 triangles.
   Each triangle in the output list is described via three indices into 
   the vertex and gradient lists. 
*/
long grab_tris_varr(long i, long j, long k)
{
  TriVertexGrp *v_triangles= (TriVertexGrp *) t_triangles;
  long iSize= t_sizes[0];
  long jSize= t_sizes[1];
  long kSize= t_sizes[2];
  long ijSize= iSize*jSize;
  double dx, dy, dz, x0, y0, z0;
  double s[8], vars[8], *var2;
  long numEdg, *triStart, *ptndx;
  long numTri, *cellIDs, iOrigin, jOrigin, kOrigin, *nTris, *ndx;
  int idxc, idxv, idxcg, idxvg, ii, case_index, mask;
  yPoint3D *xyzverts, *normals, pts[8], gradients[8];

  iOrigin= t_start[0];
  jOrigin= t_start[1];
  kOrigin= t_start[2];
  dx= t_deltas[0];
  dy= t_deltas[1];
  dz= t_deltas[2];
  x0= t_origin[0];
  y0= t_origin[1];
  z0= t_origin[2];
  cellIDs = v_triangles->cellIDs;
  xyzverts = v_triangles->xyzverts;
  normals = v_triangles->normals;
  ptndx= v_triangles->ptndx;
  var2 = v_triangles->var2;
  numTri= v_triangles->numTri;
  numEdg= v_triangles->numEdg;
  triStart= v_triangles->triStart;
  nTris= v_triangles->nTris;
  ndx= t_ptndx;
#ifdef DO_STATS
  numcross++;
#endif

  /*
     This cell spans the iso-level.
     Generate the proper triangles using MarchingCubes style
     case tables.
  */
  /* NOTE: on input i,j,k are indices into the full grid. */
  ijSize = iSize * jSize;
  /* get data values at the corners. adjust (i,j,k) to find the
     location within the chunk. */
  idxv = i-iOrigin + (j-jOrigin)*t_chunk[0] + (k-kOrigin)*t_chunk[0]*t_chunk[1];
  /* idxvg is the zero-based index of this vertex within the global array. */
  idxvg= i+iSize*j+ijSize*k;
  idxcg= i+(iSize-1)*j+(iSize-1)*(jSize-1)*k;
  if(make_strip) {
    /* record the first triangle for this cell */
    triStart[idxv]= numTri;
    nTris[idxv]= 0;
  }
  s[0] = t_var[idxvg];
  s[1] = t_var[idxvg+1];
  s[2] = t_var[idxvg+1 + iSize];
  s[3] = t_var[idxvg   + iSize];
  s[4] = t_var[idxvg           + ijSize];
  s[5] = t_var[idxvg+1         + ijSize];
  s[6] = t_var[idxvg+1 + iSize + ijSize];
  s[7] = t_var[idxvg   + iSize + ijSize];

  /* Build the case table */
  for ( ii=0, case_index = 0, mask= 1; ii < 8; ii++, mask += mask) {
    if ( s[ii] >= t_level ) 
      case_index |= mask;
  }
  if(t_var2) {
    vars[0] = t_var2[idxvg];
    vars[1] = t_var2[idxvg+1];
    vars[2] = t_var2[idxvg+1 + iSize];
    vars[3] = t_var2[idxvg   + iSize];
    vars[4] = t_var2[idxvg           + ijSize];
    vars[5] = t_var2[idxvg+1         + ijSize];
    vars[6] = t_var2[idxvg+1 + iSize + ijSize];
    vars[7] = t_var2[idxvg   + iSize + ijSize];
  }

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
  
  /* create gradients  (i,j,k are global indices) */
  ycPointGradientAll(i, j, k, iSize, jSize, kSize, t_var, dx, dy, dz, gradients);

  extract_tris_tet_ndx(case_index, idxv, idxcg, t_var2, &numTri, &numEdg, t_level, s, pts, gradients,
                       vars, v_edg_offset, cellIDs, ptndx, ndx, xyzverts, normals, var2);
 
  /* WARNING: caller must be sure triangle array doesn't overflow. */
  v_triangles->numTri= numTri;
  v_triangles->numEdg= numEdg;
  if(numTri > 0) return 1;
  else return 0;
}
