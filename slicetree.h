/*
 * $Id: slicetree.h,v 1.1 2005-09-18 22:08:03 dhmunro Exp $
 * Header file for octree used by slicing plane.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __SLICETREE__
#define __SLICETREE__

#ifdef __cplusplus

extern "C" {

#endif

typedef struct OctSpan OctSpan;
struct OctSpan {
  double xmin, xmax, ymin, ymax, zmin, zmax;
} ;

typedef struct OctSTree OctSTree;
struct OctSTree {
  long maxdepth;
  long *start;
  long *chunk;
  long *size;
  long *trsiz;
  long *offsets;
  OctSpan *ranges;
  OctSTree *next;
} ;

extern int ycMakeSliceTreeCrv(double *xyz, OctSTree *tree);
extern void firstSblk(long *start, long *sizes, long *trsiz, 
                      double *xyz, OctSpan *rng);
extern void nextSblk(long *trsiz, OctSpan *oldr, OctSpan *nrng);
extern int ycSliceTree(long maxdepth, long sizes[3], long chunk[3], long start[3], 
                       double deltas[3], double origin[3], double point[3], 
                       double normal[3], double *var, TriArrayGrp *triangles);
extern int ycSliceTreeCrv(double point[3], double normal[3], double *xyz, 
                          double *var, TriArrayGrp *triangles, OctSTree *tree);
extern long do_Sblk(long ilo, long ihi, long jlo, long jhi, long klo, 
			 long khi, long depth);
extern long do_SblkCrv(long newi, long newj, long newk, long depth);
extern long grab_Stris(long i, long j, long k, double s[8]);
extern long grab_StrisCrv(long i, long j, long k);

#ifdef __cplusplus
}
#endif

#endif /* Include/Define */
