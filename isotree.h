/*
 * $Id: isotree.h,v 1.1 2005-09-18 22:08:02 dhmunro Exp $
 * Header file for octree used by iso-surface.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __ISOTREE__
#define __ISOTREE__

#ifdef __cplusplus

extern "C" {

#endif

typedef struct OctRange OctRange;
struct OctRange {
  double lo, hi;
} ;

typedef struct OctTree OctTree;
struct OctTree {
  long maxdepth;
  long *start;
  long *chunk;
  long *size;
  long *trsiz;
  long *offsets;
  OctRange *ranges;
  OctTree *next;
} ;

extern int ycMakeContourTree(double *var, OctTree *tree);
extern void firstblk(double *data, long *start, long *sizes, long *trsiz, OctRange *rng);
extern void nextblk(long *trsiz, OctRange *oldr, OctRange *nrng);
extern int ycContourTree(double deltas[3], double origin[3], double level,
                         double *var, TriArrayGrp *triangles, OctTree *tree);
extern int ycContourTree2(double deltas[3], double origin[3], double level,
                          double *var, double *var2, TriArrayGrp *triangles,
                          OctTree *tree);
extern int ycContourTreeCrv(double level, yPoint3D *xyz, double *var,
                            TriArrayGrp *triangles, OctTree *tree);
extern int ycContourTreeCrv2(double level, yPoint3D *xyz, double *var, 
                             double *var2, TriArrayGrp *triangles, OctTree *tree);
extern int ycContourTreeZcen(double deltas[3], double origin[3], double level,
                             double *var, double *vcen, TriArrayGrp *triangles,
                             OctTree *tree);
extern int ycContourTreeZcen2(double deltas[3], double origin[3], double level,
                              double *var, double *vcen, double *var2, TriArrayGrp *triangles,
                              OctTree *tree);
extern int ycContourTreeCrvZcen(double level, yPoint3D *xyz, double *var,
                                double *vcen, TriArrayGrp *triangles, OctTree *tree);
extern int ycContourTreeCrvZcen2(double level, yPoint3D *xyz, double *var,
                                 double *vcen, double *var2, TriArrayGrp *triangles, OctTree *tree);
extern int ycContourTreeVarr(double deltas[3], double origin[3], double level,
                             double *var, TriVertexGrp *triangles, OctTree *tree, long *edgndx);
extern int ycContourTreeVarr2(double deltas[3], double origin[3], double level,
                              double *var, double *var2, TriVertexGrp *triangles, OctTree *tree,
                              long *edgndx);
extern int ycContourTreeVarrStr(double deltas[3], double origin[3], double level,
                                double *var, TriVertexGrp *triangles, OctTree *tree, long *edgndx);

extern long do_blk(long i, long j, long k, long depth);
extern long grab_tris(long i, long j, long k);
extern long grab_tris_ijk(long i, long j, long k);
extern long grab_tris_crv(long i, long j, long k);
extern long grab_tris_zcen(long i, long j, long k);
extern long grab_tris_zcen_crv(long i, long j, long k);
extern long grab_tris_varr(long i, long j, long k);
extern void ycPointGradientIntGrdAllZcen(long i, long j, long k, long iSize, 
			long jSize, double dx, double dy, double dz, 
			double *var, yPoint3D gradient[8]);
extern void ycPointGradientCrvgAllZcen(long idxg, long iSize, long jSize, 
                        yPoint3D *x, double *var, yPoint3D gradient[8]);

#ifdef __cplusplus
}
#endif

#endif /* Include/Define */
