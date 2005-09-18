/*
 * $Id: Contour3D.h,v 1.1 2005-09-18 22:07:55 dhmunro Exp $
 * Contour3D - generates isosurface(s) from volume data
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __Contour3D_h
#define __Contour3D_h

#include "TriStruct.h"


/* WARNING: the caller is responsible for freeing the storage pointed
   in the struct returned by this call, and for freeing the struct. */
extern int ycContour(long sizes[3], long maxTri, double deltas[3], double origin[3], 
			  double level, double *var, TriArray *triangles);
extern int ycContourCrv(long sizes[3], long maxTri, double level, yPoint3D *xyz, 
			  double *var, TriArray *triangles);
extern int ycContourList(long sizes[3], long maxVert, double deltas[3],
                         double origin[3], double level, double *var,
                         TriStrip *triangles);
extern int ycContourListCrv(long sizes[3], long maxVert, double level, yPoint3D *xyz,
                            double *var, TriStrip *triangles);
extern int ycContourHex(long nzone, long maxTri, double level, yPoint3D *xyz,
                        double *var, TriArray *triangles);
extern int ycContourStr(long make_strip,long sizes[3], long offsets[5], double deltas[3], 
			  double origin[3], double level, double *var, 
			  yPoint3D *grd, char *done, unsigned char *above,
			  TriArrayGrp *triangles);
extern int ycContourStr_old(long sizes[3], long start[3], long chunk[3], 
			  double deltas[3], double origin[3], 
			  double level, double *var, TriArrayGrp *triangles);
extern int ycContourCrvStr(long make_strip, long sizes[3], long offsets[5], double level, yPoint3D *xyz,
                           double *var, yPoint3D *grd, char *done,
                           unsigned char *above, TriArrayGrp *triangles);
extern int ycContourListStr(long make_strip, long sizes[3], long maxVert, double deltas[3],
                            double origin[3], double level, double *var,
                            TriStripGrp *triangles);
extern int ycContourListCrvStr(long make_strip, long sizes[3], long maxVert, double level, yPoint3D *xyz,
                               double *var, TriStripGrp *triangles);
extern int ycContourHexStr(long nzone, long maxTri, double level, yPoint3D *xyz, 
			  double *var, TriArrayGrp *triangles);
extern int ycContourStrCen(long make_strip, long sizes[3], long offsets[7], 
			  double deltas[3], double origin[3], double level, double *var, 
			  double *vcen, yPoint3D *grd, char *done, unsigned char *above,
			  TriArrayGrp *triangles);
extern int ycContourStrPt(long make_strip, long sizes[3], long offsets[7], 
			  double deltas[3], double origin[3], double level, double *var, 
			  yPoint3D *grd, char *done, unsigned char *above,
			  TriArrayGrp *triangles);
extern int ycContourPt(long make_strip, long sizes[3], long offsets[7], 
				 double deltas[3], double origin[3], double level, 
				 double *var, TriArrayGrp *triangles);
extern int ycContourPtCrv(long make_strip, long sizes[3], long offsets[5],
                          double level, yPoint3D *xyz, double *var,
                          TriArrayGrp *triangles);
extern int ycContourVarrCen(long make_strip, long sizes[3], long offsets[7], 
			  double deltas[3], double origin[3], double level, double *var, 
			  double *vcen, yPoint3D *grd, char *done, unsigned char *above,
			  long* ndx, TriVertexGrp *triangles);
extern int ycContourHexList(long base, long nzone, long npt, double level, long *zones, 
                 yPoint3D *xyz, double *var, yPoint3D *norm, 
                 TriArray *triangles);

extern int ycMakeTriStripNdx(long sizes[3], long offsets[5], char *visited,
                             TriArrayGrp *triangles, TriStripNdx *strip);
extern int ycMakeTriStripNdx_old(long sizes[3], long start[3], long chunk[3], 
				 char *visited, TriArrayGrp *triangles, TriStripNdx *strip);
extern int ycMakeTriStrip(long sizes[3], long start[3], long chunk[3],
                          char *visited, TriArrayGrp *triangles, TriStripGrp *strip);
extern int ycMakeTriStripNdxCenVertex(long sizes[3], long offsets[5], char *visited,
                                      TriVertexGrp *triangles, TriStripNdx *strip);
extern int ycMakeTriStripNdxCen(long sizes[3], long offsets[5], char *visited,
                                TriArrayGrp *triangles, TriStripNdx *strip);
extern int ycAssembleStrip(TriArrayGrp *triangles, TriStripNdx *stripNdx, 
				 TriStripGrp *strip);

extern void ycPointGradient(long i, long j, long k, long nx, long ny, long nz, double *s, 
                      double dx, double dy, double dz, yPoint3D *n);
extern void ycPointGradientAll(long i, long j, long k, long nx, long ny, long nz, double *s, 
                   double dx, double dy, double dz, yPoint3D gradient[8]);
extern void ycPointGradientCrv(long i, long j, long k, long nx, long ny, long nz, 
					   yPoint3D *x, double *s, yPoint3D *n);
extern void ycPointGradientHex(yPoint3D *pts, double *s, yPoint3D *gradients);
extern void ycPointGradientGrd(long i, long j, long k, long nx, long ny, long nz,
                               double *s, double dx, double dy, double dz, yPoint3D *n,
                               yPoint3D *grd, char *done);
extern void ycPointGradientAllGrd(long i, long j, long k, long nx, long ny, long nz,
                                  double *s, double dx, double dy, double dz,
                                  yPoint3D gradient[8], yPoint3D *grd, char *done);
extern void ycPointGradientCrvGrd(long i, long j, long k, long nx, long ny, long nz,
                                  yPoint3D *x, double *s, yPoint3D *n, yPoint3D *grd, char *done);
extern void ycPointGradientIntGrd(long i, long j, long k, long nx, long ny, long nz, double *s,
                                  double dx, double dy, double dz, yPoint3D gradient[8],
                                  yPoint3D *grd, char *done);
extern void ycPointGradientIntGrdCrv(long i, long j, long k, long nx, long ny, long nz, 
				   yPoint3D *xyz, double *s, yPoint3D gradient[8], 
				   yPoint3D *grd, char *done);
extern void ycPointGradientIntGrdCenAll(long i, long j, long k, long nx, 
			long ny, long nz, double dx, double dy, 
			double dz, yPoint3D *grd, char *done, long offsets[5],
			double *var, yPoint3D gradient[8]);
extern void ycPointGradientIntGrdPtAll(long i, long j, long k, long nx, 
			long ny, long nz, double dx, double dy, 
			double dz, yPoint3D *grd, char *done, long offsets[5],
			double *var, yPoint3D gradient[8]);
extern void ycPointGradientIntGrdPtAll2(long idxg, double dx, double dy, 
			double dz, long corners[8], long offsets[5], double *var, 
			yPoint3D gradient[8]);
extern void ycPointGradientCrvgPtAll(long idxg, long corners[8], long offsets[5],
                                     yPoint3D *x, double *var, yPoint3D gradient[8]);
extern void ycNormalize(yPoint3D *n);

extern long GetNeighbor(long sizes[3], long cellID, long prevEdg, long lastEdg,
                        long *newPrev, long *newLast);

extern void extract_tris_tet(int case_index, long idxcg, double *cntr_v2, long *numTri,
                    double lev, double s[8], yPoint3D pts[8], yPoint3D gradients[8],
                    double vars[8], long *cellIDs, yPoint3D *xyzverts,
                    yPoint3D *normals, double *v2vals);
extern void extract_tris_tet_ndx(int case_index, long idxv, long idxcg, double *cntr_v2, long *numTri,
                    long *numEdg, double lev, double s[8], yPoint3D pts[8],
                    yPoint3D gradients[8], double vars[8], long edg_offset[12], long *cellIDs,
                    long *ptndx, long *ndx, yPoint3D *xyzverts, yPoint3D *normals, double *v2vals);
extern void extract_slicetris_tet(int case_index, long idxcg, double *cntr_v2, long *numTri,
                  double s[8], yPoint3D pts[8], double vars[8],
                  long *cellIDs, yPoint3D *xyzverts, double *v2vals);
extern int ycPrepIsoTet(void);

extern int have_iso_cases;

#endif


