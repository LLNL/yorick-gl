/*
 * $Id: glWrappers.h,v 1.1.1.1 2005-09-18 22:07:48 dhmunro Exp $
 * Header file for the routines that handle mouse motion
 * in OpenGL windows.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __GLWRAPPERS__
#define __GLWRAPPERS__

#include "TriStruct.h"

#ifdef __cplusplus

extern "C" {

#endif

extern void yglInitLights3d(void);

extern void yglGrabPix3d(int nArgs);
extern float yglGetFov3d(void);
extern void yglSetFov3d(float fov);
extern void yglGetCenter3d(double *center);
extern void yglGetEye3d(double *eye);
extern void yglGetUp3d(double *up);
extern long yglGetWidth3d(void);
extern long yglGetHite3d(void);
extern int yglHasTex3d(void);
extern int yglHasTexcube3d(void);
extern void yglAlwaysShowObj3d(long flag);
extern void yglUseList3d(long flag);
extern void yglUseArray3d(long flag);
extern void yglTexTtris3d(int nArgs);
extern void yglLdtex3d(int nArgs);
extern void yglTexcells3d(int nArgs);
extern void yglTstrips_alpha3d(int nArgs);
extern void yglQstrips_alpha3d(int nArgs);
extern void yglMsmovVal3d(double val);
extern void yglOutCcw3d(long dir);
extern void yglIncSeq3d(void);
extern void yglArsum3d(long nx0, long ny0, long nz0, long nsx, long nsy, long nsz, 
            double *vin, double *vout);
extern long yglSizeTriArraysNdx3d(TriVertexGrp *list);
extern void yglCopyTriArrayNdx3d(TriVertexGrp *list, TriVertexGrp *nlist);
extern void yglCollapseTriArraysNdx3d(long colrtyp, TriVertexGrp *list, 
            TriVertexGrp *nlist);
extern void yglDoSortTriNdx3d(TriVertexGrp *oldtri, long *newptndx);
extern long yglSizeTriArrays3d(TriArrayGrp *list);
extern void yglCopyTriArray3d(long numtri, TriArrayGrp *list, TriArrayGrp *nlist);
extern void yglCollapseTriArrays3d(long colrtyp, TriArrayGrp *list, TriArrayGrp *nlist);
extern void yglDoSortTri3d(long colrtyp, TriArrayGrp *oldtri, TriArrayGrp *newtri);
extern void yglSliceTris3d(long *keep, long *nkeep, double *dp,
            TriArrayGrp *oldtri, TriArrayGrp *newtri);
extern void yglMouseFunc3d(long val);

#ifdef __cplusplus
	}
#endif

#endif /* Include/Define */
