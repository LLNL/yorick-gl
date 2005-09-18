/*
 * $Id: glfunc.h,v 1.1.1.1 2005-09-18 22:07:49 dhmunro Exp $
 * Header file for GL related functions made known to yorick
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __GLFUNC__
#define __GLFUNC__

#ifdef __cplusplus

extern "C" {

#endif

/* need glcode.h for a typedef etc. */
#include "glcode.h"

extern void YError(const char *msg);
extern int (*YputsErr)(char *s);     /* \n appended, like puts */

#ifdef _DEBUG
#define ASSERT(cond, msg) if(!(cond)) YError(msg);
#else
#define ASSERT(cond, msg) 
#endif

#define DEMAND(cond, msg) if(!(cond)) YError(msg);

extern void yglDraw3d(glWinProp *theWin3d);
extern void yglGetCenter(double *center);
extern void yglGetEye(double *eye);
extern void yglGetUp(double *up);
extern void yglGetPixels(long nx, long ny, unsigned char *pix);
extern void yglPutPixels(long nx, long ny, unsigned char *pix);
extern void yglSetGLsize(long w, long h);
extern void yglPolys(long npolys, long *len, float *xyz, float *norm, 
			  float *colr, long edge, long smooth, long do_light);
extern void yglGlyphs(long nglyph, float *origin, float *scal,
			  float *theta, float *phi, float *colr);
extern void yglCells(long nx, long ny, float xyz[3][3], float *norm, 
						 float *colr, long do_alpha);
extern void yglMap2ColorRaw3d(long ncolr, unsigned char *red,
                            unsigned char *green, 
                            unsigned char *blue,
                            double vmin, double vmax, double *var, 
                            long ntri, long *cellids, float *colr);
extern double yglGetVers3d(void);
extern int  yglQueryTex3d(glWinProp *theWin3d);
extern void yglPrepTex2d(void);
extern void yglEndTex2d(void);
extern void yglPrepTex3d(void);
extern void yglEndTex3d(void);
extern void yglPrepCubeTex(void);
extern void yglEndCubeTex(void);
extern void yglLdCubeTex(void);
extern int  yglQueryTexCube(void);
extern void yglTexPoly(long nvert, float *verts, float *texverts);
extern void yglLdTex3d(long nx, long ny, long nz, unsigned char *tex);
extern void yglSetTexPalette(unsigned char *texpal);
extern void yglTexTris(long ntri, float *verts, float *texverts);
extern void yglTexcells(long nx, long ny, long nz, double delta[3], 
				 unsigned char *data, unsigned char *ctab);
extern void yglTexcell2(long nx, long ny, long nz, double delta[3], 
				 unsigned char *texdata);
extern void yglTex3dbox(double ds, double *origin, double *len);
extern void yglPlm(long nx, long ny, float *xyz, float *colr);
extern void yglPlf(long nx, long ny, float *xyz, float *colr);
extern void yglSurf(long do_alpha, long nx, long ny, float *xyz, 
                      float *norm, float *colr);
extern void yglColrSurf(long do_alpha, long nx, long ny, float *xyz, 
                          float *norm, float *colr);
extern void yglLines(long nvert, float *xyz, float *colr);
extern void yglPoints(long nvert, float *xyz, float *colr);
extern void yglLineWidth(double width);
extern void yglPrepContext(void);
extern void yglSetPolyMode(int mode);
extern void yglSetPolySides(int sides);
extern float yglGetMatSpec(void);
extern void yglSetMatSpec(float spec);
extern void yglSetShade(int smooth);
extern void yglSetColorType(int colorType);
extern void yglSetLight(double ambient, double diffuse, double specular, 
						double spower, double *sdir);
extern void yglFrontFace3d(long dir);
extern void yglTivarray(long ntri, unsigned int *ptndx, void *ileave);
extern void yglTvarray(long do_alpha, long cpervrt, long ntri, unsigned int *ptndx,
                       float *xyz, float *norm, float *colr);

#ifdef __cplusplus
	}
#endif

#endif /* Include/Define */
