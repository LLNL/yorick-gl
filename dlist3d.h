/*
 * $Id: dlist3d.h,v 1.1.1.1 2005-09-18 22:07:44 dhmunro Exp $
 * Declare functions used for manipulating 3D display lists.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef DLIST3D_H
#define DLIST3D_H

typedef struct yBox3D yBox3D;
struct yBox3D {
  double xmin, xmax, ymin, ymax, zmin, zmax;
} ;

typedef void yglDrawFunc3d(int mode, void *data);

typedef struct yList3d_Elem yList3d_Elem;
struct yList3d_Elem {
  yBox3D box;
  yglDrawFunc3d *func;
  void *data;
  yList3d_Elem *next;
} ;

typedef struct yPoly3dData yPoly3dData;
struct yPoly3dData {
  long npolys, edge, smooth, do_light;
  long *len;
  float *xyz, *norm, *colr;
} ;

typedef struct yGlyph3dData yGlyph3dData;
struct yGlyph3dData {
  long nglyph;
  float *origin, *scal, *theta, *phi, *colr;
} ;

typedef struct yCell3dData yCell3dData;
struct yCell3dData {
  long nx, ny, do_alpha;
  float *xyz, *norm, *colr;
} ;

typedef struct yPlm3dData yPlm3dData;
struct yPlm3dData {
  long nx, ny;
  float *xyz, *colr;
} ;

typedef struct yPlf3dData yPlf3dData;
struct yPlf3dData {
  long nx, ny;
  float *xyz, *colr;
} ;

typedef struct ySurf3dData ySurf3dData;
struct ySurf3dData {
  long do_alpha, nx, ny;
  float *xyz, *norm, *colr;
} ;

typedef struct yColrSurf3dData yColrSurf3dData;
struct yColrSurf3dData {
  long do_alpha, nx, ny;
  float *xyz, *norm, *colr;
} ;

typedef struct yTstrips3dData yTstrips3dData;
struct yTstrips3dData {
  long nstrips, edge, smooth, do_light, do_alpha;
  long *len;
  float *xyz, *norm, *colr;
} ;

typedef struct yQstrips3dData yQstrips3dData;
struct yQstrips3dData {
  long nstrips, edge, smooth, do_light, do_alpha;
  long *len;
  float *xyz, *norm, *colr;
} ;

typedef struct yTstripsNdx3dData yTstripsNdx3dData;
struct yTstripsNdx3dData {
  long nstrips, ntri, nvert, numedg, edge, do_alpha;
  long *len, *ndx;
  float *xyz, *norm, *colr;
} ;

typedef struct yTarray3dData yTarray3dData;
struct yTarray3dData {
  long ntri, edge, smooth, do_light, do_alpha, cpervrt, cubemap, emit;
  float *xyz, *norm, *colr;
} ;

typedef struct yQarray3dData yQarray3dData;
struct yQarray3dData {
  long nquad, edge, smooth, do_light, do_alpha, cpervrt;
  float *xyz, *norm, *colr;
} ;

typedef struct yTivarray3dData yTivarray3dData;
struct yTivarray3dData {
  long ntri, nvert;
  unsigned int *ptndx;
  float *ileave;
} ;

typedef struct yTvarray3dData yTvarray3dData;
struct yTvarray3dData {
  long ntri, nvert, cpervrt, do_alpha;
  unsigned int *ptndx;
  float *xyz, *norm, *colr;
} ;

typedef struct yLines3dData yLines3dData;
struct yLines3dData {
  long nvert;
  float *xyz, *colr;
} ;

typedef struct yPoints3dData yPoints3dData;
struct yPoints3dData {
  long nvert;
  float *xyz, *colr;
} ;

typedef struct yTex3dData yTex3dData;
struct yTex3dData {
  double ds, *origin, *boxsiz;
} ;

typedef struct yTexcell2dData yTexcell2dData;
struct yTexcell2dData {
  long nx, ny, nz;
  double *znsiz;
  unsigned char *texval;
} ;

typedef struct yPix3dData yPix3dData;
struct yPix3dData {
  long nx, ny;
  unsigned char *pix;
} ;

extern yList3d_Elem *yglNewCachedList3dElem(void);
extern yList3d_Elem *yglNewDirectList3dElem(void);
extern void yglDrawCurr3d(void);
extern void yglClearList3d(void);
extern void yglClearCachedList3d(void);
extern void yglClearDirectList3d(void);
extern void yglDrawList3d(void);
extern void yglDrawListCache3d(void);
extern void yglDrawListDirect3d(void);
extern long yglGetBounds3d(yBox3D *box);
extern long yglGetBoundsDirectList3d(yBox3D *box);
extern long yglGetBoundsCachedList3d(yBox3D *box);
extern void yglSetLims3d(yList3d_Elem *elem, long nvert, float *xyz);
extern void yglDrawGnomon(void);
extern void yglDrawCage(void);
extern void yglPolys3d(long npolys, long *len, double *xyz, double *norm,
            double *colr, long edge, long smooth, long do_light);
extern void yglGlyphs3d(long nglyph, double *origin, double *scal,
			double *theta, double *phi, double *colr);
extern void yglCells3d(long nx, long ny, double *corners,
			double *norm, double *colr, long do_alpha);
extern void yglPlm3d(long nx, long ny, double *xyz, double *colr);
extern void yglPlf3d(long nx, long ny, double *xyz, double *colr);
extern void yglSurf3d(long do_alpha, long nx, long ny, double *xyz, 
            double *norm, double *colr);
extern void yglColrsurf3d(long do_alpha, long nx, long ny, double *xyz, 
            double *norm, double *colr);
extern void yglLines3d(long nvert, double *xyz, double *colr);
extern void yglPoints3d(long nvert, double *xyz, double *colr);
extern void yglTstrips3d(long nstrips, long *len, double *xyz,
            double *norm, double *colr, long edge,
            long smooth, long do_light, long do_alpha);
extern void yglQstrips3d(long nstrips, long *len, double *xyz,
            double *norm, double *colr, long edge,
            long smooth, long do_light, long do_alpha);
extern void yglTstripsndx3d(long nstrips, long numedg, long ntri,
            long *len, long *ndx, double *xyz,
            double *norm, double *colr, long edge, long do_alpha);
extern void yglTarray3d(long ntri, double *xyz, double *norm, double *colr, 
            long edge, long smooth, long do_light, long do_alpha,
            long cpervrt, long cubemap, long emit);
extern void yglQarray3d(long nquad, double *xyz, double *norm, double *colr, 
            long edge, long smooth, long do_light, long do_alpha,
            long cpervrt);
extern void yglTivarray3d(long ntri, long nvert, long *ptndx, double *xyz, 
            double *norm, double *colr);
extern void yglTvarray3d(long ntri, long nvert, long do_alpha, long cpervrt,
                         long *ptndx, double *xyz, double *norm, double *colr);
extern void yglTex3d(float ds, double *origin, double *boxsiz);
extern void yglTexcell2d(long nx, long ny, long nz, double *znsiz, char *texval);
extern void yglPlpix3d(long nx, long ny, char *pix);

extern yglDrawFunc3d yglDrawPolys3d;
extern yglDrawFunc3d yglDrawGlyphs3d;
extern yglDrawFunc3d yglDrawCells3d;
extern yglDrawFunc3d yglDrawPlm3d;
extern yglDrawFunc3d yglDrawPlf3d;
extern yglDrawFunc3d yglDrawSurf3d;
extern yglDrawFunc3d yglDrawColrSurf3d;
extern yglDrawFunc3d yglDrawLines3d;
extern yglDrawFunc3d yglDrawPoints3d;
extern yglDrawFunc3d yglDrawTstrips3d;
extern yglDrawFunc3d yglDrawQstrips3d;
extern yglDrawFunc3d yglDrawTstripsNdx3d;
extern yglDrawFunc3d yglDrawTarray3d;
extern yglDrawFunc3d yglDrawQarray3d;
extern yglDrawFunc3d yglDrawTivarray3d;
extern yglDrawFunc3d yglDrawTvarray3d;
extern yglDrawFunc3d yglDrawTex3d;
extern yglDrawFunc3d yglDrawTexcell2d;
extern yglDrawFunc3d yglDrawPix3d;

extern yList3d_Elem *yListDirectHead;
extern yList3d_Elem *yListCachedHead;
extern int yDrawBBox3d;
#endif
