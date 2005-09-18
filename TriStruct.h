/*
 * $Id: TriStruct.h,v 1.1.1.1 2005-09-18 22:07:58 dhmunro Exp $
 * TriStruct - defines triangle array and strip structures
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __TriStruct_h
#define __TriStruct_h


typedef struct yPoint3D yPoint3D;
struct yPoint3D {
  double x;
  double y;
  double z;
} ;

typedef struct TriArray TriArray;
struct TriArray {
  long numTri;
  long *cellIDs;
  yPoint3D *xyzverts;
  yPoint3D *normals;
  long nxti;
  long nxtj;
  long nxtk;
  TriArray *next;
} ;

typedef struct TriStrip TriStrip;
struct TriStrip {
  long nVert;
  long nStrip;
  long *cellIDs;
  yPoint3D *xyzverts;
  yPoint3D *normals;
  long *triLen;
  long nxti;
  long nxtj;
  long nxtk;
  TriStrip *next;
} ;

typedef struct TriArrayGrp TriArrayGrp;
struct TriArrayGrp {
  long numTri;
  long *cellIDs;
  yPoint3D *xyzverts;
  yPoint3D *normals;
  double *var2;
  float *colors;
  long *triEdg;
  long *triStart;
  long *nTris;
  TriArrayGrp *next;
} ;

typedef struct TriStripGrp TriStripGrp;
struct TriStripGrp {
  long nVert;
  long nStrip;
  long *cellIDs;
  yPoint3D *xyzverts;
  yPoint3D *normals;
  float *colors;
  long *triLen;
  TriStripGrp *next;
} ;

typedef struct TriStripNdx TriStripNdx;
struct TriStripNdx {
  long nVert;
  long nStrip;
  long *cellIDs;
  long *ndx;
  long *triLen;
} ;

typedef struct TriStripNdxGrp TriStripNdxGrp;
struct TriStripNdxGrp {
  long nVert;
  long nStrip;
  long numEdg;
  long numTri;
  long *cellIDs;
  long *ndx;
  yPoint3D *xyzverts;
  yPoint3D *normals;
  float *colors;
  long *triLen;
  TriStripNdxGrp *next;
} ;

typedef struct TriVertexGrp TriVertexGrp;
struct TriVertexGrp {
  long numTri;
  long numEdg;
  long *cellIDs;
  yPoint3D *xyzverts;
  yPoint3D *normals;
  float *colors;
  double *var2;
  long *ptndx;
  long *triEdg;
  long *triStart;
  long *nTris;
  TriVertexGrp *next;
} ;

typedef struct PolyArrayGrp PolyArrayGrp;
struct PolyArrayGrp {
  long numPoly;
  long *edges;
  long *cellIDs;
  yPoint3D *xyzverts;
  yPoint3D *normals;
  double *var2;
  float *colors;
  long *triEdg;
  long *triStart;
  long *nVerts;
  PolyArrayGrp *next;
} ;

#endif
