/*
 * $Id: tristruct.i,v 1.1 2005-09-18 22:08:03 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
/* NOTE: these should not be exposed to interpreted interface */

/* data structures for triangle arrays and triangle strips */

func ASSERT(cond, msg)
{
  if(!cond) error(msg);
}

struct OctSpan {
  double xmin, xmax, ymin, ymax, zmin, zmax;
}

struct OctSTree {
  long maxdepth;
  pointer start;
  pointer chunk;
  pointer size;
  pointer trsiz;
  pointer offsets;
  pointer ranges;
  pointer next;
}

struct OctRange {
  double lo, hi;
}

struct OctTree {
  long maxdepth;
  pointer start;
  pointer chunk;
  pointer size;
  pointer trsiz;
  pointer offsets;
  pointer ranges;
  OctTree *next;
}

struct TriArray{
  long numTri;
  pointer cellIDs;
  pointer xyzverts;
  pointer normals;
  long nxti;
  long nxtj;
  long nxtk;
  TriArray *next;
}

struct TriArrayGrp {
  long numTri;
  pointer cellIDs;
  pointer xyzverts;
  pointer normals;
  pointer var2;
  pointer colors;
  pointer triEdg;
  pointer triStart;
  pointer nTris;
  TriArrayGrp *next;
}

struct PolyArrayGrp {
  long numPoly;
  pointer edges;
  pointer cellIDs;
  pointer xyzverts;
  pointer normals;
  pointer var2;
  pointer colors;
  pointer triEdg;
  pointer triStart;
  pointer nVerts;
  PolyArrayGrp *next;
}

struct TriStrip {
  long nVert;
  long nStrip;
  pointer cellIDs;
  pointer xyzverts;
  pointer normals;
  pointer triLen;
  long nxti;
  long nxtj;
  long nxtk;
  TriStrip *next;
}

struct TriStripGrp {
  long nVert;
  long nStrip;
  pointer cellIDs;
  pointer xyzverts;
  pointer normals;
  pointer colors;
  pointer triLen;
  TriStripGrp *next;
}

struct TriStripNdx {
  long nVert;
  long nStrip;
  pointer cellIDs;
  pointer ptndx;
  pointer triLen;
}

struct TriStripNdxGrp {
  long nVert;
  long nStrip;
  long numEdg;
  long numTri;
  pointer cellIDs;
  pointer ptndx;
  pointer xyzverts;
  pointer normals;
  pointer colors;
  pointer triLen;
  TriStripNdxGrp *next;
}

struct TriVertexGrp {
  long numTri;
  long numEdg;
  pointer cellIDs;
  pointer xyzverts;
  pointer normals;
  pointer colors;
  pointer var2;
  pointer ptndx;
  pointer triEdg;
  pointer triStart;
  pointer nTris;
  TriVertexGrp *next;
}

struct IsoTriStrip {
  long nvert;
  long edges(12);
}

struct IsoTriStripArray {
  long nStrip;
  IsoTriStrip the_tris(6);
}
