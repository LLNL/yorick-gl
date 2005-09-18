/*
 * $Id: glGlyph.c,v 1.1.1.1 2005-09-18 22:07:45 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "glcode.h"
#include "glfunc.h"
#include "TriStruct.h"
#include <math.h>
#include <stdlib.h>

/*
#define CHEK_ERROR(x)	gl_chek_error(x)
*/
#define CHEK_ERROR(x)

typedef struct glyph glyph;
struct glyph {
  long ntri;
  yPoint3D *xyz;
  yPoint3D *nrm;
} ;

#define HT_RAT 3.0f
#define TET_DFLT

static glyph *c_glyph= 0;

static glyph pyramid_glyph;
static glyph tet_glyph;

extern void makPyrGlyph(void);
extern void makTetGlyph(void);

void yglGlyphs(long nglyph, float *origins, float *scal,
                 float *theta, float *phi, float *colr)
{
  int ng, nv;
  float x, y, z, sc, nrmx, nrmy, nrmz;
  float th, ph, snth, csth, snph, csph, csph_csth, snph_csth, csph_snth, snph_snth;
  yPoint3D xyz0, xyz, nrm;

  /* The input is a list of glyphs (ellipsoids). 
     Each glyph has a center, a scale factor, a color,
     and a direction (theta and phi).
  */

  /* draw the glyph list now */
  if(nglyph <= 0) return;

  if(!c_glyph) {
#ifdef TET_DFLT
    makTetGlyph();
    c_glyph= &tet_glyph;
#else
    makPyrGlyph();
    c_glyph= &pyramid_glyph;
#endif
  }

  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();
  CHEK_ERROR("yglGlyphs setup");
  for(ng= 0; ng < nglyph; ng++) {
    xyz0.x= origins[3*ng];
    xyz0.y= origins[3*ng+1];
    xyz0.z= origins[3*ng+2];
    sc= scal[ng];
    th= theta[ng];
    snth= (float) sin(th);
    csth= (float) cos(th);
    ph= phi[ng];
    snph= (float) sin(ph);
    csph= (float) cos(ph);
    csph_csth= csph*csth;
    snph_csth= snph*csth;
    csph_snth= csph*snth;
    snph_snth= snph*snth;
    glColor3fv(colr+3*ng);
    glBegin(GL_TRIANGLES);
    for(nv= 0; nv < 3*c_glyph->ntri; nv++) {
      xyz= (c_glyph->xyz)[nv];
      nrm= (c_glyph->nrm)[nv];
      x= (float) (xyz0.x+csph_csth*xyz.x*sc+snph_csth*xyz.y*sc-snth*xyz.z*sc);
      y= (float) (xyz0.y-snph*xyz.x*sc+csph*xyz.y*sc);
      z= (float) (xyz0.z+csph_snth*xyz.x*sc+snph_snth*xyz.y*sc+csth*xyz.z*sc);
      nrmx= (float) (csph_csth*nrm.x+snph_csth*nrm.y-snth*nrm.z);
      nrmy= (float) (snph*nrm.x+csph*nrm.y);
      nrmz= (float) (csph_snth*nrm.x+snph_snth*nrm.y+csth*nrm.z);
      glNormal3f(nrmx, nrmy, nrmz);
      glVertex3f(x, y, z);
    }
    glEnd();
  }
  CHEK_ERROR("yglGlyphs");
}

void makPyrGlyph(void)
{
  int nv;
  yPoint3D *xyz, *nrm;
  float nrz, nrt;

  /* The pyramid is made up of 6 triangles.
     compute and save the coordinates of their vertices. */
  pyramid_glyph.ntri= 6;
  xyz= (yPoint3D *) malloc(sizeof(yPoint3D)*3*pyramid_glyph.ntri);
  nrm= (yPoint3D *) malloc(sizeof(yPoint3D)*3*pyramid_glyph.ntri);
  pyramid_glyph.xyz= xyz;
  pyramid_glyph.nrm= nrm;
  nv= 0;
  /* first do the base (the normal is the same for all 6 vertices) */
  nrm[nv].x= 0.0;
  nrm[nv].y= 0.0;
  nrm[nv].z= -1.0;
  xyz[nv].x= -0.5;
  xyz[nv].y= -0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.5;
  xyz[nv].y= -0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.5;
  xyz[nv].y= 0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= -0.5;
  xyz[nv].y= -0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.5;
  xyz[nv].y= 0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= -0.5;
  xyz[nv].y= 0.5;
  xyz[nv].z= -0.5;
  nv++;

  /* do the four sides. all 3 vertices in a triangle
     have the same normal. */
  nrz= (float) (1.0/sqrt(5.0));
  nrt= (float) (2.0/sqrt(5.0));
  nrm[nv].x= 0.0;
  nrm[nv].y= -nrt;
  nrm[nv].z= nrz;
  xyz[nv].x= -0.5;
  xyz[nv].y= -0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.5;
  xyz[nv].y= -0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.0;
  xyz[nv].y= 0.0;
  xyz[nv].z= 0.5;
  nv++;

  nrm[nv].x= nrt;
  nrm[nv].y= 0.0;
  nrm[nv].z= nrz;
  xyz[nv].x= 0.5;
  xyz[nv].y= -0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.5;
  xyz[nv].y= 0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.0;
  xyz[nv].y= 0.0;
  xyz[nv].z= 0.5;
  nv++;

  nrm[nv].x= 0.0;
  nrm[nv].y= nrt;
  nrm[nv].z= nrz;
  xyz[nv].x= 0.5;
  xyz[nv].y= 0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= -0.5;
  xyz[nv].y= 0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.0;
  xyz[nv].y= 0.0;
  xyz[nv].z= 0.5;
  nv++;

  nrm[nv].x= -nrt;
  nrm[nv].y= 0.0;
  nrm[nv].z= nrz;
  xyz[nv].x= -0.5;
  xyz[nv].y= 0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= -0.5;
  xyz[nv].y= -0.5;
  xyz[nv].z= -0.5;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= 0.0;
  xyz[nv].y= 0.0;
  xyz[nv].z= 0.5;
  nv++;
}

void makTetGlyph(void)
{
  int nv;
  yPoint3D *xyz, *nrm;
  float val;

  /* The tetrahedron is made up of 4 triangles.
     Compute and save the coordinates of their vertices. */
  tet_glyph.ntri= 4;
  xyz= (yPoint3D *) malloc(sizeof(yPoint3D)*3*tet_glyph.ntri);
  nrm= (yPoint3D *) malloc(sizeof(yPoint3D)*3*tet_glyph.ntri);
  tet_glyph.xyz= xyz;
  tet_glyph.nrm= nrm;
  nv= 0;
  val= (float) (0.5*sqrt(3.0));
  /* first do the base (the normal is the same for all 3 vertices) */
  nrm[nv].x= 0.0;
  nrm[nv].y= 0.0;
  nrm[nv].z= -1.0;
  xyz[nv].x= 1.0;
  xyz[nv].y= 0.0;
  xyz[nv].z= 0.0;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= -0.5;
  xyz[nv].y= val;
  xyz[nv].z= 0.0;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv].x= -0.5;
  xyz[nv].y= -val;
  xyz[nv].z= 0.0;
  nv++;

  /* Do the three sides. All 3 vertices in a triangle
     have the same normal. */
  nrm[nv].x= HT_RAT;
  nrm[nv].y= sqrt(3.0)*HT_RAT;
  nrm[nv].z= 1.0;
  val= (float) (1.0/sqrt(nrm[nv].x*nrm[nv].x+nrm[nv].y*nrm[nv].y+nrm[nv].z*nrm[nv].z));
  nrm[nv].x *= val; nrm[nv].y *= val; nrm[nv].z *= val;
  xyz[nv].x= 0.0;
  xyz[nv].y= 0.0;
  xyz[nv].z= HT_RAT;
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv]= xyz[1];
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv]= xyz[0];
  nv++;

  nrm[nv].x= -2.0*HT_RAT;
  nrm[nv].y= 0.0;
  nrm[nv].z= 1.0;
  val= (float) (1.0/sqrt(nrm[nv].x*nrm[nv].x+nrm[nv].y*nrm[nv].y+nrm[nv].z*nrm[nv].z));
  nrm[nv].x *= val; nrm[nv].y *= val; nrm[nv].z *= val;
  xyz[nv]= xyz[1];
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv]= xyz[3];
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv]= xyz[2];
  nv++;

  nrm[nv].x= HT_RAT;
  nrm[nv].y= -sqrt(3.0)*HT_RAT;
  nrm[nv].z= 1.0;
  val= (float) (1.0/sqrt(nrm[nv].x*nrm[nv].x+nrm[nv].y*nrm[nv].y+nrm[nv].z*nrm[nv].z));
  nrm[nv].x *= val; nrm[nv].y *= val; nrm[nv].z *= val;
  xyz[nv]= xyz[2];
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv]= xyz[3];
  nv++;
  nrm[nv]= nrm[nv-1];
  xyz[nv]= xyz[0];
  nv++;
}

void yglGlyphs_old(long nglyph, float *origins, float *hite, float *base,
                     float *theta, float *phi, float *colr)
{
  int ng;
  float x, y, z, fac, x0, y0, z0, ht, bs, b;

  /* The input is a list of glyphs (ellipsoids). 
     Each (pyramidal) glyph has a center, a height, a base, a color,
     and a direction (theta and phi).
  */

  /* draw the glyph list now */
  if(nglyph <= 0) return;

  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();
  CHEK_ERROR("yglGlyphs setup");
  for(ng= 0; ng < nglyph; ng++) {
    glColor3fv(colr+3*ng);
    x0= origins[3*ng];
    y0= origins[3*ng+1];
    z0= origins[3*ng+2];
    ht= hite[ng];
    bs= base[ng];
    b= 0.5f*bs;
    /* Draw a solid pyramid with specified height 
       and base centered at the specified origin with
       the specified orientation.
    */
    x= x0-b;
    y= y0-b;
    z= z0-0.5f*ht;
    fac= (float) (1.0/sqrt(b*b+ht*ht));
    glBegin(GL_QUADS);
    glNormal3f(0.0f, 0.0f, 1.0f);
    glVertex3f(x, y, z);
    glVertex3f(x+bs, y, z);
    glVertex3f(x+bs, y+bs, z);
    glVertex3f(x, y+bs, z);
    glEnd();
    glBegin(GL_TRIANGLES);
    glNormal3f(0.0f, -ht*fac, b*fac);
    glVertex3f(x,    y, z);
    glVertex3f(x+bs, y, z);
    glVertex3f(x0,   y0, z+ht);
    glNormal3f(ht*fac, 0.0f, b*fac);
    glVertex3f(x+bs, y,    z);
    glVertex3f(x+bs, y+bs, z);
    glVertex3f(x0,   y0,   z+ht);
    glNormal3f(0.0f, ht*fac, b*fac);
    glVertex3f(x+bs, y+bs, z);
    glVertex3f(x,    y+bs, z);
    glVertex3f(x0,   y0,   z+ht);
    glNormal3f(-ht*fac, 0.0f, b*fac);
    glVertex3f(x, y+bs, z);
    glVertex3f(x, y,    z);
    glVertex3f(x0, y0,  z+ht);
    glEnd();
  }
  CHEK_ERROR("yglGlyphs");
}

void yglEllipsoids(long nglyph, float *origins, float *radii, float *ellip,
                     float *theta, float *phi, float *colr)
{
  int i, j, ng, numi, numj;
  double lat, lon, x, y, z, csth, snth, csph, snph, pi, dth;
  double csth1, snth1, fac, x0, y0, z0, eps, radius;
#define ang_step_e 15

  /* The input is a list of glyphs (ellipsoids). 
     Each ellipsoid has a center, a radius, an ellipticity, a color,
     and a direction (theta and phi).
  */

  /* draw the glyph list now */
  if(nglyph <= 0) return;

  pi= 4.0*atan(1.0);
  numi= (int)(360.0/ang_step_e)+1;
  numj= (int)(180.0/ang_step_e)+1;
  dth= pi/numj;
  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();
  CHEK_ERROR("yglEllipsoids setup");
  for(ng= 0; ng < nglyph; ng++) {
    glColor3fv(colr+3*ng);
	x0= origins[3*ng];
	y0= origins[3*ng+1];
	z0= origins[3*ng+2];
	eps= ellip[ng];
	radius= radii[ng];
    /* Draw a solid ellipsoid with specified radius 
       and ellipticity centered at the specified origin. 
       Points are every 15 degrees in latitude and longitude. 
    */
    /* do all latitudes */
    for(j= 0; j < numj; j++) {
      lat= j*pi/numj;
      csth= cos(lat);
      snth= sqrt(1.0-csth*csth);
      csth1= cos(lat+dth);
      snth1= sqrt(1.0-csth1*csth1);
      glBegin(GL_TRIANGLE_STRIP);
      /* do all longitudes */
      for(i= 0; i <= numi; i++) {
        lon= i*2.0*pi/numi;
        csph= cos(lon);
        snph= sin(lon);
        x= csph*snth;
        y= snph*snth;
        z= eps*csth;
	    fac= 1.0/sqrt(x*x+y*x+z*z);
        glNormal3f((float)(x*fac),(float)(y*fac),(float)(z*fac));
        x= x0+eps*radius*csph*snth;
        y= y0+eps*radius*snph*snth;
        z= z0+radius*csth;
        glVertex3f((float)x,(float)y,(float)z);
        x= csph*snth1;
        y= snph*snth1;
        z= eps*csth1;
	    fac= 1.0/sqrt(x*x+y*x+z*z);
        glNormal3f((float)(x*fac),(float)(y*fac),(float)(z*fac));
        x= x0+eps*radius*csph*snth1;
        y= y0+eps*radius*snph*snth1;
        z= z0+radius*csth1;
        glVertex3f((float)x,(float)y,(float)z);
      }
      glEnd();
    }
  }
  CHEK_ERROR("yglEllipsoids");
}
