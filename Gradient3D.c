/*
 * $Id: Gradient3D.c,v 1.1 2005-09-18 22:07:58 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include "Contour3D.h"
#include <stdlib.h>
#include <math.h>

void ycNormalize(yPoint3D *n)
{
  /* normalize the input vector. May fail for points larger than
     the square root of the maximum double. Answer may be wrong 
	 for point smaller than 1.0e-80; */

  double len;
  len= sqrt(n->x*n->x+n->y*n->y+n->z*n->z)+1.0e-80;
  n->x /= len;
  n->y /= len;
  n->z /= len;
}

void ycPointGradient(long i, long j, long k, long nx, long ny, long nz, double *s, 
                   double dx, double dy, double dz, yPoint3D *n)
{
  long sliceSize= nx*ny;
  long ibase= i + j*nx + k*sliceSize;

  /* x-direction */
  if ( i == 0 ) {
    n->x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( i == (nx-1) ) {
    n->x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    n->x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }

  /* y-direction */
  if ( j == 0 ) {
    n->y = (s[ibase + nx] - s[ibase]) / dy;
  } else if ( j == (ny-1) ) {
    n->y = (s[ibase] - s[ibase - nx]) / dy;
  } else {
    n->y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  }

  /* z-direction */
  if ( k == 0 ) {
    n->z = (s[ibase + sliceSize] - s[ibase]) / dz;
  } else if ( k == (nz-1) ) {
    n->z = (s[ibase] - s[ibase - sliceSize]) / dz;
  } else {
    n->z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  }
}

void ycPointGradientAll(long i, long j, long k, long nx, long ny, long nz, double *s, 
                   double dx, double dy, double dz, yPoint3D gradient[8])
{
  long ii, jj, kk;
  long sliceSize= nx*ny;
  long ibase, vrt= 0;

  /* guaranteed that in j and k directions the point is in the
     interior of the mesh. Have to test edge cases for i direction */
  ii= i; jj= j; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  /* x-direction */
  if ( ii == 0 ) {
    gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( ii == (nx-1) ) {
    gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }
  /* y-direction */
  gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  /* z-direction */
  gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  vrt++;

  ii= i+1; jj= j; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  /* x-direction */
  if ( ii == 0 ) {
    gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( ii == (nx-1) ) {
    gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }
  /* y-direction */
  gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  /* z-direction */
  gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  vrt++;

  ii= i+1; jj= j+1; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  /* x-direction */
  if ( ii == 0 ) {
    gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( ii == (nx-1) ) {
    gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }
  /* y-direction */
  gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  /* z-direction */
  gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  vrt++;

  ii= i; jj= j+1; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  /* x-direction */
  if ( ii == 0 ) {
    gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( ii == (nx-1) ) {
    gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }
  /* y-direction */
  gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  /* z-direction */
  gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  vrt++;

  ii= i; jj= j; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  /* x-direction */
  if ( ii == 0 ) {
    gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( ii == (nx-1) ) {
    gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }
  /* y-direction */
  gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  /* z-direction */
  gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  vrt++;

  ii= i+1; jj= j; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  /* x-direction */
  if ( ii == 0 ) {
    gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( ii == (nx-1) ) {
    gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }
  /* y-direction */
  gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  /* z-direction */
  gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  vrt++;

  ii= i+1; jj= j+1; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  /* x-direction */
  if ( ii == 0 ) {
    gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( ii == (nx-1) ) {
    gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }
  /* y-direction */
  gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  /* z-direction */
  gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  vrt++;

  ii= i; jj= j+1; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  /* x-direction */
  if ( ii == 0 ) {
    gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( ii == (nx-1) ) {
    gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }
  /* y-direction */
  gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  /* z-direction */
  gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  vrt++;
}

void ycPointGradientCrv(long i, long j, long k, long nx, long ny, long nz, 
					   yPoint3D *x, double *s, yPoint3D *n)
{
  double sp, sm, del2;
  long sliceSize= nx*ny;
  long ibase= i + j*nx + k*sliceSize;
  yPoint3D delta;

  /* i-direction */
  if ( i == 0 ) {
    sp = s[ibase+1];
    sm = s[ibase];
	delta.x= x[ibase+1].x-x[ibase].x;
	delta.y= x[ibase+1].y-x[ibase].y;
	delta.z= x[ibase+1].z-x[ibase].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x = (sp - sm)*delta.x / del2;
    n->y = (sp - sm)*delta.y / del2;
    n->z = (sp - sm)*delta.z / del2;
  } else if ( i == (nx-1) ) {
    sp = s[ibase];
    sm = s[ibase-1];
	delta.x= x[ibase].x-x[ibase-1].x;
	delta.y= x[ibase].y-x[ibase-1].y;
	delta.z= x[ibase].z-x[ibase-1].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x = (sp - sm)*delta.x / del2;
    n->y = (sp - sm)*delta.y / del2;
    n->z = (sp - sm)*delta.z / del2;
  } else {
    sp = s[ibase+1];
    sm = s[ibase-1];
	delta.x= x[ibase+1].x-x[ibase-1].x;
	delta.y= x[ibase+1].y-x[ibase-1].y;
	delta.z= x[ibase+1].z-x[ibase-1].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x = (sp - sm)*delta.x / del2;
    n->y = (sp - sm)*delta.y / del2;
    n->z = (sp - sm)*delta.z / del2;
  }

  /* j-direction */
  if ( j == 0 ) {
    sp = s[ibase + nx];
    sm = s[ibase];
	delta.x= x[ibase+nx].x-x[ibase].x;
	delta.y= x[ibase+nx].y-x[ibase].y;
	delta.z= x[ibase+nx].z-x[ibase].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  } else if ( j == (ny-1) ) {
    sp = s[ibase];
    sm = s[ibase - nx];
	delta.x= x[ibase].x-x[ibase-nx].x;
	delta.y= x[ibase].y-x[ibase-nx].y;
	delta.z= x[ibase].z-x[ibase-nx].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  } else {
    sp = s[ibase + nx];
    sm = s[ibase - nx];
	delta.x= x[ibase+nx].x-x[ibase-nx].x;
	delta.y= x[ibase+nx].y-x[ibase-nx].y;
	delta.z= x[ibase+nx].z-x[ibase-nx].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  }

  /* k-direction */
  if ( k == 0 ) {
    sp = s[ibase + sliceSize];
    sm = s[ibase];
	delta.x= x[ibase+sliceSize].x-x[ibase].x;
	delta.y= x[ibase+sliceSize].y-x[ibase].y;
	delta.z= x[ibase+sliceSize].z-x[ibase].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  } else if ( k == (nz-1) ) {
    sp = s[ibase];
    sm = s[ibase - sliceSize];
	delta.x= x[ibase].x-x[ibase-sliceSize].x;
	delta.y= x[ibase].y-x[ibase-sliceSize].y;
	delta.z= x[ibase].z-x[ibase-sliceSize].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  } else {
    sp = s[ibase + sliceSize];
    sm = s[ibase - sliceSize];
	delta.x= x[ibase+sliceSize].x-x[ibase-sliceSize].x;
	delta.y= x[ibase+sliceSize].y-x[ibase-sliceSize].y;
	delta.z= x[ibase+sliceSize].z-x[ibase-sliceSize].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  }
}

#define del_2(a) {del2=a.x*a.x+a.y*a.y+a.z*a.z+1.0e-80;}
#define do_diff(a,n,m) {a.x= pts[n].x-pts[m].x;a.y= pts[n].y-pts[m].y;a.z= pts[n].z-pts[m].z;}
#define bump_grad(n,a) {gradients[n].x += sp*a.x/del2;gradients[n].y += sp*a.y/del2;gradients[n].z += sp*a.z/del2;}
#define set_grad(n,a) {gradients[n].x= sp*a.x/del2;gradients[n].y= sp*a.y/del2;gradients[n].z= sp*a.z/del2;}

void ycPointGradientHex(yPoint3D *pts, double *s, yPoint3D *gradients)
{
  double sp, del2;
  yPoint3D d100_000, d110_010, d101_001, d111_011;
  yPoint3D d010_000, d110_100, d011_001, d111_101;
  yPoint3D d001_000, d101_100, d011_010, d111_110;

  do_diff(d100_000,1,0)
  do_diff(d110_010,3,2)
  do_diff(d101_001,5,4)
  do_diff(d111_011,7,6)

  do_diff(d010_000,2,0)
  do_diff(d110_100,3,1)
  do_diff(d011_001,6,4)
  do_diff(d111_101,7,5)

  do_diff(d001_000,4,0)
  do_diff(d101_100,5,1)
  do_diff(d011_010,6,2)
  do_diff(d111_110,7,3)

  /* x=0, y=0, z=0 corner */
  sp = s[1]-s[0];
  del_2(d100_000);
  set_grad(0,d100_000);
  sp = s[2]-s[0];
  del_2(d010_000);
  bump_grad(0,d010_000);
  sp = s[4]-s[0];
  del_2(d001_000);
  bump_grad(0,d001_000);

  /* x=1, y=0, z=0 corner */
  sp = s[1]-s[0];
  del_2(d100_000);
  set_grad(1,d100_000);
  sp = s[3]-s[1];
  del_2(d110_100);
  bump_grad(1,d110_100);
  sp = s[5]-s[1];
  del_2(d101_100);
  bump_grad(1,d101_100);

  /* x=0, y=1, z=0 corner */
  sp = s[3]-s[2];
  del_2(d110_010);
  set_grad(3,d110_010);
  sp = s[2]-s[0];
  del_2(d010_000);
  bump_grad(3,d010_000);
  sp = s[6]-s[2];
  del_2(d011_010);
  bump_grad(3,d011_010);

  /* x=1, y=1, z=0 corner */
  sp = s[3]-s[2];
  del_2(d110_010);
  set_grad(2,d110_010);
  sp = s[3]-s[1];
  del_2(d110_100);
  bump_grad(2,d110_100);
  sp = s[7]-s[3];
  del_2(d111_110);
  bump_grad(2,d111_110);

  /* x=0, y=0, z=1 corner */
  sp = s[5]-s[4];
  del_2(d101_001);
  set_grad(4,d101_001);
  sp = s[6]-s[4];
  del_2(d011_001);
  bump_grad(4,d011_001);
  sp = s[4]-s[0];
  del_2(d001_000);
  bump_grad(4,d001_000);

  /* x=1, y=0, z=1 corner */
  sp = s[5]-s[4];
  del_2(d101_001);
  set_grad(5,d101_001);
  sp = s[7]-s[5];
  del_2(d111_101);
  bump_grad(5,d111_101);
  sp = s[5]-s[1];
  del_2(d101_100);
  bump_grad(5,d101_100);

  /* x=0, y=1, z=1 corner */
  sp = s[7]-s[6];
  del_2(d111_011);
  set_grad(7,d111_011);
  sp = s[6]-s[4];
  del_2(d011_001);
  bump_grad(7,d011_001);
  sp = s[6]-s[2];
  del_2(d011_010);
  bump_grad(7,d011_010);

  /* x=1, y=1, z=1 corner */
  sp = s[7]-s[6];
  del_2(d111_011);
  set_grad(6,d111_011);
  sp = s[7]-s[5];
  del_2(d111_101);
  bump_grad(6,d111_101);
  sp = s[7]-s[3];
  del_2(d111_110);
  bump_grad(6,d111_110);
}

void ycPointGradientGrd(long i, long j, long k, long nx, long ny, long nz, 
						double *s, double dx, double dy, double dz, yPoint3D *n, 
						yPoint3D *grd, char *done)
{
  long sliceSize= nx*ny;
  long ibase= i + j*nx + k*sliceSize;

  if(done[ibase]) {
	*n = grd[ibase];
	return;
  }
  /* x-direction */
  if ( i == 0 ) {
    n->x = (s[ibase+1] - s[ibase]) / dx;
  } else if ( i == (nx-1) ) {
    n->x = (s[ibase] - s[ibase-1]) / dx;
  } else {
    n->x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
  }

  /* y-direction */
  if ( j == 0 ) {
    n->y = (s[ibase + nx] - s[ibase]) / dy;
  } else if ( j == (ny-1) ) {
    n->y = (s[ibase] - s[ibase - nx]) / dy;
  } else {
    n->y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
  }

  /* z-direction */
  if ( k == 0 ) {
    n->z = (s[ibase + sliceSize] - s[ibase]) / dz;
  } else if ( k == (nz-1) ) {
    n->z = (s[ibase] - s[ibase - sliceSize]) / dz;
  } else {
    n->z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
  }
  grd[ibase] = *n;
  done[ibase]= 1;
}

void ycPointGradientAllGrd(long i, long j, long k, long nx, long ny, long nz, double *s, 
                   double dx, double dy, double dz, yPoint3D gradient[8], 
				   yPoint3D *grd, char *done)
{
  long ii, jj, kk;
  long sliceSize= nx*ny;
  long ibase, vrt= 0;

  /* guaranteed that in j and k directions the point is in the
     interior of the mesh. Have to test edge cases for i direction */
  ii= i; jj= j; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    if ( ii == 0 ) {
      gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
	} else if ( ii == (nx-1) ) {
      gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
	} else {
      gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
	}
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    if ( ii == 0 ) {
      gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
	} else if ( ii == (nx-1) ) {
      gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
	} else {
      gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
	}
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j+1; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    if ( ii == 0 ) {
      gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
	} else if ( ii == (nx-1) ) {
      gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
	} else {
      gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
	}
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j+1; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    if ( ii == 0 ) {
      gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
	} else if ( ii == (nx-1) ) {
      gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
	} else {
      gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
	}
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    if ( ii == 0 ) {
      gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
	} else if ( ii == (nx-1) ) {
      gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
	} else {
      gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
	}
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    if ( ii == 0 ) {
      gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
	} else if ( ii == (nx-1) ) {
      gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
	} else {
      gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
	}
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j+1; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    if ( ii == 0 ) {
      gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
	} else if ( ii == (nx-1) ) {
      gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
	} else {
      gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
	}
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j+1; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    if ( ii == 0 ) {
      gradient[vrt].x = (s[ibase+1] - s[ibase]) / dx;
	} else if ( ii == (nx-1) ) {
      gradient[vrt].x = (s[ibase] - s[ibase-1]) / dx;
	} else {
      gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
	}
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;
}

void ycPointGradientCrvGrd(long i, long j, long k, long nx, long ny, long nz, 
					   yPoint3D *x, double *s, yPoint3D *n, yPoint3D *grd,
					   char *done)
{
  double sp, sm, del2;
  long sliceSize= nx*ny;
  long ibase= i + j*nx + k*sliceSize;
  yPoint3D delta;

  if(done[ibase]) {
	*n = grd[ibase];
	return;
  }
  /* i-direction */
  if ( i == 0 ) {
    sp = s[ibase+1];
    sm = s[ibase];
	delta.x= x[ibase+1].x-x[ibase].x;
	delta.y= x[ibase+1].y-x[ibase].y;
	delta.z= x[ibase+1].z-x[ibase].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x = (sp - sm)*delta.x / del2;
    n->y = (sp - sm)*delta.y / del2;
    n->z = (sp - sm)*delta.z / del2;
  } else if ( i == (nx-1) ) {
    sp = s[ibase];
    sm = s[ibase-1];
	delta.x= x[ibase].x-x[ibase-1].x;
	delta.y= x[ibase].y-x[ibase-1].y;
	delta.z= x[ibase].z-x[ibase-1].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x = (sp - sm)*delta.x / del2;
    n->y = (sp - sm)*delta.y / del2;
    n->z = (sp - sm)*delta.z / del2;
  } else {
    sp = s[ibase+1];
    sm = s[ibase-1];
	delta.x= x[ibase+1].x-x[ibase-1].x;
	delta.y= x[ibase+1].y-x[ibase-1].y;
	delta.z= x[ibase+1].z-x[ibase-1].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x = (sp - sm)*delta.x / del2;
    n->y = (sp - sm)*delta.y / del2;
    n->z = (sp - sm)*delta.z / del2;
  }

  /* j-direction */
  if ( j == 0 ) {
    sp = s[ibase + nx];
    sm = s[ibase];
	delta.x= x[ibase+nx].x-x[ibase].x;
	delta.y= x[ibase+nx].y-x[ibase].y;
	delta.z= x[ibase+nx].z-x[ibase].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  } else if ( j == (ny-1) ) {
    sp = s[ibase];
    sm = s[ibase - nx];
	delta.x= x[ibase].x-x[ibase-nx].x;
	delta.y= x[ibase].y-x[ibase-nx].y;
	delta.z= x[ibase].z-x[ibase-nx].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  } else {
    sp = s[ibase + nx];
    sm = s[ibase - nx];
	delta.x= x[ibase+nx].x-x[ibase-nx].x;
	delta.y= x[ibase+nx].y-x[ibase-nx].y;
	delta.z= x[ibase+nx].z-x[ibase-nx].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  }

  /* k-direction */
  if ( k == 0 ) {
    sp = s[ibase + sliceSize];
    sm = s[ibase];
	delta.x= x[ibase+sliceSize].x-x[ibase].x;
	delta.y= x[ibase+sliceSize].y-x[ibase].y;
	delta.z= x[ibase+sliceSize].z-x[ibase].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  } else if ( k == (nz-1) ) {
    sp = s[ibase];
    sm = s[ibase - sliceSize];
	delta.x= x[ibase].x-x[ibase-sliceSize].x;
	delta.y= x[ibase].y-x[ibase-sliceSize].y;
	delta.z= x[ibase].z-x[ibase-sliceSize].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  } else {
    sp = s[ibase + sliceSize];
    sm = s[ibase - sliceSize];
	delta.x= x[ibase+sliceSize].x-x[ibase-sliceSize].x;
	delta.y= x[ibase+sliceSize].y-x[ibase-sliceSize].y;
	delta.z= x[ibase+sliceSize].z-x[ibase-sliceSize].z;
	del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    n->x += (sp - sm) *delta.x / del2;
    n->y += (sp - sm) *delta.y / del2;
    n->z += (sp - sm) *delta.z / del2;
  }
  grd[ibase]= *n;
  done[ibase]= 1;
}

void ycPointGradientIntGrd(long i, long j, long k, long nx, long ny, long nz, double *s, 
                   double dx, double dy, double dz, yPoint3D gradient[8], 
				   yPoint3D *grd, char *done)
{
  long ii, jj, kk;
  long sliceSize= nx*ny;
  long ibase, vrt= 0;

  /* guaranteed that in i, j and k directions the point is in the
     interior of the mesh. */
  ii= i; jj= j; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j+1; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j+1; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j+1; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j+1; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* x-direction */
    gradient[vrt].x = 0.5 * (s[ibase+1] - s[ibase-1]) / dx;
    /* y-direction */
    gradient[vrt].y = 0.5 * (s[ibase + nx] - s[ibase - nx]) / dy;
    /* z-direction */
    gradient[vrt].z = 0.5 * (s[ibase + sliceSize] - s[ibase - sliceSize]) / dz;
	grd[ibase]= gradient[vrt];
	done[ibase]= 1;
  }
  vrt++;
}

void ycPointGradientIntGrdCrv(long i, long j, long k, long nx, long ny, long nz, 
				   yPoint3D *xyz, double *s, yPoint3D gradient[8], 
				   yPoint3D *grd, char *done)
{
  long ii, jj, kk;
  long sliceSize= nx*ny;
  long ibase, vrt= 0;
  double fac, del2;
  yPoint3D delta;

  /* guaranteed that in i, j and k directions the point is in the
     interior of the mesh. */
  ii= i; jj= j; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* i-direction */
    delta.x= xyz[ibase+1].x-xyz[ibase-1].x;
    delta.y= xyz[ibase+1].y-xyz[ibase-1].y;
    delta.z= xyz[ibase+1].z-xyz[ibase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+1] - s[ibase-1]) / del2;
    gradient[vrt].x = delta.x*fac;
    gradient[vrt].y = delta.y*fac;
    gradient[vrt].z = delta.z*fac;

    /* j-direction */
    delta.x= xyz[ibase+nx].x-xyz[ibase-nx].x;
    delta.y= xyz[ibase+nx].y-xyz[ibase-nx].y;
    delta.z= xyz[ibase+nx].z-xyz[ibase-nx].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+nx] - s[ibase-nx]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;

    /* k-direction */
    delta.x= xyz[ibase+sliceSize].x-xyz[ibase-sliceSize].x;
    delta.y= xyz[ibase+sliceSize].y-xyz[ibase-sliceSize].y;
    delta.z= xyz[ibase+sliceSize].z-xyz[ibase-sliceSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+sliceSize] - s[ibase-sliceSize]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;
    grd[ibase]= gradient[vrt];
    done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* i-direction */
    delta.x= xyz[ibase+1].x-xyz[ibase-1].x;
    delta.y= xyz[ibase+1].y-xyz[ibase-1].y;
    delta.z= xyz[ibase+1].z-xyz[ibase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+1] - s[ibase-1]) / del2;
    gradient[vrt].x = delta.x*fac;
    gradient[vrt].y = delta.y*fac;
    gradient[vrt].z = delta.z*fac;

    /* j-direction */
    delta.x= xyz[ibase+nx].x-xyz[ibase-nx].x;
    delta.y= xyz[ibase+nx].y-xyz[ibase-nx].y;
    delta.z= xyz[ibase+nx].z-xyz[ibase-nx].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+nx] - s[ibase-nx]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;

    /* k-direction */
    delta.x= xyz[ibase+sliceSize].x-xyz[ibase-sliceSize].x;
    delta.y= xyz[ibase+sliceSize].y-xyz[ibase-sliceSize].y;
    delta.z= xyz[ibase+sliceSize].z-xyz[ibase-sliceSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+sliceSize] - s[ibase-sliceSize]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;
    grd[ibase]= gradient[vrt];
    done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j+1; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* i-direction */
    delta.x= xyz[ibase+1].x-xyz[ibase-1].x;
    delta.y= xyz[ibase+1].y-xyz[ibase-1].y;
    delta.z= xyz[ibase+1].z-xyz[ibase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+1] - s[ibase-1]) / del2;
    gradient[vrt].x = delta.x*fac;
    gradient[vrt].y = delta.y*fac;
    gradient[vrt].z = delta.z*fac;

    /* j-direction */
    delta.x= xyz[ibase+nx].x-xyz[ibase-nx].x;
    delta.y= xyz[ibase+nx].y-xyz[ibase-nx].y;
    delta.z= xyz[ibase+nx].z-xyz[ibase-nx].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+nx] - s[ibase-nx]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;

    /* k-direction */
    delta.x= xyz[ibase+sliceSize].x-xyz[ibase-sliceSize].x;
    delta.y= xyz[ibase+sliceSize].y-xyz[ibase-sliceSize].y;
    delta.z= xyz[ibase+sliceSize].z-xyz[ibase-sliceSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+sliceSize] - s[ibase-sliceSize]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;
    grd[ibase]= gradient[vrt];
    done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j+1; kk= k;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* i-direction */
    delta.x= xyz[ibase+1].x-xyz[ibase-1].x;
    delta.y= xyz[ibase+1].y-xyz[ibase-1].y;
    delta.z= xyz[ibase+1].z-xyz[ibase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+1] - s[ibase-1]) / del2;
    gradient[vrt].x = delta.x*fac;
    gradient[vrt].y = delta.y*fac;
    gradient[vrt].z = delta.z*fac;

    /* j-direction */
    delta.x= xyz[ibase+nx].x-xyz[ibase-nx].x;
    delta.y= xyz[ibase+nx].y-xyz[ibase-nx].y;
    delta.z= xyz[ibase+nx].z-xyz[ibase-nx].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+nx] - s[ibase-nx]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;

    /* k-direction */
    delta.x= xyz[ibase+sliceSize].x-xyz[ibase-sliceSize].x;
    delta.y= xyz[ibase+sliceSize].y-xyz[ibase-sliceSize].y;
    delta.z= xyz[ibase+sliceSize].z-xyz[ibase-sliceSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+sliceSize] - s[ibase-sliceSize]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;
    grd[ibase]= gradient[vrt];
    done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* i-direction */
    delta.x= xyz[ibase+1].x-xyz[ibase-1].x;
    delta.y= xyz[ibase+1].y-xyz[ibase-1].y;
    delta.z= xyz[ibase+1].z-xyz[ibase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+1] - s[ibase-1]) / del2;
    gradient[vrt].x = delta.x*fac;
    gradient[vrt].y = delta.y*fac;
    gradient[vrt].z = delta.z*fac;

    /* j-direction */
    delta.x= xyz[ibase+nx].x-xyz[ibase-nx].x;
    delta.y= xyz[ibase+nx].y-xyz[ibase-nx].y;
    delta.z= xyz[ibase+nx].z-xyz[ibase-nx].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+nx] - s[ibase-nx]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;

    /* k-direction */
    delta.x= xyz[ibase+sliceSize].x-xyz[ibase-sliceSize].x;
    delta.y= xyz[ibase+sliceSize].y-xyz[ibase-sliceSize].y;
    delta.z= xyz[ibase+sliceSize].z-xyz[ibase-sliceSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+sliceSize] - s[ibase-sliceSize]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;
    grd[ibase]= gradient[vrt];
    done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* i-direction */
    delta.x= xyz[ibase+1].x-xyz[ibase-1].x;
    delta.y= xyz[ibase+1].y-xyz[ibase-1].y;
    delta.z= xyz[ibase+1].z-xyz[ibase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+1] - s[ibase-1]) / del2;
    gradient[vrt].x = delta.x*fac;
    gradient[vrt].y = delta.y*fac;
    gradient[vrt].z = delta.z*fac;

    /* j-direction */
    delta.x= xyz[ibase+nx].x-xyz[ibase-nx].x;
    delta.y= xyz[ibase+nx].y-xyz[ibase-nx].y;
    delta.z= xyz[ibase+nx].z-xyz[ibase-nx].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+nx] - s[ibase-nx]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;

    /* k-direction */
    delta.x= xyz[ibase+sliceSize].x-xyz[ibase-sliceSize].x;
    delta.y= xyz[ibase+sliceSize].y-xyz[ibase-sliceSize].y;
    delta.z= xyz[ibase+sliceSize].z-xyz[ibase-sliceSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+sliceSize] - s[ibase-sliceSize]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;
    grd[ibase]= gradient[vrt];
    done[ibase]= 1;
  }
  vrt++;

  ii= i+1; jj= j+1; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* i-direction */
    delta.x= xyz[ibase+1].x-xyz[ibase-1].x;
    delta.y= xyz[ibase+1].y-xyz[ibase-1].y;
    delta.z= xyz[ibase+1].z-xyz[ibase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+1] - s[ibase-1]) / del2;
    gradient[vrt].x = delta.x*fac;
    gradient[vrt].y = delta.y*fac;
    gradient[vrt].z = delta.z*fac;

    /* j-direction */
    delta.x= xyz[ibase+nx].x-xyz[ibase-nx].x;
    delta.y= xyz[ibase+nx].y-xyz[ibase-nx].y;
    delta.z= xyz[ibase+nx].z-xyz[ibase-nx].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+nx] - s[ibase-nx]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;

    /* k-direction */
    delta.x= xyz[ibase+sliceSize].x-xyz[ibase-sliceSize].x;
    delta.y= xyz[ibase+sliceSize].y-xyz[ibase-sliceSize].y;
    delta.z= xyz[ibase+sliceSize].z-xyz[ibase-sliceSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+sliceSize] - s[ibase-sliceSize]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;
    grd[ibase]= gradient[vrt];
    done[ibase]= 1;
  }
  vrt++;

  ii= i; jj= j+1; kk= k+1;
  ibase= ii + jj*nx + kk*sliceSize;
  if(done[ibase]) {
	gradient[vrt] = grd[ibase];
  } else {
    /* i-direction */
    delta.x= xyz[ibase+1].x-xyz[ibase-1].x;
    delta.y= xyz[ibase+1].y-xyz[ibase-1].y;
    delta.z= xyz[ibase+1].z-xyz[ibase-1].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+1] - s[ibase-1]) / del2;
    gradient[vrt].x = delta.x*fac;
    gradient[vrt].y = delta.y*fac;
    gradient[vrt].z = delta.z*fac;

    /* j-direction */
    delta.x= xyz[ibase+nx].x-xyz[ibase-nx].x;
    delta.y= xyz[ibase+nx].y-xyz[ibase-nx].y;
    delta.z= xyz[ibase+nx].z-xyz[ibase-nx].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+nx] - s[ibase-nx]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;

    /* k-direction */
    delta.x= xyz[ibase+sliceSize].x-xyz[ibase-sliceSize].x;
    delta.y= xyz[ibase+sliceSize].y-xyz[ibase-sliceSize].y;
    delta.z= xyz[ibase+sliceSize].z-xyz[ibase-sliceSize].z;
    del2= delta.x*delta.x+delta.y*delta.y+delta.z*delta.z+1.0e-80;
    fac= (s[ibase+sliceSize] - s[ibase-sliceSize]) / del2;
    gradient[vrt].x += delta.x*fac;
    gradient[vrt].y += delta.y*fac;
    gradient[vrt].z += delta.z*fac;
    grd[ibase]= gradient[vrt];
    done[ibase]= 1;
  }
  vrt++;
}
