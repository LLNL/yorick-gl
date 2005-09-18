/*
 * $Id: testgl.i,v 1.1.1.1 2005-09-18 22:07:52 dhmunro Exp $
 * Demonstration of 3D plots using the OpenGL based
 * 3D graphics package in yorick.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

require, "dlist3d.i"
use_interleave= 0;
require, "testisotree.i"

require,"slicenew.i";
require, "pnm.i";
require, "contour.i";

palette3d, "earth.gp";
palette3d, "stern.gp";

n_pass= 1;
keepview= 1;
do_rotate= 1;
sho_scene= 1;

have_ylm= 0;

func testgl
{
  "testgl.i - exercises the OpenGL based 3D graphics package in Yorick\n";
  "testall3d,nn  runs thru all the samples with nn-by-nn-by-nn grid";
  "and reports performance numbers";
  "demo3d,nn shows all the different kinds of plots and pauses at each";
  " ";
  "mak_sphere,ntheta,nphi  makes a sphere with the specified grid.";
  "mak_ellipsoid,ntheta,nphi,eps  makes an ellipsoid with the specified grid.";
  "rsphere rotates the sphere created by mak_sphere.";
  "mak_slice,ngrid draws a combination of iso-surfaces and a";
  "  slicing plane. ngrid is the size of the grid.";
  "dslice rotates the scene made by mak_slice";
  "rlines draws and rotates polylines in 3D";
  "sho_isotree,ngrid draws a single iso-surface using cell centered data";
  "  tri arrays and ghost cells. Uses an octree";
  "sho_isohex,ngrid draws a single iso-surface using point centered data";
  "  arbitrarily connected grid";
  "sho_volviz3,ngrid draws a volume visualization of a 3D cell array";
  "sho_texiso,ngrid draws an iso-surface inside a volume visualization";
  "sho_drum,ngrid draws an oscillating drum head";
  "Other functions may be of interest - examine the file";
}

func testall3d(size)
{
  extern make_strip;

  if(is_void(size)) size= 32;
  make_strip= 0;
  // rotate each scene the same number of times
  num_rotn= 120;
  target= 2.0;  // target run time for each "movie"

  start_tim= fin_tim= array(0.0, 3);
  timer, start_tim;

  "sho_isoreg: iso-surface on regular grid";
  sho_isoreg,size;
  stdview3d;
  do_rot,num_rotn,targ_tim=target;

  "sho_iso2: transparent iso-surfaces on a regular grid";
  sho_iso2,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isozcenreg: iso-surface on regular grid with zone-centered data";
  sho_isozcenreg,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isozcenregngrd: iso-surface on regular grid with zone-centered data, no guard";
  sho_isozcenregngrd,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isoregndx: iso-surface on regular grid returning indices into point list";
  sho_isoregndx,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isozcenregndx: iso-surface on regular grid with zone-centered data returning indices into point list";
  sho_isozcenregndx,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isocrv: iso-surface on curvilinear grid";
  sho_isocrv,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isozcencrv: iso-surface on curvilinear grid with zone-centered data";
  sho_isozcencrv,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isocrvndx: iso-surface on curvilinear grid returning indices into point list";
  sho_isocrvndx,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isozcencrvndx: iso-surface on curvilinear grid returning indices into point list";
  sho_isozcencrvndx,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isohex: iso-surface on a pile of hexahedra";
  sho_isohex,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_isotree: octtree based iso-surface on a regular grid with a color per vertex";
  sho_isotree,size,cpervrt=1;
  do_rot,num_rotn,targ_tim=target;

  "sho_isotreevarr: octtree based iso-surface on a regular grid with a color per vertex and indexed triangles";
  sho_isotreevarr,size,cpervrt=1;
  do_rot,num_rotn,targ_tim=target;

  "sho_volviz3: volume visualization";
  sho_volviz3,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_texiso: volume visualization using 3D texture hardware";
  sho_texiso,size;
  do_rot,num_rotn,targ_tim=target;

  "sho_drum: oscillating drum head";
  drum_size= size;
  sho_drum,targ_tim=target;

  "mak_slice: slice plane and iso-surface";
  mak_slice,size;
  do_rot,num_rotn,targ_tim=target;

  "mak_sphere: sphere made with quad strips and a color per zone";
  mak_sphere,4*size,4*size;
  do_rot,num_rotn,targ_tim=target;

  "mak_colrsurf: sphere made with quad strips and a color at each vertex";
  mak_colrsurf,4*size,4*size;
  do_rot,num_rotn,targ_tim=target;

  "mak_ellipsoid: ellipsoid made with tri strips";
  mak_ellipsoid,4*size,4*size,0.5;
  do_rot,num_rotn,targ_tim=target;

  "rlines: polylines";
  rlines,128*size;
  do_rot,num_rotn,targ_tim=target;

  "sho_points3d: 3D points";
  sho_points3d,size*3000;
  do_rot,num_rotn,targ_tim=target;

  "sho_glyph: 3D glyphs";
  sho_glyph,size;
  do_rot,num_rotn,targ_tim=target;

  "slicecrv_movie: slicing plane on a curvilinear grid using an octtree";
  slicecrv_movie,size;

  "slice_movie: slicing plane on a regular grid using an octtree";
  slice_movie,size;

  "isotree_movie: iso-surface on a regular grid using an octtree";
  isotree_movie,size;

  timer,fin_tim;
  write,"****************************************";
  tot_tim= (fin_tim - start_tim)(3);
  write,"time for all plots is", tot_tim;
  write,"  ";
}

func demo3d(size)
{
  extern make_strip, n_pass;

  if(is_void(size)) size= 32;
  make_strip= 0;
  old_npass= n_pass;
  n_pass= 1;

  "Iso-surfaces can be computed on Cartesian and curvilinear grids";
  sho_isoreg,size;
  stdview3d;
  pause,2000;

  "Transparent iso-surfaces can be be drawn with additional effort";
  sho_iso2,size;
  pause,2000;

  "Iso-surfaces may have color variation";
  sho_isotree,size,cpervrt=1;
  pause,2000;

  "Volume visualizations can be helpful for 3D fields";
  sho_volviz3,size;
  pause,2000;

  "This is a 3D object made from triangles computed by yorick.";
  mak_sphere,4*size,4*size;
  draw3d,1;
  pause,2000;

  "This is the same sphere, but with a backdrop";
  cage_style3d,-1;
  cage_rgb3d,[1.0,1.0,0.0];
  grid_rgb3d,[0.0,0.0,1.0];
  mak_sphere,4*size,4*size;
  draw3d,1;
  pause,2000;
  cage_style3d,0;

  "This is the same sphere, but with a backdrop and user chosen location";
  cage_style3d,1;
  cage_rgb3d,[1.0,1.0,0.0];
  grid_rgb3d,[0.0,0.0,1.0];
  cage_limits3d,[-2.0,2.0,-1.0,1.0,-1.0,3.0];
  mak_sphere,4*size,4*size;
  draw3d,1;
  pause,2000;
  cage_style3d,0;

  "This drum head is a more complex surface computed by yorick.";
  drum_size= size;
  sho_drum,targ_tim=1.0;
  pause,2000;

  "Slice planes and iso-surfaces can be combined in a single image.";
  "This image also involves clipping an iso-surface against a plane.";
  mak_slice,size;
  eye= [-2.0,2.0,0.0];
  center= [0.0,0.0,0.0];
  up= [0.0,0.0,1.0];
  lookat3d,eye,center,up;
  draw3d,1;
  pause,2000;

  "Polylines in 3D might be used to show streamlines.";
  rlines,128*size;
  draw3d,1;
  pause,2000;

  "Points in 3D might come from a particle-in-cell simulation";
  sho_points3d,size*3000;
  draw3d,1;
  pause,2000;

  "Surface plots are used to show heights above a 2D grid";
  sho_surf3d,size;
  draw3d,1;
  pause,2000;

  "Volume visualizations allow all zonesss to contribute to the final image";
  sho_texiso,size;
  lookat3d,[2.0,0.0,0.0],[0.0,0.0,0.0],[0.0,1.0,0.0];
  draw3d,1;
  pause,2000;

  n_pass= old_npass;
}


stern_red= [0x00,0x15,0x2b,0x40,0x56,0x6c,0x82,0x98,0xb8,0xd5,0xeb,0xfc,0xf6,0xf0,0xea,
0xe4,0xdd,0xd7,0xd1,0xcb,0xc4,0xbe,0xb8,0xb2,0xac,0xa0,0x9a,0x94,0x8e,0x87,
0x81,0x7b,0x75,0x6f,0x68,0x62,0x5c,0x56,0x4f,0x49,0x43,0x3c,0x31,0x2b,0x25,
0x1f,0x18,0x12,0x0c,0x06,0x40,0x41,0x42,0x43,0x44,0x46,0x47,0x48,0x4a,0x4b,
0x4d,0x4e,0x4f,0x50,0x51,0x53,0x54,0x55,0x56,0x57,0x59,0x5a,0x5b,0x5c,0x5d,
0x60,0x61,0x62,0x63,0x64,0x66,0x67,0x68,0x69,0x6a,0x6c,0x6d,0x6e,0x6f,0x70,
0x72,0x73,0x75,0x76,0x77,0x79,0x7a,0x7b,0x7c,0x7d,0x80,0x81,0x82,0x83,0x84,
0x86,0x87,0x88,0x8a,0x8b,0x8d,0x8e,0x8f,0x90,0x91,0x93,0x94,0x95,0x96,0x97,
0x99,0x9a,0x9b,0x9c,0x9d,0xa0,0xa1,0xa2,0xa3,0xa4,0xa6,0xa7,0xa8,0xa9,0xaa,
0xac,0xad,0xae,0xaf,0xb0,0xb2,0xb3,0xb5,0xb6,0xb7,0xb9,0xba,0xbb,0xbc,0xbd,
0xc0,0xc1,0xc2,0xc3,0xc4,0xc6,0xc7,0xc8,0xca,0xcb,0xcd,0xce,0xcf,0xd0,0xd1,
0xd3,0xd4,0xd5,0xd6,0xd7,0xd9,0xda,0xdb,0xdc,0xdd,0xe0,0xe1,0xe2,0xe3,0xe4,
0xe6,0xe7,0xe8,0xe9,0xea,0xec,0xed,0xee,0xef,0xf0,0xf2,0xf3,0xf5,0xf6,0xf7,
0xf9,0xfa,0xfb,0xfc,0xfd];

stern_green= [0x00,0x01,0x02,0x03,0x04,0x06,0x07,0x08,0x0a,0x0b,0x0d,0x0e,0x0f,0x10,0x11,
0x13,0x14,0x15,0x16,0x17,0x19,0x1a,0x1b,0x1c,0x1d,0x20,0x21,0x22,0x23,0x24,
0x26,0x27,0x28,0x29,0x2a,0x2c,0x2d,0x2e,0x2f,0x30,0x32,0x33,0x35,0x36,0x37,
0x39,0x3a,0x3b,0x3c,0x3d,0x40,0x41,0x42,0x43,0x44,0x46,0x47,0x48,0x4a,0x4b,
0x4d,0x4e,0x4f,0x50,0x51,0x53,0x54,0x55,0x56,0x57,0x59,0x5a,0x5b,0x5c,0x5d,
0x60,0x61,0x62,0x63,0x64,0x66,0x67,0x68,0x69,0x6a,0x6c,0x6d,0x6e,0x6f,0x70,
0x72,0x73,0x75,0x76,0x77,0x79,0x7a,0x7b,0x7c,0x7d,0x80,0x81,0x82,0x83,0x84,
0x86,0x87,0x88,0x8a,0x8b,0x8d,0x8e,0x8f,0x90,0x91,0x93,0x94,0x95,0x96,0x97,
0x99,0x9a,0x9b,0x9c,0x9d,0xa0,0xa1,0xa2,0xa3,0xa4,0xa6,0xa7,0xa8,0xa9,0xaa,
0xac,0xad,0xae,0xaf,0xb0,0xb2,0xb3,0xb5,0xb6,0xb7,0xb9,0xba,0xbb,0xbc,0xbd,
0xc0,0xc1,0xc2,0xc3,0xc4,0xc6,0xc7,0xc8,0xca,0xcb,0xcd,0xce,0xcf,0xd0,0xd1,
0xd3,0xd4,0xd5,0xd6,0xd7,0xd9,0xda,0xdb,0xdc,0xdd,0xe0,0xe1,0xe2,0xe3,0xe4,
0xe6,0xe7,0xe8,0xe9,0xea,0xec,0xed,0xee,0xef,0xf0,0xf2,0xf3,0xf5,0xf6,0xf7,
0xf9,0xfa,0xfb,0xfc,0xfd];

stern_blue= [0x00,0x01,0x03,0x06,0x08,0x0b,0x0d,0x0f,0x13,0x16,0x19,0x1b,0x1d,0x20,0x22,
0x25,0x27,0x29,0x2c,0x2e,0x31,0x33,0x35,0x38,0x3a,0x3f,0x41,0x43,0x46,0x48,
0x4b,0x4d,0x4f,0x52,0x54,0x57,0x59,0x5b,0x5e,0x60,0x63,0x65,0x69,0x6c,0x6e,
0x71,0x73,0x75,0x78,0x7a,0x7f,0x81,0x83,0x86,0x88,0x8b,0x8d,0x8f,0x93,0x96,
0x99,0x9b,0x9d,0xa0,0xa2,0xa5,0xa7,0xa9,0xac,0xae,0xb1,0xb3,0xb5,0xb8,0xba,
0xbf,0xc1,0xc3,0xc6,0xc8,0xcb,0xcd,0xcf,0xd2,0xd4,0xd7,0xd9,0xdb,0xde,0xe0,
0xe3,0xe6,0xe9,0xec,0xee,0xf1,0xf3,0xf5,0xf8,0xfa,0xfe,0xf9,0xf4,0xef,0xea,
0xe5,0xe0,0xda,0xd3,0xcc,0xc7,0xc2,0xbd,0xb8,0xb3,0xae,0xa8,0xa3,0x9f,0x99,
0x94,0x8f,0x8a,0x85,0x80,0x76,0x71,0x6c,0x67,0x62,0x5d,0x58,0x52,0x4e,0x49,
0x43,0x3e,0x39,0x34,0x2f,0x2a,0x23,0x1b,0x17,0x11,0x0c,0x07,0x02,0x02,0x06,
0x0f,0x14,0x18,0x1d,0x21,0x26,0x2a,0x2f,0x37,0x3b,0x40,0x45,0x4a,0x4e,0x52,
0x57,0x5c,0x60,0x65,0x69,0x6e,0x73,0x77,0x7c,0x80,0x89,0x8d,0x92,0x97,0x9b,
0x9f,0xa4,0xa9,0xae,0xb1,0xb6,0xbb,0xc0,0xc4,0xc8,0xcd,0xd3,0xda,0xdf,0xe3,
0xe8,0xec,0xf1,0xf6,0xfa];

func sho_surf3d(ngrid)
{
  extern sho_scene;

  xx= span(0.0,3.0*pi,ngrid)(,-:1:ngrid);
  yy= span(0.0,4.0*pi,ngrid)(-:1:ngrid,);
  f=  sin(xx)+cos(yy);
  xyz= transpose([xx,yy,f], [3,1,2]);
  // get zone centered normals
  normzc= get_normal3d(xyz);
  // point center the normals
  norm= xyz;
  norm(1,..)= ptcen(normzc(1,..));
  norm(2,..)= ptcen(normzc(2,..));
  norm(3,..)= ptcen(normzc(3,..));
  color= [0.3,1.0,0.0];
  if(!sho_scene) return;
  clear3d;
  plsurf3d, xyz, norm, color;
  center= [0.5*(min(xx)+max(xx)), 0.5*(min(yy)+max(yy)), 0.0];
  up= [0.0,0.0,1.0];
  dist= 1.5*(min(xx)+max(xx));
  eye= center+dist*[1.0, 0.0, 0.4];
  lookat3d,eye,center,up;
  draw3d,1;
}

func sho_pic3d(npt,use_dir)
{
  clr_count;
  res= mak_pic3d(npt,use_dir);
  if(!is_void(res)) {
    write,"number of points is ",long(res(1));
    write,"number of quads is ",long(res(2));
    write,"time to draw points is ",res(3);
    write,"time to build display list is ",res(4);
  }
}

func mak_pic3d(npoint,use_dir)
{
  extern n_poly_3d, n_tri_3d, n_point_3d;
  extern n_pass, stern_red, stern_green, stern_blue;
  extern sho_scene;

  // generate a set of random points inside the -1 to 1 cube.
  // assign each a velocity based on a "rotating fluid"
  // color them based on the magnitude or direction of
  // the velocity
  if(is_void(npoint)) npoint= 10000;
  pxyz= array(0.0, 3, npoint);
  pxyz(1,)= 2.0*random(npoint)-1.0;
  pxyz(2,)= 2.0*random(npoint)-1.0;
  pxyz(3,)= 2.0*random(npoint)-1.0;
  r= abs(pxyz(1,..),pxyz(2,..));
  rp= r+!r;
  phi= atan(pxyz(2,..),pxyz(1,..)+(!r));
  vxyz= array(0.0, 3, npoint);
  vxyz(1,)= sin(phi)*r;
  vxyz(2,)= -cos(phi)*r;
  vxyz(3,)= 0.0;
  if(use_dir) {
    byt_y= bytscl(phi, top=199);
    pcolor= array(0.0, 3, npoint);
    pcolor(1,)= stern_red(1+byt_y)/255.0;
    pcolor(2,)= stern_green(1+byt_y)/255.0;
    pcolor(3,)= stern_blue(1+byt_y)/255.0;
  } else {
    vel= abs(vxyz(1,..),vxyz(2,..),vxyz(3,..));
    byt_y= bytscl(vel, top=199);
    pcolor= array(0.0, 3, npoint);
    pcolor(1,)= stern_red(1+byt_y)/255.0;
    pcolor(2,)= stern_green(1+byt_y)/255.0;
    pcolor(3,)= stern_blue(1+byt_y)/255.0;
  }

  // create a sphere with radius r0
  nth= 37;
  nphi= 73;
  phi= span(0.0, 2*pi, nphi);
  csph= cos(phi);
  snph= sin(phi);
  theta= span(0.0, pi, nth);
  csth= cos(theta);
  snth= sin(theta);
  r0= 0.125;
  origin= [0.0, 0.0, 0.0];
  // set vertex order to get correct outer side for GL
  // by reversing on the phi index
  xx= r0*snth(-,)*csph(0:1:-1);
  yy= r0*snth(-,)*snph(0:1:-1);
  zz= r0*csth(-,)*array(1.0, nphi);
  sxyz= [xx, yy, zz];
  sxyz= transpose(sxyz, [3, 1, 2]);
  sxyz += origin;
  xx= yy= zz= [];
  scolor= [0.5, 1.0, 0.0];
  nquad= (nth-1)*(nphi-1);

  tstart= tfin= array(0.0, 3);
  list_tim= draw_tim= 0.0;
  back_rgb3d, [0.0,0.0,0.0];
  if(!sho_scene) return [0, 0, 0.0, 0.0];
  clear3d;
  /* draw all points in the lists */
  timer, tstart;
  plpoints3d, pxyz, pcolor;
  // plot a sphere using OpenGL Quad Strips
  for(i= 1; i < nphi; i++) {
    txyz= sxyz(,i:i+1,);
    norm= txyz/r0;
    plqstrips3d, numberof(txyz(1,1,)), txyz, scolor, norm;
  }
  timer,tfin;
  list_tim += (tfin-tstart)(3);

  center= [0.0, 0.0, 0.0];
  up= [0.0,1.0,0.0];
  dist= 3.0;

  theta= pi/18;
  view= [sin(theta), 0.0, cos(theta)];
  eye= center+dist*view;
  lookat3d,eye,center,up;
  draw3d,1;

  n_point_3d= npoint;
  n_poly_3d= nquad;
  n_tri_3d= 2*nquad;
  return [npoint, nquad, draw_tim, list_tim];
}

func sho_points3d(npt)
{
  clr_count;
  res= mak_points3d(npt);
  if(!is_void(res)) {
    write,"number of points is ",long(res(1));
    write,"number of quads is ",long(res(2));
    write,"time to draw points is ",res(3);
    write,"time to build display list is ",res(4);
  }
}

func mak_points3d(npoint)
{
  extern n_poly_3d, n_tri_3d, n_point_3d;
  extern n_pass, stern_red, stern_green, stern_blue;
  extern sho_scene;

  // generate a set of random points inside the volume
  // color them based on the value of y32
  if(is_void(npoint)) npoint= 10000;
  pxyz= array(0.0, 3, npoint);
  pxyz(1,)= 2.0*random(npoint)-1.0;
  pxyz(2,)= 2.0*random(npoint)-1.0;
  pxyz(3,)= 2.0*random(npoint)-1.0;
  r= abs(pxyz(1,..),pxyz(2,..),pxyz(3,..));
  rp= r+!r;
  theta= acos(pxyz(3,..)/rp);
  phi= atan(pxyz(2,..),pxyz(1,..)+(!r));
  y32= sin(theta)^2*cos(theta)*cos(2*phi);
  byt_y= bytscl(y32, top=199);
  pcolor= array(0.0, 3, npoint);
  pcolor(1,)= stern_red(1+byt_y)/255.0;
  pcolor(2,)= stern_green(1+byt_y)/255.0;
  pcolor(3,)= stern_blue(1+byt_y)/255.0;

  // create a sphere with radius r
  nth= 37;
  nphi= 73;
  phi= span(0.0, 2*pi, nphi);
  csph= cos(phi);
  snph= sin(phi);
  theta= span(0.0, pi, nth);
  csth= cos(theta);
  snth= sin(theta);
  r= 0.125;
  origin= [0.0, 0.0, 0.0];
  // set vertex order to get correct outer side for GL
  // by reversing on the phi index
  xx= r*snth(-,)*csph(0:1:-1);
  yy= r*snth(-,)*snph(0:1:-1);
  zz= r*csth(-,)*array(1.0, nphi);
  sxyz= [xx, yy, zz];
  sxyz= transpose(sxyz, [3, 1, 2]);
  sxyz += origin;
  xx= yy= zz= [];
  scolor= [0.5, 1.0, 0.0];
  nquad= (nth-1)*(nphi-1);

  tstart= tfin= array(0.0, 3);
  list_tim= draw_tim= 0.0;
  back_rgb3d, [0.0,0.0,0.0];
  if(!sho_scene) return [0, 0, 0.0, 0.0];
  clear3d;
  /* draw all points in the lists */
  timer, tstart;
  plpoints3d, pxyz, pcolor;
  // plot a sphere using OpenGL Quad Strips
  for(i= 1; i < nphi; i++) {
    txyz= sxyz(,i:i+1,);
    norm= txyz/r;
    plqstrips3d, numberof(txyz(1,1,)), txyz, scolor, norm;
  }
  timer,tfin;
  list_tim += (tfin-tstart)(3);

  center= [0.0, 0.0, 0.0];
  up= [0.0,1.0,0.0];
  dist= 3.0;

  theta= pi/18;
  view= [sin(theta), 0.0, cos(theta)];
  eye= center+dist*view;
  lookat3d,eye,center,up;
  draw3d,1;

  n_point_3d= npoint;
  n_poly_3d= nquad;
  n_tri_3d= 2*nquad;
  return [npoint, nquad, draw_tim, list_tim];
}

func sho_glyph(npt)
{
  clr_count;
  res= mak_glyph(npt);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    write,"number of quads is ",long(res(2));
    write,"time to draw glyphs is ",res(3);
    write,"time to build display list is ",res(4);
  }
}

func mak_glyph(npts)
{
  extern n_poly_3d, n_tri_3d, n_point_3d;
  extern n_pass, stern_red, stern_green, stern_blue;
  extern sho_scene;

  // generate a set of points on a regular cylindrical grid
  // the grid runs from -zmx to zmx in z, 0 to 2*pi in phi
  // and 0 to rmx in r.
  // The vector function is the normal to families of
  // curves with sine variations:
  // 0 = r-a*(1+b*sin(0.5*pi*z/zmx))
  // so the function is:
  // r-hat + (0.5*pi*b*r/zmx)*cos(0.5*pi*z/zmx)/(1+b*sin(0.5*pi*z/zmx))
  if(is_void(npts)) npts= 32;
  nz= max(3, npts/5);
  nrad= max(3, nz/2);
  nphi= max(3, npts/3);
  nglyph= nz*nphi*nrad;
  if(is_void(zmx))  zmx= 1.0;
  if(is_void(amx))  amx= zmx;
  if(is_void(b))    b= 0.5;
  phi= span(0.0, 2*pi, nphi)(-:1:nz, -:1:nrad, );
  z= span(-zmx, zmx, nz)(,-:1:nrad, -:1:nphi);
  a= (span(1.0, nrad, nrad)*amx/nrad)(-:1:nz,,-:1:nphi);
  r= (1.0+b*sin(0.5*pi*z/zmx))*a;
  pxyz= array(0.0, 3, nz, nrad, nphi);
  pxyz(1,..)= r*cos(phi);
  pxyz(2,..)= r*sin(phi);
  pxyz(3,..)= z;
  f= pxyz;
  f(1,..)= cos(phi);
  f(2,..)= sin(phi);
  f(3,..)= 0.5*pi*b/zmx*a*cos(0.5*pi*z/zmx);
  fsz= abs(f(1,..), f(2,..), f(3,..));
  f /= fsz(-,..);
  byt_y= char(10)+bytscl(fsz, top=189);
  pcolor= pxyz;
  pcolor(1,..)= stern_red(1+byt_y)/255.0;
  pcolor(2,..)= stern_green(1+byt_y)/255.0;
  pcolor(3,..)= stern_blue(1+byt_y)/255.0;
  // set the height and base of a glyph
  scal= array(0.07, nglyph);
  // set theta and phi (the glyph orientation)
  thglyph= atan(abs(f(1,..),f(2,..)), f(3,..));
  phiglyph= phi;
  ntri= nglyph*4;

  // create a sphere with radius r0
  nth= 13;
  nphi= 25;
  phi= span(0.0, 2*pi, nphi);
  csph= cos(phi);
  snph= sin(phi);
  theta= span(0.0, pi, nth);
  csth= cos(theta);
  snth= sin(theta);
  r0= 0.125;
  origin= [0.0, 0.0, 0.0];
  // set vertex order to get correct outer side for GL
  // by reversing on the phi index
  xx= r0*snth(-,)*csph(0:1:-1);
  yy= r0*snth(-,)*snph(0:1:-1);
  zz= r0*csth(-,)*array(1.0, nphi);
  sxyz= [xx, yy, zz];
  sxyz= transpose(sxyz, [3, 1, 2]);
  sxyz += origin;
  xx= yy= zz= [];
  scolor= [0.5, 1.0, 0.0];
  nquad= (nth-1)*(nphi-1);

  tstart= tfin= array(0.0, 3);
  list_tim= draw_tim= 0.0;
  back_rgb3d, [0.0,0.0,0.0];
  if(!sho_scene) return [0, 0, 0.0, 0.0];
  clear3d;
  /* draw all points in the lists */
  timer, tstart;
  plglyphs3d, pxyz(,*), scal(*), thglyph(*), phiglyph(*), pcolor(,*);
  // plot a sphere using OpenGL Quad Strips
  for(i= 1; i < nphi; i++) {
    txyz= sxyz(,i:i+1,);
    norm= txyz/r0;
    plqstrips3d, numberof(txyz(1,1,)), txyz, scolor, norm;
  }
  timer,tfin;
  list_tim += (tfin-tstart)(3);

  center= [0.0, 0.0, 0.0];
  up= [0.0,1.0,0.0];
  dist= 3.0;

  theta= pi/18;
  view= [sin(theta), 0.0, cos(theta)];
  eye= center+dist*view;
  lookat3d,eye,center,up;
  draw3d,1;

  n_point_3d= 0;
  n_poly_3d= nquad+ntri;
  n_tri_3d= 2*nquad+ntri;
  return [ntri, nquad, draw_tim, list_tim];
}

func one_glyph(npts)
{
  extern sho_scene;

  // create a sphere with radius r0
  nth= 13;
  nphi= 25;
  phi= span(0.0, 2*pi, nphi);
  csph= cos(phi);
  snph= sin(phi);
  theta= span(0.0, pi, nth);
  csth= cos(theta);
  snth= sin(theta);
  r0= 0.125;
  origin= [0.0, 0.0, 0.0];
  // set vertex order to get correct outer side for GL
  // by reversing on the phi index
  xx= r0*snth(-,)*csph(0:1:-1);
  yy= r0*snth(-,)*snph(0:1:-1);
  zz= r0*csth(-,)*array(1.0, nphi);
  sxyz= [xx, yy, zz];
  sxyz= transpose(sxyz, [3, 1, 2]);
  sxyz += origin;
  xx= yy= zz= [];
  scolor= [0.5, 1.0, 0.0];
  nquad= (nth-1)*(nphi-1);

  pxyz= [[2*r0, 0.0, 0.0]];
  scal= [0.1];
  thglyph= [0.0];
  phiglyph= [0.0];
  pcolor= [[0.0, 0.3, 1.0]];

  back_rgb3d, [0.0,0.0,0.0];
  if(!sho_scene) return;
  clear3d;
  plglyphs3d, pxyz, scal, thglyph, phiglyph, pcolor;
  // plot a sphere using OpenGL Quad Strips
  for(i= 1; i < nphi; i++) {
    txyz= sxyz(,i:i+1,);
    norm= txyz/r0;
    plqstrips3d, numberof(txyz(1,1,)), txyz, scolor, norm;
  }

  center= [0.0, 0.0, 0.0];
  up= [0.0,1.0,0.0];
  dist= 3.0;

  theta= pi/18;
  view= [sin(theta), 0.0, cos(theta)];
  eye= center+dist*view;
  lookat3d,eye,center,up;
}


func sho_iso2(ngrid)
{
  clr_count;
  res= mak_iso2(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_iso2(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern tris2;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  var= get_isoreg(nx,ny,nz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  colr2= [0.6,0.5,0.0,0.3];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3cenreg(origin,delta,var,0.6*level,colr);
  for(np= 1; np <= n_pass; np++) tris2= iso3cenreg(origin,delta,var,level,colr2);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  back_rgb3d, [0.5,0.5,0.9];
  clear3d;
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    timer, tstart;
    /* draw all opaque triangles in the lists */
    if(tris) {
      res= pltrilists3d(tris);
      if(!is_void(res)) {
        numTri += res(1);
        nStrip += res(2);
        nVert += res(3);
      }
    }
    /* draw all the translucent triangles */
    if(tris2) {
      tris3= CollapseTri(tris2);
      tris4= SortTri(tris3);
      res= pltrilists3d(tris4);
      if(!is_void(res)) {
        numTri += res(1);
        nStrip += res(2);
        nVert += res(3);
      }
    }
    timer, tfin;
    list_tim += (tfin-tstart)(3);

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}

func sho_isoreg(ngrid)
{
  clr_count;
  res= mak_isoreg(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isoreg(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  var= get_isoreg(nx,ny,nz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3cenreg(origin,
      delta,var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}

func sho_isozcenreg(ngrid)
{
  clr_count;
  res= mak_isozcenreg(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isozcenreg(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  var= get_isoreg(nx,ny,nz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3zcenreg(origin,delta,
      var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}

func sho_isozcenregngrd(ngrid)
{
  clr_count;
  res= mak_isozcenregngrd(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isozcenregngrd(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  var= get_isoreg(nx,ny,nz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3zcenregngrd(origin,delta,
      var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}

func sho_isoregndx(ngrid)
{
  clr_count;
  res= mak_isoregndx(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isoregndx(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  var= get_isoreg(nx,ny,nz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3cenregndx(origin,
      delta,var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}

func sho_isozcenregndx(ngrid)
{
  clr_count;
  res= mak_isozcenregndx(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isozcenregndx(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  var= get_isoreg(nx,ny,nz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3zcenregndx(origin,
      delta,var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}



func sho_isocrv(ngrid)
{
  clr_count;
  res= mak_isocrv(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isocrv(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  xyz= get_crvgrid(nx,ny,nz);
  var= get_isocrv(xyz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3cencrv(xyz,var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}

func sho_isozcencrv(ngrid)
{
  clr_count;
  res= mak_isozcencrv(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isozcencrv(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  xyz= get_crvgrid(nx,ny,nz);
  var= get_isocrv(xyz);
  var= var(zcen,zcen,zcen);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3zcencrv(xyz,var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}



func sho_isocrvndx(ngrid)
{
  clr_count;
  res= mak_isocrvndx(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isocrvndx(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  xyz= get_crvgrid(nx,ny,nz);
  var= get_isocrv(xyz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3cencrvndx(xyz,var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}

func sho_isozcencrvndx(ngrid)
{
  clr_count;
  res= mak_isozcencrvndx(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isozcencrvndx(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  xyz= get_crvgrid(nx,ny,nz);
  var= get_isocrv(xyz);
  var= var(zcen,zcen,zcen);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3zcencrvndx(xyz,var,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}


func sho_isohex(ngrid)
{
  clr_count;
  res= mak_isohex(ngrid);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to extract iso-surface is ",res(3);
    write,"time to draw iso-surface is ",res(4);
    write,"time to build display list is ",res(5);
  }
}

func mak_isohex(ngrid)
{
  extern n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 3) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(n_pass)) n_pass= 1;
  xyz= array(0.0, 3, nx,ny,nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,ny)(-,);
  xyz(3,..)= span(-1,1,nz)(-,-,);
  r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
  rp= r+!r;
  theta= acos(xyz(3,..)/rp);
  phi= atan(xyz(2,..),xyz(1,..)+(!r));
  y32= sin(theta)^2*cos(theta)*cos(2*phi);
  f= r*(1.+y32);
  /* compute the gradient */
  dfx= xyz(1,..)/rp*(1.0+2*xyz(3,..)*(xyz(3,..)^2+2*xyz(2,..)^2)/rp^3);
  dfy= xyz(2,..)/rp*(1.0-2*xyz(3,..)*(xyz(3,..)^2+2*xyz(1,..)^2)/rp^3);
  dfz= xyz(3,..)/rp+(xyz(1,..)^2-xyz(2,..)^2)*(xyz(1,..)^2+xyz(2,..)^2-xyz(3,..)^2)/rp^4;
  grad= xyz;
  grad(1,..)= dfx;
  grad(2,..)= dfy;
  grad(3,..)= dfz;
  r= theta= phi= dfx= dfy= dfz= [];
  /* Convert from 3D to 1D coord etc. arrays */
  xyz= xyz(,*);
  grad= grad(,*);
  f= f(*);

  /* The iso-surface function wants the grid to be specified by indices 
     into point, gradient, and variable arrays. All zones are made up
	 of the same offsets from the lowest numbered corner of the zone. */
  nzone= (nx-1)*(ny-1)*(nz-1);
  hexndx= array(0, 8, nzone);
  offsets= [0, 1, nx+1, nx, nx*ny, nx*ny+1, nx*ny+nx+1, nx*ny+nx];
  for(k= 1; k < nz; k++) {
    for(j= 1; j < ny; j++) {
      for(i= 1; i < nx; i++) {
        idzn= 1+(i-1)+(nx-1)*(j-1)+(nx-1)*(ny-1)*(k-1);
        idpt= 1+(i-1)+nx*(j-1)+nx*ny*(k-1);
        hexndx(,idzn)= idpt+offsets;
      }
    }
  }

  write, "   test uses "+pr1(nzone)+" cells";

  numTri= nStrip= nVert= 0;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  level= 1.0;
  colr= [0.0,1.0,0.2];
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tris= iso3hex(xyz,grad,hexndx,
      f,0.6*level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tristrips in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim += (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  return [numTri/n_pass, nStrip/n_pass, iso_tim/n_pass, draw_tim/n_pass, 
          list_tim/n_pass];
}


func draw_polys(polys, color)
{
  nv= *(polys.nVerts);
  xyz= *(polys.xyzverts);
  norm= *(polys.normals);
  plpoly3d, nv, xyz, color, norm;
}

func sho_isotree(ngrid,cpervrt=)
{
  clr_count;
  res= do_isotree(ngrid,cpervrt=cpervrt);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to build octtree is ",res(3);
    write,"time to extract iso-surface is ",res(4);
    write,"time to draw iso-surface is ",res(5);
    write,"time to build display list is ",res(6);
  }
}

func do_isotree(ngrid,cpervrt=)
{
  extern tree, n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(cpervrt)) cpervrt= 0;
  if(is_void(n_pass)) n_pass= 1;
  var= get_isoreg(nx,ny,nz);
  if(cpervrt) {
    var2= get_var2reg(nx,ny,nz);
    vmax= max(var2);
    vmin= min(var2);
  }

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides. This means that the volume
     being contoured is a bit less than -1 to 1.
  */
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 0.6;
  colr= [0.0,1.0,0.2];
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tree= mak_isotree(var);
  timer,tfin;
  tree_tim= (tfin-tstart)(3);
  timer, tstart;
  if(cpervrt) {
    for(np= 1; np <= n_pass; np++) tris= iso3_tree(origin,delta,
        var,level,colr,tree,var2=var2);
  } else {
    for(np= 1; np <= n_pass; np++) tris= iso3_tree(origin,delta,
        var,level,colr,tree);
  }
  timer,tfin;
  iso_tim= (tfin-tstart)(3);
  if(tris && cpervrt) {
    triptr= &tris;
    while(1) {
      if(!triptr) break;
      if(is_void(*triptr)) break;
      if(triptr->numTri) {
        /* one color per per vertex based on auxiliary variable */
        ndx= 1 + (long) (gl_ncolr*(*(triptr->var2)-vmin)/(vmax-vmin));
        ndx= min(ndx, gl_ncolr);
        nclr= numberof(ndx)/3;
        colr= array(float, 3, 3, nclr);
        colr(1,,)= gl_rr(ndx)/256.0;
        colr(2,,)= gl_gg(ndx)/256.0;
        colr(3,,)= gl_bb(ndx)/256.0;
        triptr->colors= &colr;
      }
      triptr= triptr->next;
    } ;
  }
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tri arrays in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim= (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  sizes= [nx, ny, nz];
  return [numTri/n_pass, nStrip/n_pass, tree_tim/n_pass, iso_tim/n_pass, 
	        draw_tim/n_pass, list_tim/n_pass];
}

func sho_isotreevarr(ngrid,cpervrt=)
{
  clr_count;
  res= do_isotreevarr(ngrid,cpervrt=cpervrt);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    if(res(1) != res(2)) {
      write,"number of triangles per strip is",res(1)/double(res(2));
    }
    write,"time to build octtree is ",res(3);
    write,"time to extract iso-surface is ",res(4);
    write,"time to draw iso-surface is ",res(5);
    write,"time to build display list is ",res(6);
  }
}

func do_isotreevarr(ngrid,cpervrt=)
{
  extern tree, n_poly_3d, n_tri_3d, n_strip_3d, n_pass;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ny= nz= ngrid;
  }
  if(is_void(cpervrt)) cpervrt= 0;
  if(is_void(n_pass)) n_pass= 1;
  var= get_isoreg(nx,ny,nz);
  if(cpervrt) {
    var2= get_var2reg(nx,ny,nz);
    vmax= max(var2);
    vmin= min(var2);
  }

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides. This means that the volume
     being contoured is a bit less than -1 to 1.
  */
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  level= 0.6;
  colr= [0.0,1.0,0.2];
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  for(np= 1; np <= n_pass; np++) tree= mak_isotree(var);
  timer,tfin;
  tree_tim= (tfin-tstart)(3);
  timer, tstart;
  if(cpervrt) {
    for(np= 1; np <= n_pass; np++) tris= iso3_treevarr(origin,delta,
        var,level,colr,tree,var2=var2);
  } else {
    for(np= 1; np <= n_pass; np++) tris= iso3_treevarr(origin,delta,
        var,level,colr,tree);
  }
  timer,tfin;
  iso_tim= (tfin-tstart)(3);
  if(tris && cpervrt) {
    triptr= &tris;
    while(1) {
      if(!triptr) break;
      if(is_void(*triptr)) break;
      if(triptr->numTri) {
        /* one color per per vertex based on auxiliary variable */
        ndx= 1 + (long) (gl_ncolr*(*(triptr->var2)-vmin)/(vmax-vmin));
        ndx= min(ndx, gl_ncolr);
        nclr= numberof(ndx);
        colr= array(float, 3, nclr);
        colr(1,)= gl_rr(ndx(*))/256.0;
        colr(2,)= gl_gg(ndx(*))/256.0;
        colr(3,)= gl_bb(ndx(*))/256.0;
        triptr->colors= &colr;
      }
      triptr= triptr->next;
    } ;
  }
  if(!sho_scene) return [1.0, 1.0, 0.0, 0.0, 0.0];
  for(np= 1; np <= n_pass; np++) {
    clear3d;
    /* draw all tri arrays in the lists */
    timer, tstart;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer,tfin;
    list_tim= (tfin-tstart)(3);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }

    timer, tstart;
    draw3d, 1;
    if(is_void(keepview)) {
      stdview3d;
    } else if(!keepview) {
      stdview3d;
    }
    timer, tfin;
    draw_tim= (tfin-tstart)(3);
  }
  n_tri_3d= numTri/n_pass;
  n_poly_3d= numTri/n_pass;
  sizes= [nx, ny, nz];
  return [numTri/n_pass, nStrip/n_pass, tree_tim/n_pass, iso_tim/n_pass, 
	        draw_tim/n_pass, list_tim/n_pass];
}

func sho_volviz3(ngrid)
{
  extern sho_scene;

  if(is_void(ngrid)) ngrid= 6;
  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ngrid;
    ny= ngrid;
    nz= ngrid;
  }
  var= get_isoreg(nx,ny,nz);
  maxvar= max(var);
  /* create functions that are large where f is near
     interesting values */
  var1= 0.25*maxvar;
  var2= 0.5*maxvar;
  var3= 0.75*maxvar;
  dv= 0.07*maxvar;
  v1= exp(-((var-var1)/dv)^2);
  v2= exp(-((var-var2)/dv)^2);
  v3= exp(-((var-var3)/dv)^2);
  texout= array(char,4,nx,ny,nz);
  alpha= v1+v2+v3;
  alpha /= max(alpha);
  /* a guess at the alpha needed to make each surface
     just barely opaque */
  bigalpha= 0.07/dv*14.0/nx*255;
  /* alpha is large where any of the three functions are large */
  texout(4,..)= char(bigalpha*alpha);
  /* turn on specific colors where one of the funcions is large */
  texout(1,..)= char(255*v1);
  texout(2,..)= char(255*v2);
  texout(3,..)= char(255*v3);
  f= [];

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  back_rgb3d, [0.0,0.0,0.0];
  lookat3d, [0.0, 0.0, 2.0], [0.5, 0.5, 0.5], [0.0, 1.0, 0.0];
  clear3d;

  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  delta= [1.0/(nx-1.0), 1.0/(ny-1.0), 1.0/(nz-1.0)];
  if(!sho_scene) return [0.0, 0.0];
  /* draw the volume visualization */
  pltex2dvol,delta,texout;
  timer, tfin;
  list_tim= (tfin-tstart)(3);

  timer, tstart;
  draw3d, 1;
  draw3d_trigger;
  timer, tfin;
  draw_tim= (tfin-tstart)(3);
  write,"time to build display list is ",list_tim;
  write,"time to draw is ",draw_tim;
  return [list_tim, draw_tim];
}

func sho_texiso(ngrid,frac=)
{
  extern sho_scene;

  ASSERT, has_tex3d(), "3D textures are not available.";
  make_strip= 0;
  if(is_void(ngrid)) ngrid= 10;
  if(ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ngrid;
    ny= ngrid;
    nz= ngrid;
  }
  nslab= max(nx,ny,nz);
  xyz= array(0.0, 3, nx,ny,nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,ny)(-,);
  xyz(3,..)= span(-1,1,nz)(-,-,);
  r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
  theta= acos(xyz(3,..)/(r+1.0e-10));
  phi= atan(xyz(2,..),xyz(1,..)+(!r));
  y32= 1.5*sqrt(3.0)*sin(theta)^2*cos(theta)*cos(2*phi);
  y21= 2.0*sin(theta)*cos(theta)*cos(phi);
  y22= sin(theta)^2;
  exr= exp(-3.0*r/max(r));
  r1= 0.2*max(r);
  r2= 0.4*max(r)
  dr= 0.07*max(r);
  f32= (1.0+y32)*exp( -((r-r1)/dr)^2 );
  f21= (1.0+y21)*exp( -((r-r2)/dr)^2 );
  mxr= max(r);
  f22= 4.0*r/mxr*(1.0-r/mxr)*y22;
  f= f32+f21;
  df= 0.07*max(f);
  r= theta= phi= [];
  maxvar= max(f);
  texout= array(char,4,nx,ny,nz);
  alpha= f/max(f);
  /* a guess at the alpha needed to make each surface
     just barely opaque */
  bigalpha= 0.07/df*14.0/nslab*255;
  /* alpha is large where any of the three functions are large */
  texout(4,..)= char(bigalpha*alpha);
  /* turn on specific colors where one of the funcions is large */
  texout(3,..)= char(255*f32/max(f32));
  texout(1,..)= 0;
  texout(2,..)= char(255*f21/max(f21));
  f= [];

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  len= [2.0,2.0,2.0];
  delta= [2.0/(nx-1.0),2.0/(ny-1.0),2.0/(nz-1.0)];
  origin= [-1.0,-1.0,-1.0];  /* location of the first real cell */
  if(is_void(frac)) {
    level= 0.12*max(f22);
  } else {
    level= frac*max(f22);
  }

  back_rgb3d, [0.0,0.0,0.0];
  lookat3d, [0.0, 0.0, 2.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0];
  clear3d;

  numTri= nStrip= nVert= 0;
  colr= [0.0, 0.0, 1.0];
  iso_tim= list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  tris= iso3cenreg(origin,delta,f22,level,colr);
  timer,tfin;
  iso_tim += (tfin-tstart)(3);

  tstart= tfin= array(0.0, 3);
  timer, tstart;
  if(!sho_scene) return [0.0, 0.0];
  /* draw iso-surface using all triangles */
  res= pltrilists3d(tris);
  if(!is_void(res)) {
    numTri += res(1);
    nStrip += res(2);
    nVert += res(3);
  }
  timer, tfin;
  list_tim += (tfin-tstart)(3);

  timer, tstart;
  /* draw the volume visualization */
  pltex3dvol, nslab, len, texout, origin=origin;
  timer, tfin;
  list_tim += (tfin-tstart)(3);

  timer, tstart;
  draw3d, 1;
  timer, tfin;
  draw_tim= (tfin-tstart)(3);
  write,"time to build display list is ",list_tim;
  write,"time to draw is ",draw_tim;
  return [list_tim, draw_tim];
}

func mak_ghost(f)
{
  dimf= dimsof(f);
  nx= dimf(2);
  ny= dimf(3);
  nz= dimf(4);
  /* wrap ghost cells around the data */
  fg= array(0.0, nx+2, ny+2, nz+2);
  fg(2:-1,2:-1,2:-1)= f;
  /* extrapolate from the 6 faces to get more accurate ghosts */
  fg(1,2:-1,2:-1)= 2*f(1,..)-f(2,..);
  fg(0,2:-1,2:-1)= 2*f(0,..)-f(-1,..);
  fg(2:-1,1,2:-1)= 2*f(,1,)-f(,2,);
  fg(2:-1,0,2:-1)= 2*f(,0,)-f(,-1,);
  fg(2:-1,2:-1,1)= 2*f(..,1)-f(..,2);
  fg(2:-1,2:-1,0)= 2*f(..,0)-f(..,-1);
  /* extrapolate to the 12 edges */
  fg(1,1,2:-1)= 2.0*f(1,1,)-f(2,2);
  fg(0,1,2:-1)= 2.0*f(0,1,)-f(-1,2,);
  fg(1,0,2:-1)= 2.0*f(1,0,)-f(2,-1,);
  fg(0,0,2:-1)= 2.0*f(0,0,)-f(-1,-1,);
  fg(1,2:-1,1)= 2.0*f(1,,1)-f(2,,2);
  fg(0,2:-1,1)= 2.0*f(0,,1)-f(-1,,2);
  fg(1,2:-1,0)= 2.0*f(1,,-1)-f(2,,-1);
  fg(0,2:-1,0)= 2.0*f(-1,,0)-f(-1,,-1);
  fg(2:-1,1,1)= 2.0*f(,1,1)-f(,2,2);
  fg(2:-1,0,1)= 2.0*f(,0,1)-f(,-1,2);
  fg(2:-1,1,0)= 2.0*f(,1,0)-f(,2,-1);
  fg(2:-1,0,0)= 2.0*f(,0,0)-f(,-1,-1);
  /* interpolate to the eight corners */
  fg(0,0,0)= (fg(-1,0,0)+fg(0,-1,0)+fg(0,0,-1))/3.0;
  fg(1,0,0)= (fg(2,0,0)+fg(1,-1,0)+fg(1,0,-1))/3.0;
  fg(0,1,0)= (fg(-1,1,0)+fg(0,2,0)+fg(0,1,-1))/3.0;
  fg(1,1,0)= (fg(2,1,0)+fg(1,2,0)+fg(1,1,-1))/3.0;
  fg(0,0,1)= (fg(-1,0,1)+fg(0,-1,1)+fg(0,0,2))/3.0;
  fg(1,0,1)= (fg(2,0,1)+fg(1,-1,1)+fg(1,0,2))/3.0;
  fg(0,1,1)= (fg(-1,1,1)+fg(0,2,1)+fg(0,1,2))/3.0;
  fg(1,1,1)= (fg(2,1,1)+fg(1,2,1)+fg(1,1,2))/3.0;
  return fg;
}

func clr_count
{
  extern n_tri_3d, n_poly_3d, n_strip_3d, n_point_3d, n_line_3d;

  n_tri_3d= n_poly_3d= n_strip_3d= n_point_3d= n_line_3d= 0;
}

func do_rot(num,theta,phi,dtheta,dphi,targ_tim=)
{
  extern n_tri_3d, n_poly_3d, n_strip_3d, n_point_3d, n_line_3d, do_rotate;

  if(!do_rotate) return;
  if(is_void(num)) num= 50;
  if(is_void(theta)) theta= pi/2;
  if(is_void(phi)) phi= 0.0;
  if(is_void(dtheta)) dtheta= pi/10.0;
  if(is_void(dphi)) dphi= pi/18;

  eye= center= array(0.0, 3);
  get_center3d,center;
  get_eye3d,eye;
  up= [0.0,0.0,1.0];
  dist= sqrt( sum( (eye-center)^2 ) );

  elapsed= [0.,0.,0.];
  timer, elapsed;
  elapsed0= elapsed;
  for(i= 1; i <= num; i++) {
    phi += dphi;
    thet= theta+dtheta*sin(2.3*phi);
    csth= cos(thet);
    snth= sin(thet);
    view= [snth*cos(phi), snth*sin(phi), csth];
    eye= center+dist*view;
    lookat3d,eye,center,up;
    draw3d, 1;
    if(targ_tim) {
      timer, elapsed;
      if(elapsed(3)-elapsed0(3) > targ_tim) {
        num= i;
        break;
      }
    }
  }

  timer, elapsed;
  rate= num/(elapsed(3)-elapsed0(3)+1.0e-10);
  write,"display rate is "+pr1(rate)+" frames per second";
  if(n_poly_3d) write,pr1(n_poly_3d*rate)+" polygons per second";
  if(n_tri_3d) write,pr1(n_tri_3d*rate)+" equiv. triangles per second";
  if(n_poly_3d) write,pr1(n_poly_3d)+" polygons drawn per frame";
  if(n_tri_3d && n_strip_3d) write,"average number of triangles per strip is "+
      pr1(n_tri_3d/double(n_strip_3d));
  if(n_line_3d) write,pr1(n_line_3d*rate)+" line segments per second";
  if(n_point_3d) write,pr1(n_point_3d*rate)+" points per second";
}

func mak_slice(ngrid)
{
  extern gl_rr, gl_gg, gl_bb, gl_ncolr, gl_ctab;
  extern n_poly_3d, n_tri_3d;
  extern sho_scene;

  if(is_void(ngrid) || ngrid < 6) {
    nx= ny= nz= 64;
  } else {
    nx= ngrid;
    ny= ngrid;
    nz= ngrid;
  }
  numTri= 0;
  nStrip= 0;
  nVert= 0;

  clr_count;
  palette3d,"stern.gp";
  f= get_isoreg(nx,ny,nz);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];

  level= 1.0;
  colr1= [0.0,1.0,0.2];
  colrtyp= -3;  /* three component scalar color */
  level1= 0.5;
  tris_i= iso3cenreg(origin,delta,f,level1,colr1)
  if(tris_i) {
    num= SizeTriArrays3d(&tris_i);
    newtris= TriArrayGrp(
      numTri= 0,
      cellIDs= &array(0, num),
      xyzverts= &array(0.0, 3, 3, num),
      normals= &array(0.0, 3, 3, num),
      colors= &array(float, 3, num),
      var2= &nulvar,
      triEdg= &nulvar,
      triStart= &nulvar,
      nTris= &nulvar,
      next= &nulvar
    );
    if(tris_i.var2) newtris.var2= &array(0.0, 3, num);
    CollapseTriArrays3d, colrtyp, &tris_i, &newtris;
    tris_i= newtris;
    newtris= [];
  }
  colr2= [1.0,0.2,0.0];
  colrtyp= -3;  /* three component scalar color */
  level2= 1.0;
  tris_o= iso3cenreg(origin,delta,f,level2,colr2)
  if(tris_o) {
    num= SizeTriArrays3d(&tris_o);
    newtris= TriArrayGrp(
      numTri= 0,
      cellIDs= &array(0, num),
      xyzverts= &array(0.0, 3, 3, num),
      normals= &array(0.0, 3, 3, num),
      colors= &array(float, 3, num),
      var2= &nulvar,
      triEdg= &nulvar,
      triStart= &nulvar,
      nTris= &nulvar,
      next= &nulvar
    );
    if(tris_o.var2) newtris.var2= &array(0.0, 3, num);
    CollapseTriArrays3d, colrtyp, &tris_o, &newtris;
    tris_o= newtris;
    newtris= [];
  }

  pxy= [0, 0, 1.0, 0];
  pyz= [1.0, 0, 0, 0];

  tnew_il= slice2only(pyz, tris_i);
  tnew_ir= slice2only(-pyz, tris_i);
  tnew_ir= slice2only(-pxy, tnew_ir);

  tnew_ol= slice2only(pyz, tris_o);

  norm= pyz(1:3);
  pnt= norm*pyz(4);
  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides if guard is non-zero. 
     This means that the volume being contoured is a bit less 
     than -1 to 1.
  */
  tris= slice_tree([nx,ny,nz], origin, delta, pnt, norm, f, guard=1);
  vmin= 0.0;
  vmax= max(f);
  triptr= &tris;
  while(1) {
    if(is_void(*triptr)) break;
    if(!triptr) break;
    ntri= triptr->numTri;
    if(ntri) {
      /* one color per vertex based on var2 */
      vx= *(triptr->var2);
      vx= min(vmax, max(vmin, vx));
      ndx= 1 + (long) (gl_ncolr*(vx-vmin)/(vmax-vmin));
      ndx= min(ndx, gl_ncolr);
      colr= array(float, 3, 3, numberof(ndx(1,)));
      colr(1,,)= gl_rr(ndx)/256.0;
      colr(2,,)= gl_gg(ndx)/256.0;
      colr(3,,)= gl_bb(ndx)/256.0;
      triptr->colors= &colr;
    }
    triptr= triptr->next;
  } ;

  if(!sho_scene) return;
  clear3d;
  if(tris) {
    res= pltrilists3d(tris,emit=1);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }
  }
  if(tnew_il) {
    res= pltrilists3d(tnew_il);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }
  }
  if(tnew_ir) {
    res= pltrilists3d(tnew_ir);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }
  }
  if(tnew_ol) {
    res= pltrilists3d(tnew_ol);
    if(!is_void(res)) {
      numTri += res(1);
      nStrip += res(2);
      nVert += res(3);
    }
  }
  draw3d, 1;
  stdview3d;
  n_poly_3d= numTri;
  n_tri_3d=numTri;
}

func dslice(size, numrot, dphi)
{
  if(is_void(size)) size= 40;
  mak_slice,size;
  theta= pi/3;
  phi= 0.0;
  if(is_void(dphi)) dphi= pi/18;
  if(is_void(numrot)) numrot= 50;
  do_rot,numrot,theta,phi,0.0,dphi;
}


func rsphere(num, dtheta)
{
  mak_sphere;
  theta= pi/6;
  phi= 0.0;
  if(is_void(dphi)) dphi= pi/18;
  if(is_void(num)) num= 50;
  do_rot,num,theta,phi,0.0,dphi;
}

func makppm(name)
{
  extern pix;
  nx= get_width3d();
  ny= get_hite3d();
  pix= array(char, 3, nx, ny);
  grabpix3d, nx, ny, pix;
  if(is_void(name)) {
    pnm_write, pix, "test.ppm";
  } else {
    pnm_write, pix, name;
  }
}

func getput
{
  extern pix;
  nx= get_width3d();
  ny= get_hite3d();
  pix= array(char, 3, nx, ny);
  grabpix3d, nx, ny, pix;
  clear3d;
  pix= pix(,,0:1:-1);
  putpix3d,pix;
}

func mak_sphere(nth, nphi, triangle=, alpha=, origin=, radius=, color=, kpfrm=)
{
  extern n_poly_3d, n_tri_3d;
  extern do_smooth, one_color;
  extern sho_scene;

  if(is_void(alpha)) {
    alpha= -1.0;
  } else {
    alpha= min(1.0, max(0.0, alpha));
  }
  if(is_void(do_smooth)) do_smooth= 1;
  if(is_void(triangle)) triangle= 0;
  if(is_void(one_color)) one_color= 0;
  if(is_void(color)) color= [0.5, 1.0, 0.0];
  if(is_void(kpfrm)) kpfrm= 0;
  if(!kpfrm) clear3d;

  clr_count;
  if(is_void(nth)) nth= 37;
  if(is_void(nphi)) nphi= 73;
  // make phi decrease as index increases so zones will be "walked"
  // in the correct order for OpenGL
  phi= span(2*pi, 0.0, nphi);
  csph= cos(phi);
  snph= sin(phi);
  theta= span(0.0, pi, nth);
  csth= cos(theta);
  snth= sin(theta);
  if(is_void(radius)) r= 1.0;
  else r= radius;
  if(is_void(origin)) origin= [0.0, 0.0, 0.0];
  if(one_color) {
    cbase= array(color, nphi-1);
  } else {
    if(alpha >= 0) {
      cbase= array([0.0, 0.0, 0.0, 1.0], nphi-1);
      cbase(4,)= alpha;
    } else {
      cbase= array([0.0, 0.0, 0.0], nphi-1);
    }
    cbase(2,1:nphi/2)= span(0, 1.0, nphi/2);
    cbase(2,nphi/2:nphi-1)= span(1.0, 0.0, nphi-1-nphi/2+1);
    cbase(1,nphi/2:nphi-1)= span(0, 1.0, nphi-1-nphi/2+1);
    cbase(3,1:nphi/2)= span(1.0, 0.0, nphi/2);
  }
  xx= r*snth(-,)*csph;
  yy= r*snth(-,)*snph;
  zz= r*csth(-,)*array(1.0, nphi);
  xyz= [xx, yy, zz];
  xyz= transpose(xyz, [3, 1, 2]);
  xyz += origin;
  xx= yy= zz= [];

  if(!sho_scene) return;
  // plot a sphere using OpenGL Quad Strips
  for(i= 1; i < nphi; i++) {
    if(triangle) {
      pxyz= xyz(,i:i+1,)(,*)/r;
      if(do_smooth) {
        norm= pxyz;
      } else {
        norm= pxyz(,zcen)(,*);
      }
      color= array(cbase(,i), 2*nth-2);
      pltstrips3d, numberof(pxyz(1,)), pxyz, color, norm;
    } else {
      pxyz= xyz(,i:i+1,);
      if(do_smooth) {
        norm= pxyz;
      } else {
        norm= pxyz(,zcen,zcen)(,*);
      }
      color= array(cbase(,i), nth-1);
      plqstrips3d, numberof(pxyz(1,1,)), pxyz, color, norm;
    }
  }
  if(!kpfrm) stdview3d;
  draw3d,1;

  if(triangle) {
    n_poly_3d += 2*(nth-1)*(nphi-1);
    n_tri_3d += 2*(nth-1)*(nphi-1);
  } else {
    n_poly_3d += (nth-1)*(nphi-1);
    n_tri_3d += 2*(nth-1)*(nphi-1);
  }
}

func mak_colrsurf(nth, nphi, triangle=, kpfrm=)
{
  extern n_poly_3d, n_tri_3d;
  extern sho_scene;

  // draw a sphere colored by the value of Y32 on the sphere.
  if(is_void(nth)) nth= 37;
  if(is_void(nphi)) nphi= 73;
  if(is_void(triangle)) triangle= 0;
  if(is_void(kpfrm)) kpfrm= 0;
  if(!kpfrm) clear3d;

  clr_count;

  // make phi decrease as index increases so zones will be "walked"
  // in the correct order for OpenGL
  phi= span(2*pi, 0.0, nphi);
  csph= cos(phi);
  snph= sin(phi);
  theta= span(0.0, pi, nth);
  csth= cos(theta);
  snth= sin(theta);
  // set vertex order to get correct outer side for GL
  // by reversing on the phi index
  xx= snth(-,)*csph(0:1:-1);
  yy= snth(-,)*snph(0:1:-1);
  zz= csth(-,)*array(1.0, nphi);
  xyz= [xx, yy, zz];
  xyz= transpose(xyz, [3, 1, 2]);
  xx= yy= zz= [];

  /* The surface is a unit radius sphere, so the normal is the same as
     the ray from the origin to the point on the surface. */
  norm= xyz;

  y32= snth^2*csth*cos(2*phi);
  var= 1.+y32;

  vmax= max(var);
  vmin= min(var);
  ndx= 1 + (long) (gl_ncolr*(var-vmin)/(vmax-vmin));
  ndx= min(ndx, gl_ncolr);
  colr= array(float, 3, nth, nphi);
  colr(1,,)= gl_rr(ndx)/256.0;
  colr(2,,)= gl_gg(ndx)/256.0;
  colr(3,,)= gl_bb(ndx)/256.0;

  if(!sho_scene) return;
  plcolrsurf3d, xyz, norm, colr;
  if(!kpfrm) stdview3d;
  draw3d,1;

  if(triangle) {
    n_poly_3d += 2*(nth-1)*(nphi-1);
    n_tri_3d += 2*(nth-1)*(nphi-1);
  } else {
    n_poly_3d += (nth-1)*(nphi-1);
    n_tri_3d += 2*(nth-1)*(nphi-1);
  }
}

func mak_ellipsoid(nth, nphi, eps, origin=, radius=, color=, kpfrm=)
{
  extern n_poly_3d, n_tri_3d;
  extern do_smooth, one_color;
  extern sho_scene;

  if(is_void(do_smooth)) do_smooth= 1;
  if(is_void(one_color)) one_color= 0;
  if(is_void(color)) color= [0.5, 1.0, 0.0];
  if(is_void(kpfrm)) kpfrm= 0;
  if(!kpfrm) clear3d;

  clr_count;
  if(is_void(nth)) nth= 37;
  if(is_void(nphi)) nphi= 73;
  phi= span(0.0, 2*pi, nphi);
  csph= cos(phi);
  snph= sin(phi);
  theta= span(0.0, pi, nth);
  csth= cos(theta);
  snth= sin(theta);
  if(is_void(radius)) r= 1.0;
  else r= radius;
  if(is_void(eps)) eps= 1.0;
  else if(eps <= 0.0) eps= 1.0;
  if(is_void(origin)) origin= [0.0, 0.0, 0.0];
  if(one_color) {
    cbase= array(color, nphi-1);
  } else {
    cbase= array([0.0, 0.0, 0.0], nphi-1);
    cbase(2,1:nphi/2)= span(0, 1.0, nphi/2);
    cbase(2,nphi/2:nphi-1)= span(1.0, 0.0, nphi-1-nphi/2+1);
    cbase(1,nphi/2:nphi-1)= span(0, 1.0, nphi-1-nphi/2+1);
    cbase(3,1:nphi/2)= span(1.0, 0.0, nphi/2);
  }
  // set vertex order to get correct outer side for GL
  // by reversing on the phi index
  xx= r*snth(-,)*csph(0:1:-1);
  yy= r*snth(-,)*snph(0:1:-1);
  zz= r*csth(-,)*array(1.0, nphi);
  xyz= [eps*xx, eps*yy, zz];
  xyz= transpose(xyz, [3, 1, 2]);
  norm= [xx, yy, eps*zz];
  norm= transpose(norm, [3, 1, 2]);
  ln= abs(norm(1,..), norm(2,..), norm(3,..));
  norm /= (ln+!ln)(-:1:3,..);
  xyz += origin;
  xx= yy= zz= [];

  if(!sho_scene) return;
  // plot a sphere using OpenGL Quad Strips
  for(i= 1; i < nphi; i++) {
    pxyz= xyz(,i:i+1,)(,*)/r;
    nrm= norm(,i:i+1,)(,*)/r;
    color= array(cbase(,i), 2*nth-2);
    pltstrips3d, numberof(pxyz(1,)), pxyz, color, nrm;
  }
  if(!kpfrm) stdview3d;
  draw3d,1;

  n_poly_3d += 2*(nth-1)*(nphi-1);
  n_tri_3d += 2*(nth-1)*(nphi-1);
}

func rlines(nseg,npos)
{
  extern n_line_3d;
  extern sho_scene;

  if(is_void(nseg)) nseg= 200;
  if(is_void(npos)) npos= 8;
  if(is_void(dtheta)) dtheta= pi/18;
  if(is_void(dphi)) dphi= 2.0*pi/3.0;

  clr_count;
  // generate a number of streamlines with 
  // various initial positions
  zhi= 2.0;
  sp= span(0.0, zhi, npos);
  // set angular shift between points on a streamline
  phi_hi= 1.3*dphi;
  phi= span(0.0, phi_hi, nseg);
  xyz= array(0.0, 3, nseg, 3*npos);
  // compute all streamlines
  phi0= 0.0;
  r0= 1.0;
  x= r0*cos(phi+phi0)(,-:1:npos);
  y= r0*sin(phi+phi0)(,-:1:npos);
  z= (r0*sp)(-,)+r0*phi/phi_hi;
  txyz= [x, y, z];
  xyz(,,1:npos)= transpose(txyz, [3,1,2]);

  phi1= dphi;
  r1= 1.5;
  ibase= npos;
  x= r1*cos(phi+phi1)(,-:1:npos);
  y= r1*sin(phi+phi1)(,-:1:npos);
  z= (r1*sp)(-,)+r1*phi/phi_hi;
  txyz= [x, y, z];
  xyz(,,ibase+1:ibase+npos)= transpose(txyz, [3,1,2]);

  phi2= 2*dphi;
  r2= 2.0;
  ibase= 2*npos;
  x= r2*cos(phi+phi2)(,-:1:npos);
  y= r2*sin(phi+phi2)(,-:1:npos);
  z= (r2*sp)(-,)+r2*phi/phi_hi;
  txyz= [x, y, z];
  xyz(,,ibase+1:ibase+npos)= transpose(txyz, [3,1,2]);

  if(!sho_scene) return [0.0, 0.0, 3];
  /* put the lines into the display list */
  color= [[0.4, 0.9, 0.0], [1.0, 0.0, 0.0], [0.0, 0.7, 0.7]];
  clear3d;
  list_time= draw_time= 0.0;
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  for(i= 1; i <= npos; i++) {
    pllines3d, xyz(..,i), color(,1);
  }
  for(i= npos+1; i <= 2*npos; i++) {
    pllines3d, xyz(..,i), color(,2);
  }
  for(i= 2*npos+1; i <= 3*npos; i++) {
    pllines3d, xyz(..,i), color(,3);
  }
  timer,tfin;
  list_time += tfin(3)-tstart(3);

  center= [0.0, 0.0, 0.0];
  up= [0.0,1.0,0.0];
  dist= 12.0;

  theta= 0.5*dtheta;
  view= [sin(theta), 0.0, cos(theta)];
  eye= center+dist*view;
  lookat3d,eye,center,up;
  draw3d,1;

  n_line_3d= 3*npos*nseg;
  return [list_time, draw_time, 3*npos*nseg];
}

func sho_drum(time_limit,targ_tim=)
/* DOCUMENT sho_drum
     Display a drum head moving in 3D.
     The drumhead is initially stationary, but has a bump near one
     edge.  Yorick is solving a 2D wave equation to compute the
     evolution of this bump.

     There is an optional argument to sho_drum: the 
     time limit on the duration of the movie in seconds (default
     is 60 seconds).
     
     The size of the grid is set by the external variable drum_size.
 */
{
  require, "movie.i";
  extern sho_scene;

  maxfr= 200;
  /* generate a num-by-num cell mesh on the [-1,1] square */
  if(is_void(num)) num= 30;
  x= span(-1, 1, num+1)(,-:1:num+1);
  y= transpose(x);
  /* map the square mesh into a mesh on the unit circle
     this mesh has more nearly equal area cells than a polar
     coordinate circle */
  scale= max(abs(y),abs(x))/(abs(y,x)+1.e-30);
  /* note that abs(y,x)=sqrt(x^2+y^2) */
  x*= scale;
  y*= scale;

  f= f0= exp(-8.*abs(y+.67,x+.25)^2)*(1.-abs(y,x)^2);
  fdot= 0.0*f(2:-1,2:-1);

  af= abs(f(2:-1,2:-1));
  lf= laplacian(f, y,x);
  /* set the "Courant" lmit */
  xdz= x(dif,zcen);
  xzd= x(zcen,dif);
  ydz= y(dif,zcen);
  yzd= y(zcen,dif);
  dt= 0.375*sqrt(min(abs(xdz*yzd - xzd*ydz)));

  mcolor= [1.0,1.0,0.0];

  lookat3d,[1.75,0.0,0.77],[0.0,0.0,0.0], [-0.401,-0.043,0.915];
  clear3d;

  elapsed= [0.,0.,0.];
  timer, elapsed;
  elapsed0= elapsed;

  /* roll the movie */
  f= f0;
  for(i= 0; i < maxfr; i++) {
    /* display first */
    xyz= [x, y, f];
    xyz= transpose(xyz,2);
    /* compute surface normals at the vertices */
    norm= get_normal3d(xyz);
    norm= norm(,pcen,pcen);

    if(sho_scene) {
      clear3d;
      plsurf3d,xyz,norm,mcolor;
      draw3d,1;
      if(targ_tim) {
        timer, elapsed;
        if(elapsed(3)-elapsed0(3) > targ_tim) {
          break;
        }
      }
    }
  
    /* then take a step forward in time */
    lf= laplacian(f, y, x);
    af= abs(f(2:-1,2:-1));
    fdot+= lf*dt;
    f(2:-1,2:-1)+= fdot*dt;
  }
  timer, elapsed;
  rate= i/(elapsed(3)-elapsed0(3)+1.0e-10);
  write,"display rate is "+pr1(rate)+" frames segments per second";
}

func laplacian(f, y,x)
{
  /* There are many ways to form the Laplacian as a finite difference.
     This one is nice in Yorick.  */
  /* Start with the two median vectors across each zone.  */
  fdz= f(dif,zcen);
  fzd= f(zcen,dif);
  xdz= x(dif,zcen);
  xzd= x(zcen,dif);
  ydz= y(dif,zcen);
  yzd= y(zcen,dif);

  /* Estimate the gradient at the center of each cell.  */
  area= xdz*yzd - xzd*ydz;
  gradfx= (fdz*yzd - fzd*ydz)/area;
  gradfy= (xdz*fzd - xzd*fdz)/area;

  /* Now consider the mesh formed by the center points of the original.  */
  x= x(zcen,zcen);
  y= y(zcen,zcen);
  xdz= x(dif,);
  xzd= x(,dif);
  ydz= y(dif,);
  yzd= y(,dif);
  area= xdz(,zcen)*yzd(zcen,) - xzd(zcen,)*ydz(,zcen);

  return ((xdz*gradfy(zcen,)-ydz*gradfx(zcen,))(,dif) +
	  (yzd*gradfx(,zcen)-xzd*gradfy(,zcen))(dif,)) / area;
}


func get_isoreg(nx,ny,nz)
{
  // in some special cases the variable may have
  // been pre-computed and saved in a file
  if(nx == 80 && ny == 80 && nz == 80) {
    var= get_y32(80);
  } else if(nx == 100 && ny == 100 && nz == 100) {
    var= get_y32(100);
  } else if(nx == 120 && ny == 120 && nz == 120) {
    var= get_y32(120);
  } else if(nx == 150 && ny == 150 && nz == 150) {
    var= get_y32(150);
  } else {
    xyz= array(0.0, 3, nx,ny,nz);
    xyz(1,..)= span(-1,1,nx);
    xyz(2,..)= span(-1,1,ny)(-,);
    xyz(3,..)= span(-1,1,nz)(-,-,);
    r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
    rp= r+!r;
    theta= acos(xyz(3,..)/rp);
    phi= atan(xyz(2,..),xyz(1,..)+(!r));
    y32= sin(theta)^2*cos(theta)*cos(2*phi);
    var= r*(1.+y32);
  }
  return var;
}

func get_var2reg(nx,ny,nz)
{
  xyz= array(0.0, 3, nx,ny,nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,ny)(-,);
  xyz(3,..)= span(-1,1,nz)(-,-,);
  r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
  rp= r+!r;
  theta= acos(xyz(3,..)/rp);
  phi= atan(xyz(2,..),xyz(1,..)+(!r));
  var2= (sin(theta)^2+cos(theta))*cos(4*phi);
  return var2;
}

func get_crvgrid(nx,ny,nz)
{
  xyz= array(0.0, 3, nx,ny,nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,ny)(-,);
  xyz(3,..)= span(-1,1,nz)(-,-,);
  return xyz;
}

func get_isocrv(xyz)
{
  r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
  rp= r+!r;
  theta= acos(xyz(3,..)/rp);
  phi= atan(xyz(2,..),xyz(1,..)+(!r));
  y32= sin(theta)^2*cos(theta)*cos(2*phi);
  var= r*(1.+y32);
  return var;
}

func mak_y32(nx, ny, nz)
{
  if(is_void(nx) || nx < 2) nx= 40;
  if(is_void(ny) || ny < 2) ny= nx;
  if(is_void(nz) || nz < 2) nz= nx;
  xyz= array(0.0, 3, nx,ny,nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,ny)(-,);
  xyz(3,..)= span(-1,1,nz)(-,-,);
  r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
  rp= r+!r;
  theta= acos(xyz(3,..)/rp);
  phi= atan(xyz(2,..),xyz(1,..)+(!r));
  y32= sin(theta)^2*cos(theta)*cos(2*phi);
  f= r*(1.+y32);
  r= theta= phi= [];
  return f;
}

func get_y32(ngrid)
{
  extern f80, f100, f120, f150;
  extern have_ylm;

  // calling this function is faster than calling mak_y32
  // if the file already exists or sav_y32 has previously 
  // been called
  if(have_ylm) {
    if(ngrid == 80) var= f80;
    else if(ngrid == 100) var= f100;
    else if(ngrid == 120) var= f120;
    else if(ngrid == 150) var= f150;
    else var= mak_y32(ngrid);
  } else if( open("y32vals.pdb", "rb", 1) ) {
    /* the file exists, so restore function values */
    restore,openb("y32vals.pdb");
    have_ylm= 1;
    if(ngrid == 80) var= f80;
    else if(ngrid == 100) var= f100;
    else if(ngrid == 120) var= f120;
    else if(ngrid == 150) var= f150;
    else var= mak_y32(ngrid);
  } else {
    var= mak_y32(ngrid);
  }
  return var;
}

func sav_y32
{
  extern f80, f100, f120, f150;

  f80= mak_y32(80);
  f100= mak_y32(100);
  f120= mak_y32(120);
  f150= mak_y32(150);
  have_ylm= 1;
  f= createb("y32vals.pdb");
  save,f,f80,f100,f120,f150;
  close,f;
}

grab_32= 0;
if(grab_32) get_y32;

test3d_n= [20,20,20];

// create or activate the first OpenGL window
win3d,0;

thetalight= pi/4.0;
light3d, diffuse=.7, specular=1, sdir=[cos(thetalight),.25,sin(thetalight)];
mak_sphere,73,73;
stdview3d;
// sho_texiso,16;

/* show objects while rotating with the mouse if non-zero */
always_show_obj3d,0;

cage_style3d,0;

write,"mak_sphere,ntheta,nphi to make a new sphere";
write,"mak_slice to slice a pair of iso-surfaces and add a slicing plane";
write,"do_rot to rotate the existing scene";

if(is_void(bequiet)) {
  testgl;
  bequiet= 1;
}

