/*
 * $Id: glprofile.i,v 1.1 2005-09-18 22:07:50 dhmunro Exp $
 * run this file to generate profiling information for the OpenGL
 * based 3D graphics package in yorick
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

require, "dlist3d.i"
use_interleave= 0;
require, "contour.i";

palette3d, "earth.gp";

// gl_win_size,1024,820;  /* width,height */

n_pass= 5;

func mak_slicecrv(ngrid)
{
  extern tree;
  extern gl_rr,gl_gg,gl_bb, gl_ncolr, gl_ctab;
  extern tot_slice_tri;

  normal= [1.0, 0.0, 0.0];
  point0= [0.0, 0.0, 0.0];
  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ngrid;
    ny= ngrid;
    nz= ngrid;
  }
  npts= 20;

  xyz= array(0.0, 3, nx,ny,nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,nx)(-,);
  xyz(3,..)= span(-1,1,nx)(-,-,);
  /* Now rotate each plane by phi degrees around the z-axis
     relative to the one below it. */
  phi= 0.5*pi/(nz-1.0);
  snph= sin(phi);
  csph= cos(phi);
  for(i= 2; i <= nz; i++) {
    xyz(1,,,i)= csph*xyz(1,,,i-1)-snph*xyz(2,,,i-1);
    xyz(2,,,i)= snph*xyz(1,,,i-1)+csph*xyz(2,,,i-1);
  }
  r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
  theta= acos(xyz(3,..)/(r+!r));
  phi= atan(xyz(2,..),xyz(1,..)+(!r));
  y32= sin(theta)^2*cos(theta)*cos(2*phi);
  f= r*(1.+y32);   /* point centered */
  r= theta= phi= [];
  xyzmin= min(xyz);
  xyzmax= max(xyz)
  sizes= (dimsof(f))(2:4);

  vcen= f(zcen,zcen,zcen);
  varmax= max(vcen);
  varfac= gl_ncolr/varmax;

  /* start the slice plane at the back, sweep through
     to the middle, rotate it 180 degrees, and sweep
     to the front. */
  posns1= span(xyzmin, 0, npts);
  posns2= span(0, xyzmax, npts);
  norms= array(0.0, 3, npts);
  norms(1,)= cos(span(0.0, pi, npts));
  norms(2,)= sin(span(0.0, pi, npts));
  norms(3,)= 0.0;

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  slice_tim= colr_tim= 0.0;
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  tree= mak_slice_treecrv(xyz,guard=0);
  timer,tfin;
  tree_tim= (tfin-tstart)(3);

  tris_sav= array(pointer, 3*npts);
  tris_sav(*)= &[];
  tot_tri= 0;

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides. This means that the volume
     being contoured is a bit less than -1 to 1.
  */
  for(jj= 1; jj <= npts; jj++) {
    timer, tstart;
    point= point0+posns1(jj)*normal;
    tris= slice_treecrv(xyz,point,normal,tree);
    timer,tfin;
    slice_tim += (tfin-tstart)(3);

    timer, tstart;
    if(tris) {
      if(!tris_sav(jj)) {
        tris_sav(jj)= &tris;
      } else {
        tris.next= tris_sav(jj);
        tris_sav(jj)= &tris;
      }

      triptr= tris_sav(jj);
      while(1) {
        if(!triptr) break;
        if(is_void(*triptr)) break;
        ntri= triptr->numTri;
        tot_tri += ntri;
        if(ntri) {
          /* one color per triangle based on cellID */
          colr= map2color(ntri, *(triptr->cellIDs), vcen, vmin=0.0, 
		                vmax=varmax);
          triptr->colors= &colr;
        }
        triptr= triptr->next;
      } ;
      timer, tfin;
      colr_tim += (tfin-tstart)(3);
    }
  }

  for(jj= 1; jj <= npts; jj++) {
    timer, tstart;
    tris= slice_treecrv(xyz,point0,norms(,jj),tree);
    timer,tfin;
    slice_tim += (tfin-tstart)(3);

    jp= jj+npts;
    timer, tstart;
    if(tris) {
      if(!tris_sav(jp)) {
        tris_sav(jp)= &tris;
      } else {
        tris.next= tris_sav(jp);
        tris_sav(jp)= &tris;
      }

      triptr= tris_sav(jp);
      while(1) {
        if(!triptr) break;
        if(is_void(*triptr)) break;
        ntri= triptr->numTri;
        tot_tri += ntri;
        if(ntri) {
          /* one color per triangle based on cellID */
          colr= map2color(ntri, *(triptr->cellIDs), vcen, vmin=0.0, 
		                vmax=varmax);
          triptr->colors= &colr;
        }
        triptr= triptr->next;
      } ;
      timer, tfin;
      colr_tim += (tfin-tstart)(3);
    }
  }

  for(jj= 1; jj <= npts; jj++) {
    timer, tstart;
    point= point0+posns2(jj)*normal;
    tris= slice_treecrv(xyz,point,normal,tree);
    timer,tfin;
    slice_tim += (tfin-tstart)(3);

    timer, tstart;
    jp= jj+2*npts;
    if(tris) {
      if(!tris_sav(jp)) {
        tris_sav(jp)= &tris;
      } else {
        tris.next= tris_sav(jp);
        tris_sav(jp)= &tris;
      }

      triptr= tris_sav(jp);
      while(1) {
        if(!triptr) break;
        if(is_void(*triptr)) break;
        ntri= triptr->numTri;
        tot_tri += ntri;
        if(ntri) {
          /* one color per triangle based on cellID */
          colr= map2color(ntri, *(triptr->cellIDs), vcen, vmin=0.0, 
		                vmax=varmax);
          triptr->colors= &colr;
        }
        triptr= triptr->next;
      } ;
      timer, tfin;
      colr_tim += (tfin-tstart)(3);
    }
  }

  write,"Slice planes computed using an Octree";
  write,"Grid is ", nx, "by", ny, "by", nz;
  write,"positions range from ",xyzmin," to ",xyzmax;
  write,"time to build octree is ",tree_tim;
  write,"time to extract extract plane is ",slice_tim;
  write,"time to compute colors is ",colr_tim;
  write,"number of triangles for all frames is ",tot_tri;
  write,"Storage for results is roughly ",8*(3*3+3*3+3)*tot_tri," bytes";

  tot_slice_tri= tot_tri;
  return tris_sav;
}

func slicecrv_movie(tris_sav)
{
  extern tot_slice_tri, n_pass;

  palette3d,"stern.gp";

  normal= [1.0, 0.0, 0.0];
  lookat3d,2.9*normal,[0.0,0.0,0.0],[0.0,0.0,1.0];

  list_tim= draw_tim= 0.0;
  tstart= tfin= array(0.0, 3);

  /* Draw all slice planes */
  nplanes= numberof(tris_sav);
  for(ii= 1; ii <= n_pass; ii++) {
    for(jj= 1; jj <= nplanes; jj++) {
      triptr= tris_sav(jj);
      tris= *triptr;

      timer,tstart;
      clear3d;
      /* draw all tri arrays in the lists */
      if(tris) {
        res= pltrilists3d(tris,emit=1);
      }
      timer,tfin;
      list_tim += (tfin-tstart)(3);

      timer, tstart;
      draw3d, !making_movie;
      timer, tfin;
      draw_tim += (tfin-tstart)(3);
    }
  }

  write,"number of passes is ",n_pass;
  write,"time to draw slicing plane is ", draw_tim;
  write,"time to build display list is ",list_tim;
  write,"frames per second is",n_pass*nplanes/(draw_tim+list_tim);
  write,"triangles per second is ",n_pass*tot_slice_tri/(draw_tim+list_tim);
  write,"triangle rendering rate is ",n_pass*tot_slice_tri/draw_tim;
}

func mak_isos(ngrid)
{
  extern tot_iso_tri, do_collapse;

  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ngrid;
    ny= ngrid;
    nz= ngrid;
  }
  nlev= 40;
  iso_lo= 0.1;
  iso_hi= 1.0;

  xyz= array(0.0, 3, nx,ny,nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,nx)(-,);
  xyz(3,..)= span(-1,1,nx)(-,-,);
  r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
  theta= acos(xyz(3,..)/(r+!r));
  phi= atan(xyz(2,..),xyz(1,..)+(!r));
  y32= sin(theta)^2*cos(theta)*cos(2*phi);
  f= r*(1.+y32);   /* point centered */
  r= theta= phi= [];

  origin= [-1.0, -1.0, -1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  colr= [0.0,1.0,0.2];
  levels= span(iso_lo, iso_hi, nlev);

  tstart= tfin= array(0.0, 3);
  timer, tstart;
  tree= mak_isotree(f);
  timer,tfin;
  tree_tim= (tfin-tstart)(3);
  iso_tim= append_tim= 0.0;

  tris_sav= array(pointer, nlev);
  tris_sav(*)= &[];
  tot_tri= 0;

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides. This means that the volume
     being contoured is a bit less than -1 to 1.
  */
  for(jj= 1; jj <= nlev; jj++) {
    timer, tstart;
    tris= iso3_tree(origin,delta,f,levels(jj),colr,tree);
    timer,tfin;
    iso_tim += (tfin-tstart)(3);

    timer, tstart;
    if(tris) {
      if(!tris_sav(jj)) {
        tris_sav(jj)= &tris;
      } else {
        tris.next= tris_sav(jj);
        tris_sav(jj)= &tris;
      }

      triptr= tris_sav(jj);
      while(1) {
        if(!triptr) break;
        if(is_void(*triptr)) break;
        ntri= triptr->numTri;
        tot_tri += ntri;
        triptr= triptr->next;
      } ;
    }
    if(do_collapse) {
      /* collapse into a single triangle array if requested */
      newtris= &tris;
      num= SizeTriArrays3d(tris_sav(jj));
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
      if(tris_sav(jj)->var2) newtris.var2= &array(0.0, 3, num);
      CollapseTriArrays3d, -3, tris_sav(jj), &newtris;
      tris_sav(jj)= &newtris;
    }
    timer, tfin;
    append_tim += (tfin-tstart)(3);
  }

  write,"Iso-surfaces computed using an Octree";
  write,"Grid is ", nx, "by", ny, "by", nz;
  write,"number of iso-levels is ",nlev;
  write,"iso-levels range from ",iso_lo," to ",iso_hi;
  write,"time to append is ",append_tim;
  write,"time to extract iso-surfaces is ",iso_tim;
  write,"number of triangles for all frames is ",tot_tri;
  write,"Storage for results is roughly ",8*(3*3+3*3+3)*tot_tri," bytes";

  tot_iso_tri= tot_tri;
  return tris_sav;
}

func isotree_movie(tris_sav)
{
  extern tot_iso_tri, n_pass;

  nlev= numberof(tris_sav);

  tstart= tfin= array(0.0, 3);
  list_tim= draw_tim= clr_tim= 0.0;
  lookat3d,[2.2,-0.03,-0.22],[0.0,-0.03,-0.22],[0.0,0.0,1.0];

  /* Draw all iso-surfaces */
  for(ii= 1; ii <= n_pass; ii++) {
    for(jj= 1; jj <= nlev; jj++) {
      triptr= tris_sav(jj);
      tris= *triptr;

      timer,tstart;
      clear3d;
      timer,tfin;
      clr_tim += (tfin-tstart)(3);
      timer,tstart;
      /* draw all tri arrays in the lists */
      if(tris) {
        res= pltrilists3d(tris);
      }
      timer,tfin;
      list_tim += (tfin-tstart)(3);

      timer, tstart;
      draw3d, !making_movie;
      timer, tfin;
      draw_tim += (tfin-tstart)(3);
    }
  }

  write,"number of passes is ",n_pass;
  write,"time to draw slicing plane is ", draw_tim;
  write,"time to build display list is ",list_tim;
  write,"time to clear screen is ",clr_tim;
  write,"frames per second is",n_pass*nlev/(draw_tim+list_tim);
  write,"triangles per second is ",n_pass*tot_iso_tri/(draw_tim+list_tim);
  write,"triangle rendering rate is ",n_pass*tot_iso_tri/draw_tim;
}

func runprof(size)
{
  extern make_strip;

  if(is_void(size)) size= 32;
  make_strip= 0;
  // rotate each scene the same number of times
  num_rotn= 12;

  "sho_isoreg: iso-surface on regular grid";
  sho_isoreg,size
  do_rot,num_rotn;
  "sho_isoregzcen: iso-surface on regular grid with zone-centered data";
  sho_isoregzcen,size
  do_rot,num_rotn;
  "sho_isoregndx: iso-surface on regular grid returning indices into point list";
  sho_isoregndx,size
  do_rot,num_rotn;
//  "sho_isoregzcenndx: iso-surface on regular grid with zone-centered data returning indices into point list";
//  sho_isoregzcenndx,size
//  do_rot,num_rotn;
  "sho_isocrv: iso-surface on curvilinear grid";
  sho_isocrv,size
  do_rot,num_rotn;
//  "sho_isocrvzcen: iso-surface on curvilinear grid with zone-centered data";
//  sho_isocrvzcen,size
//  do_rot,num_rotn;
//  "sho_isocrvndx: iso-surface on curvilinear grid returning indices into point list";
//  sho_isocrvndx,size
//  do_rot,num_rotn;
//  "sho_isoregzcrvndx: iso-surface on curvilinear grid returning indices into point list";
//  sho_isoregzcrvndx,size
//  do_rot,num_rotn;
  "sho_isohex: iso-surface on a pile of hexahedra";
  sho_isohex,size
  do_rot,num_rotn;
  "sho_iso: iso-surface on a regular grid";
  sho_iso,size
  do_rot,num_rotn;
//  "sho_isondx: iso-surface on a regular grid returning indices into point list";
//  sho_isondx,size
//  do_rot,num_rotn;

  "sho_isotree: iso-surface on regular grid using an octtree";
  sho_isotree,size;
  do_rot,num_rotn;
  "sho_drum: oscillating drum head";
  drum_size= size;
  sho_drum;
  "mak_slice: slice plane and iso-surface";
  mak_slice,size;
  do_rot,num_rotn;
  "mak_sphere: sphere made with quad strips";
  mak_sphere,4*size,4*size;
  do_rot,num_rotn;
  "rlines: polylines";
  rlines,12,1024*size;
  "slicecrv_movie: slicing plane on a curvilinear grid using an octtree";
  slicecrv_movie,size;
  "slice_movie: slicing plane on a regular grid using an octtree";
  slice_movie,size;
  "isotree_movie: iso-surface on a regular grid using an octtree";
  isotree_movie,size;
}

// create the OpenGL window
// make_3d;  // sometimes it dramatically reduces performance to call this function!!

system,"netstat -s |grep etra | grep -v term | grep -v TCP"
/* runprof,32; */
if(is_void(ngrid)) ngrid= 60;
tri_arr= mak_slicecrv(ngrid);
win3d,0;
slicecrv_movie,tri_arr;
tri_arr= [];
do_collapse= 1;
tris_sav= mak_isos(ngrid);
isotree_movie,tris_sav;
system,"netstat -s |grep etra | grep -v term | grep -v TCP"

// quit;
