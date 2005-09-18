/*
 * $Id: testisotree.i,v 1.1 2005-09-18 22:07:55 dhmunro Exp $
 * Test the version of the iso-surface program that uses octrees.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

require, "dlist3d.i"
use_interleave= 0;

require, "pnm.i";
require, "contour.i";

palette3d, "earth.gp";

func sho_slicecrv(ngrid,point=,normal=)
{
  if(is_void(point)) point= [0.0, 0.0, 0.0];
  if(is_void(normal)) normal= [1.0, 0.0, 0.0];
  res= do_slicetreecrv(ngrid,point,normal);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    write,"time to build tree is ",res(2);
    write,"time to extract slicing plane is ",res(3);
    write,"time to set colors is ",res(4);
    write,"time to draw slicing plane is ",res(5);
    write,"time to build display list is ",res(6);
  }
}

func do_slicetreecrv(ngrid,point,normal)
{
  extern tree;
  extern gl_rr,gl_gg,gl_bb, gl_ncolr, gl_ctab;
  extern n_poly_3d, n_tri_3d;

  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ngrid;
    ny= ngrid;
    nz= ngrid;
  }
  xyz= array(0.0, 3, nx, ny, nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,ny)(-,);
  xyz(3,..)= span(-1,1,nz)(-,-,);
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

  palette3d,"stern.gp";
  vcen= f(zcen,zcen,zcen);  /* zone centered */
  varmax= max(vcen);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;
  clear3d;

  if(is_void(point)) point= [0.0, 0.0, 0.0];
  if(is_void(pl_normal)) pl_normal= [1.0, 0.0, 0.0];

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides fi guard is non-zero. 
     This means that the volume being contoured is a bit less 
     than -1 to 1.
  */
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  tree= mak_slice_treecrv(xyz,guard=0);
  timer,tfin;
  tree_tim= (tfin-tstart)(3);

  timer, tstart;
  tris= slice_treecrv(xyz,point,pl_normal,tree);
  timer,tfin;
  slice_tim= (tfin-tstart)(3);

  timer, tstart;
  triptr= &tris;
  while(1) {
    if(is_void(*triptr)) break;
    if(!triptr) break;
    ntri= triptr->numTri;
    if(ntri) {
      /* one color per triangle based on cellID */
      colr= map2color(ntri, *(triptr->cellIDs), vcen, vmin=0.0, vmax=varmax);
      triptr->colors= &colr;
    }
    triptr= triptr->next;
  } ;
  timer, tfin;
  colr_tim= (tfin-tstart)(3);

  timer,tstart;
  /* draw all tri arrays in the lists */
  if(tris) {
      res= pltrilists3d(tris,emit=1);
  }
  timer,tfin;
  list_tim= (tfin-tstart)(3);
  if(!is_void(res)) {
    numTri += res;
    nStrip += res;
    nVert += 3*res;
  }

  timer, tstart;
  draw3d, !making_movie;
  if(is_void(keepview)) {
    stdview3d;
  } else if(!keepview) {
    stdview3d;
  }
  timer, tfin;
  draw_tim= (tfin-tstart)(3);
  n_tri_3d= numTri;
  n_poly_3d= numTri;
  sizes= [nx, ny, nz];
//  chek_slice, point, normal, sizes, tris;
  write,"Grid is ", nx, "by", ny, "by", nz;
  return [numTri, tree_tim, slice_tim, colr_tim, draw_tim, list_tim];
}

func slicecrv_movie(ngrid)
{
  extern tree;
  extern gl_rr,gl_gg,gl_bb, gl_ncolr, gl_ctab;

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

  palette3d,"stern.gp";
  vcen= f(zcen,zcen,zcen);
  varmax= max(vcen);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  slice_tim= list_tim= draw_tim= colr_tim= 0.0;
  lookat3d,2.9*normal,[0.0,0.0,0.0],[0.0,0.0,1.0];

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

  tstart= tfin= array(0.0, 3);
  timer, tstart;
  tree= mak_slice_treecrv(xyz,guard=0);
  timer,tfin;
  tree_tim= (tfin-tstart)(3);

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
    triptr= &tris;
    while(1) {
      if(is_void(*triptr)) break;
      if(!triptr) break;
      ntri= triptr->numTri;
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

  for(jj= 1; jj <= npts; jj++) {
    timer, tstart;
    tris= slice_treecrv(xyz,point0,norms(,jj),tree);
    timer,tfin;
    slice_tim += (tfin-tstart)(3);

    timer, tstart;
    triptr= &tris;
    while(1) {
      if(is_void(*triptr)) break;
      if(!triptr) break;
      ntri= triptr->numTri;
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

  for(jj= 1; jj <= npts; jj++) {
    timer, tstart;
    point= point0+posns2(jj)*normal;
    tris= slice_treecrv(xyz,point,normal,tree);
    timer,tfin;
    slice_tim += (tfin-tstart)(3);

    timer, tstart;
    triptr= &tris;
    while(1) {
      if(is_void(*triptr)) break;
      if(!triptr) break;
      ntri= triptr->numTri;
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

  write,"Slice planes computed using an Octree";
  write,"Grid is ", nx, "by", ny, "by", nz;
  write,"positions range from ",xyzmin," to ",xyzmax;
  write,"time to build octree is ",tree_tim;
  write,"time to extract extract plane is ",slice_tim;
  write,"time to set colors is ", colr_tim;
  write,"time to draw slicing plane is ", draw_tim;
  write,"time to build display list is ",list_tim;
  write,"frames per second is",(3*npts)/(slice_tim+colr_tim+draw_tim+list_tim);
}

func sho_slice(ngrid,point=,normal=)
{
  if(is_void(point)) point= [0.0, 0.0, 0.0];
  if(is_void(normal)) normal= [1.0, 0.0, 0.0];
  res= do_slicetree(ngrid,point,normal);
  if(!is_void(res)) {
    write,"number of triangles is ",long(res(1));
    write,"time to extract slicing plane is ",res(2);
    write,"time to draw slicing plane is ",res(3);
    write,"time to build display list is ",res(4);
  }
}

func do_slicetree(ngrid,point,normal)
{
  extern gl_rr,gl_gg,gl_bb, gl_ncolr, gl_ctab;
  extern tree;
  extern n_poly_3d, n_tri_3d;
  extern have_slice;
  extern have_sphere;

  have_sphere= 0;
  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ngrid;
    ny= ngrid;
    nz= ngrid;
  }
  xyz= array(0.0, 3, nx,ny,nz);
  xyz(1,..)= span(-1,1,nx);
  xyz(2,..)= span(-1,1,nx)(-,);
  xyz(3,..)= span(-1,1,nx)(-,-,);
  r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
  theta= acos(xyz(3,..)/(r+!r));
  phi= atan(xyz(2,..),xyz(1,..)+(!r));
  y32= sin(theta)^2*cos(theta)*cos(2*phi);
  f= r*(1.+y32);
  r= theta= phi= [];
  sizes= (dimsof(f))(2:4);

  palette3d,"stern.gp";
  vcen= f(zcen,zcen,zcen);
  varmax= max(vcen);
  varfac= gl_ncolr/varmax;

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  numTri= nStrip= nVert= 0;
  clear3d;

  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];
  if(is_void(point)) point= [0.0, 0.0, 0.0];
  if(is_void(normal)) normal= [1.0, 0.0, 0.0];

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides. This means that the volume
     being contoured is a bit less than -1 to 1.
  */
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  tris= slice_tree(sizes,origin,delta,point,normal);
  timer,tfin;
  slice_tim= (tfin-tstart)(3);

  triptr= &tris;
  while(1) {
    if(is_void(*triptr)) break;
    if(!triptr) break;
    ntri= triptr->numTri;
    if(ntri) {
      /* one color per triangle based on cellID */
      cellids= (*(triptr->cellIDs))(1:ntri);
      /* cellIds are C-style so add one before using as yorick indices */
      vals= vcen(*)(cellids+1);
      nd= long(vals*varfac)+1;
      nd= min(nd, gl_ncolr);
      colr= array(0.0, 3, ntri);
      colr(1,)= gl_rr(nd)/256.0;
      colr(2,)= gl_gg(nd)/256.0;
      colr(3,)= gl_bb(nd)/256.0;
      triptr->colors= &colr;
    }
    triptr= triptr->next;
  } ;

  timer,tstart;
  /* draw all tri arrays in the lists */
  if(tris) {
    res= pltrilists3d(tris,emit=1);
  }
  timer,tfin;
  list_tim= (tfin-tstart)(3);
  if(!is_void(res)) {
    numTri += res;
    nStrip += res;
    nVert += 3*res;
  }

  timer, tstart;
  draw3d, !making_movie;
  if(is_void(keepview)) {
    stdview3d;
  } else if(!keepview) {
    stdview3d;
  }
  timer, tfin;
  draw_tim= (tfin-tstart)(3);
  n_tri_3d= numTri;
  n_poly_3d= numTri;
  sizes= [nx, ny, nz];
//  chek_slice, point, normal, sizes, tris;
  write,"Grid is ", nx, "by", ny, "by", nz;
  return [numTri, slice_tim, draw_tim, list_tim];
}

func slice_movie(ngrid)
{
  extern tree;
  extern gl_rr,gl_gg,gl_bb, gl_ncolr, gl_ctab;

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

  palette3d,"stern.gp";
  vcen= f(zcen,zcen,zcen);
  varmax= max(vcen);

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  slice_tim= list_tim= draw_tim= colr_tim= 0.0;
  lookat3d,2.9*normal,[0.0,0.0,0.0],[0.0,0.0,1.0];

  /* start the slice plane at the back, sweep through
     to the middle, rotate it 180 degrees, and sweep
     to the front. */
  posns1= span(xyzmin, 0, npts);
  posns2= span(0, xyzmax, npts);
  norms= array(0.0, 3, npts);
  norms(1,)= cos(span(0.0, pi, npts));
  norms(2,)= sin(span(0.0, pi, npts));
  norms(3,)= 0.0;
  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(ny-1.0), 2.0/(nz-1.0)];

  write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

  tstart= tfin= array(0.0, 3);

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides. This means that the volume
     being contoured is a bit less than -1 to 1.
  */
  for(jj= 1; jj <= npts; jj++) {
    timer, tstart;
    point= point0+posns1(jj)*normal;
    tris= slice_tree(sizes,origin,delta,point,normal);
    timer,tfin;
    slice_tim += (tfin-tstart)(3);

    timer, tstart;
    triptr= &tris;
    while(1) {
      if(is_void(*triptr)) break;
      if(!triptr) break;
      ntri= triptr->numTri;
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

  for(jj= 1; jj <= npts; jj++) {
    timer, tstart;
    tris= slice_tree(sizes,origin,delta,point0,norms(,jj));
    timer,tfin;
    slice_tim += (tfin-tstart)(3);

    timer, tstart;
    triptr= &tris;
    while(1) {
      if(is_void(*triptr)) break;
      if(!triptr) break;
      ntri= triptr->numTri;
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

  for(jj= 1; jj <= npts; jj++) {
    timer, tstart;
    point= point0+posns2(jj)*normal;
    tris= slice_tree(sizes,origin,delta,point,normal);
    timer,tfin;
    slice_tim += (tfin-tstart)(3);

    timer, tstart;
    triptr= &tris;
    while(1) {
      if(is_void(*triptr)) break;
      if(!triptr) break;
      ntri= triptr->numTri;
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

  write,"Slice planes computed using an Octree";
  write,"Grid is ", nx, "by", ny, "by", nz;
  write,"positions range from ",xyzmin," to ",xyzmax;
  write,"time to extract extract plane is ",slice_tim;
  write,"time to set colors is ", colr_tim;
  write,"time to draw slicing plane is ", draw_tim;
  write,"time to build display list is ",list_tim;
  write,"frames per second is",(3*npts)/(slice_tim+colr_tim+draw_tim+list_tim);
}

func print_tree(tree)
{
  trptr= &tree;
  while(trptr && *trptr) {
    write,"maxdepth is",trptr->maxdepth;
    write,"start is";*(trptr->start);
    write,"chunk is";*(trptr->chunk);
    write,"triz is";*(trptr->trsiz);
    write,"ranges are";*(trptr->ranges);
    trptr= trptr->next;
  }
}

func chek_iso(level, sizes, tris)
{
  /* The input is an iso-surface level, the dimensions of the
     data array, and an array of triangles.
     Evaluate the function at every triangle vertex to see
     if the iso-surface is correct */

  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);
  triptr= &tris;
  while(triptr && *triptr) {
    xyz= *(triptr->xyzverts);
    r= abs(xyz(1,..), xyz(2,..), xyz(3,..));
    theta= acos(xyz(3,..)/r);
    phi= atan(xyz(2,..), xyz(1,..)+(!r));
    y32= sin(theta)^2*cos(theta)*cos(2*phi);
    f= r*(1.+y32);
    err= f-level;
    rerr= abs(err)/(abs(f)+(!f));
    write,"minimum of function minus iso-level is",min(err);
    write,"maximum of function minus iso-level is",max(err);
    write,"maximum relative error is",max(rerr);
    write,"RMS fractional error is", rerr(*)(rms);
    triptr= triptr->next;
  }
}

func chek_slice(point, normal, sizes, tris)
{
  /* The input is a slicing plane, the dimensions of the
     data array, and an array of triangles.
     Evaluate the function at every triangle vertex to see
     if it is on the plane. */

  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);
  triptr= &tris;
  while(triptr && *triptr) {
    xyz= *(triptr->xyzverts);
    dist= ((xyz-point)*normal)(sum,..);
    write,"minimum distance to plane is",min(dist);
    write,"maximum distance to plane is",max(dist);
    triptr= triptr->next;
  }
}

func iso_movie(ngrid)
{
  extern tree;

  if(is_void(ngrid) || ngrid < 6) {
    nx= test3d_n(1);
    ny= test3d_n(2);
    nz= test3d_n(3);
  } else {
    nx= ngrid;
    ny= ngrid;
    nz= ngrid;
  }
  nlev= 20;
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
  f= r*(1.+y32);
  r= theta= phi= [];

  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(nx-1.0), 2.0/(nx-1.0)];
  colr= [0.0,1.0,0.2];
  levels= span(iso_lo, iso_hi, nlev);
  tstart= tfin= array(0.0, 3);
  iso_tim= list_tim= draw_tim= 0.0;
  lookat3d,[2.2,-0.03,-0.22],[0.0,-0.03,-0.22],[0.0,0.0,1.0];

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides. This means that the volume
     being contoured is a bit less than -1 to 1.
  */
  for(jj= 1; jj <= nlev; jj++) {
    timer, tstart;
    tris= iso3cenreg(origin,delta,f,levels(jj),colr);
    timer,tfin;
    iso_tim += (tfin-tstart)(3);
    /* draw all tri arrays in the lists */
    clear3d;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer, tstart;
    list_tim += (tstart-tfin)(3);
    timer, tstart;
    draw3d, !making_movie;
    timer, tfin;
    draw_tim += (tfin-tstart)(3);
  }
  write,"Iso-surfaces computed using a Marching Cubes-like scheme";
  write,"Grid is ", nx, "by", ny, "by", nz;
  write,"iso-levels range from ",iso_lo," to ",iso_hi;
  write,"time to extract iso-surface is ",iso_tim;
  write,"time to draw iso-surface is ", draw_tim;
  write,"time to build display list is ",list_tim;
  write,"frames per second is",nlev/(iso_tim+draw_tim+list_tim);
}

func isotree_movie(ngrid)
{
  extern tree;

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
  f= r*(1.+y32);
  r= theta= phi= [];

  origin= [-1.0,-1.0,-1.0];
  delta= [2.0/(nx-1.0), 2.0/(nx-1.0), 2.0/(nx-1.0)];
  colr= [0.0,1.0,0.2];
  levels= span(iso_lo, iso_hi, nlev);
  tstart= tfin= array(0.0, 3);
  timer, tstart;
  tree= mak_isotree(f);
  timer,tfin;
  tree_tim= (tfin-tstart)(3);
  iso_tim= list_tim= draw_tim= 0.0;
  lookat3d,[2.2,-0.03,-0.22],[0.0,-0.03,-0.22],[0.0,0.0,1.0];

  /* NOTE: the input array is assumed to include a layer of
     guard cells on all sides. This means that the volume
     being contoured is a bit less than -1 to 1.
  */
  for(jj= 1; jj <= nlev; jj++) {
    timer, tstart;
    tris= iso3_tree(origin,delta,f,levels(jj),colr,tree);
    timer,tfin;
    iso_tim += (tfin-tstart)(3);
    /* draw all tri arrays in the lists */
    clear3d;
    if(tris) {
      res= pltrilists3d(tris);
    }
    timer, tstart;
    list_tim += (tstart-tfin)(3);
    timer, tstart;
    draw3d, !making_movie;
    timer, tfin;
    draw_tim += (tfin-tstart)(3);
  }
  write,"Iso-surfaces computed an Octree";
  write,"Grid is ", nx, "by", ny, "by", nz;
  write,"iso-levels range from ",iso_lo," to ",iso_hi;
  write,"number of triangles is ",long(numTri);
  write,"time to build octtree is ",tree_tim;
  write,"time to extract iso-surface is ",iso_tim;
  write,"time to draw iso-surface is ", draw_tim;
  write,"time to build display list is ",list_tim;
  write,"frames per second is",nlev/(iso_tim+draw_tim+list_tim);
}

#if 0
test3d_n= [20,20,20];

/* show objects while rotating with the mouse if non-zero */
always_show_obj3d,0;

// create or activate the first OpenGL window
win3d,0;

thetalight= pi/4.0;
light3d, diffuse=.7, specular=1, sdir=[cos(thetalight),.25,sin(thetalight)];
isotree_movie,30;
stdview3d;
#endif

