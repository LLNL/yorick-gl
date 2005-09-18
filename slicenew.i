/*
 * $Id: slicenew.i,v 1.1 2005-09-18 22:07:54 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

/* ------------------------------------------------------------------------ */

func nuslice2x(plane, tris)
/* DOCUMENT slice2, plane, tris

     Slice a polygon list, retaining only those triangles or
     parts of triangles on the positive side of PLANE, that is,
     the side where xyz(+)*PLANE(+:1:3)-PLANE(4) > 0.0.

   SEE ALSO: slice2, slice2_precision
 */
{
  _slice2x= 1;
  return nuslice2(plane, tris);
}

func nuslice2(plane, tris)
/* DOCUMENT slice2, plane, tris

     Slice an array of triangles, retaining only those triangles or
     parts of triangles on the positive side of PLANE, that is,
     the side where xyz(+)*PLANE(+:1:3)-PLANE(4) > 0.0.

     In order to plot two intersecting slices, one could
     slice (for example) the horizontal plane twice (slice2x) -
     first with the plane of the vertical slice, then with minus
     that same plane.  Then, plot first the back part of the
     slice, then the vertical slice, then the front part of the
     horizontal slice.  Of course, the vertical plane could
     be the one to be sliced, and "back" and "front" vary
     depending on the view point, but the general idea always
     works.

   SEE ALSO: slice2x, slice2_precision
 */
{
  have_values= !is_void(values);

  ntri= tris.numTris;
  /* make a list of the triangle number for each points in xyzverts */
  ndxs= (long(span(1, ntri, ntri))(-:1:3,))(*);

  /* form dot products of all the points with the given plane */
  xyzverts= *(tris.xyzverts);
  normals= *(tris.normals);
  ids= *(tris.cellIDs);
  colr= *(tris.colors);
  var2= *(tris.var2);
  hasVar2= !is_void(var2);
  dp= xyzverts(+,)*plane(+:1:3) - plane(4);

  /* separate into lists of unclipped and partially clipped polys */
  if (!slice2_precision) {
    /* if precision is not set, slice exactly at dp==0.0, with
     * any points exactly at dp==0.0 treated as if they had dp>0.0 */
    keep= (dp>=0.0);
  } else {
    /* if precision is set, polygons are clipped to +-precision,
     * so that any poly crossing +precision is clipped to dp>=+precision,
     * any poly crossing -precision is clipped to dp<=-precision, and
     * any poly lying entirely between +-precision is discarded entirely */
    keep= (dp>=slice2_precision);
  }
  nkeep= long(histogram(ndxs, keep));
  mask0= (nkeep==nverts);  /* all vertices on "right" side of plane */
  mask1= (nkeep!=0 & !mask0);  /* part of the vertices on the "right" side */
  list1= where(mask1);
  if (numberof(list1)) {
    normc= normals(,,list1);
    colrc= colr(,list1);
    if(hasVar2) var2c= var2(,list1);
    idsc= ids(list1);
    list= where(mask1(ndxs));
    xyzc= xyzverts(,,list);
  }
  if (_slice2x) {
    if (!slice2_precision) {
      mask2= !nkeep;
      normc0= normc;
      colrc0= colrc;
      if(hasVar2) var2c0= var2c;
      idsc0= idsc;
      xyzc0= xyzc;
    } else {
      keep2= (dp>-slice2_precision);
      nkeep2= long(histogram(ndxs, keep2));
      mask2= (nkeep!=0 & nkeep<nverts);
      list2= where(mask2);
      if (numberof(list2)) {
        normc0= normals(,,list2);
        colrc0= colr(,list2);
        if(hasVar2) var2c0= var2(,list2);
        idsc0= ids(list2);
        listc= where(mask2(ndxs));
        xyzc0= xyzverts(,,listc);
      }
      mask2= !nkeep2;
    }
    list2= where(mask2);
    if (numberof(list2)) {
      normcb= normals(,,list2);
      colrcb= colr(,list2);
      if(hasVar2) var2cb= var2(,list2);
      idscb= ids(list2);
      listcb= where(mask2(ndxs));
      xyzcb= xyzverts(,,listcb);
    } else {
      normcb= colrcb= var2cb= idscb= xyzcb= [];
    }
  }
  list0= where(mask0);
  if (numberof(list0)<numberof(ntri)) {
    xyznew= xyzverts(,,list0);
  }

  /* done if no partially clipped polys */
  if (!numberof(list1) && !numberof(listc)) return;
  if (!numberof(list1)) goto skip;

  /* get dot products and keep list for the clipped polys */
  dp= dp(list);
  if (slice2_precision) dp-= slice2_precision;
  keep= (dp>=0.0);

  /* get the indices of the first and last points in each clipped poly */
  last= nvertc(psum);
  frst= last - nvertc + 1;

  /* get indices of previous and next points in cyclic order */
  prev= next= indgen(0:numberof(keep)-1);
  prev(frst)= last;
  next(prev)= indgen(numberof(keep));

  _nuslice2_part;

  grow, nverts, nvertc;
  if (have_values) grow, values, valuec;
  grow, xyzverts, xyzc;

  if (_slice2x) {
  skip:
    nvertc= nvertc0;
    valuec= valuec0;
    xyzc= xyzc0;
    if (!slice2_precision) {
      keep= !keep;
    } else {
      dp= dp(listc)+slice2_precision;
      keep= (dp>=0.0);
    }

    _nuslice2_part;

    grow, nvertb, nvertc;
    if (have_values) grow, valueb, valuec;
    grow, xyzvertb, xyzc;
  }
}

func slice2only(plane, tris)
/* DOCUMENT slice2, plane, tris

     Slice an array of triangles, retaining only those triangles or
     parts of triangles on the positive side of PLANE, that is,
     the side where xyz(+)*PLANE(+:1:3)-PLANE(4) > 0.0.

     In order to plot two intersecting slices, one could
     slice (for example) the horizontal plane twice (slice2x) -
     first with the plane of the vertical slice, then with minus
     that same plane.

   SEE ALSO: slice2x, slice2_precision
 */
{
  if(!tris) return [];
  ntri= tris.numTri;
  /* make a list of the triangle number for each point in xyzverts */
  ndxs= (long(span(1, ntri, ntri))(-:1:3,))(*);

  /* form dot products of all the points with the given plane */
  xyzverts= *(tris.xyzverts);
  normals= *(tris.normals);
  ids= *(tris.cellIDs);
  colr= *(tris.colors);
  var2= *(tris.var2);
  dp= xyzverts(+,)*plane(+:1:3) - plane(4);

  /* separate into lists of unclipped and partially clipped polys */
  if (!slice2_precision) {
    /* if precision is not set, slice exactly at dp==0.0, with
     * any points exactly at dp==0.0 treated as if they had dp>0.0 */
    keep= (dp>=0.0);
  } else {
    /* if precision is set, polygons are clipped to +-precision,
     * so that any poly crossing +precision is clipped to dp>=+precision,
     * any poly crossing -precision is clipped to dp<=-precision, and
     * any poly lying entirely between +-precision is discarded entirely */
    keep= (dp>=slice2_precision);
  }
  nkeep= long(histogram(ndxs, keep));
  /* Create storage for the output triangles.
     There will be one output triangle for each entry in list0.
     There will be one triangle for everywhere that nkeep==1
     (one vertex is "to the right" of the plane) and two
     triangles where nkeep==2 (two vertices to the "right").
  */
  nulvar= [];
  numnew= sum(nkeep)-2*sum( nkeep==3 );
  if(numnew <= 0) return [];
  newxyz= array(0.0, 3, 3, numnew);
  newnorm= array(0.0, 3, 3, numnew);
  newcolr= array(float, 3, numnew);
  hasVar2= !is_void(var2);
  if(hasVar2) {
    newvar2= array(double, 3, numnew);
  } else {
    newvar2= nulvar;
  }
  newids= array(long, numnew);
  newtri= TriArrayGrp(numTri= numnew,
        xyzverts= &newxyz,
        normals= &newnorm,
        colors= &newcolr,
        cellIDs= &newids,
        var2= &newvar2,
        /* the next 3 are not used except while forming strips */
        nTris= &nulvar,  
        triEdg= &nulvar,
        triStart= &nulvar,
        next= &nulvar);
  /* call the routine that clips triangles and copies
     unaltered triangles */
  SliceTris3d,keep,nkeep,dp,&tris,&newtri;
  return newtri;
}

func test_slice2(offset)
{
    require,"testgl.i";

    extern n_poly_3d, n_tri_3d;
    extern have_slice;
    extern have_sphere;

    if(is_void(offset)) offset= 0.2;
    plane= [1.0,0.0,0.0,offset];
    write,format="slice plane is in direction (%f,%f,%f) with offset %f\n",
          plane(1),plane(2),plane(3),plane(4);

    ngrid= 40;
    nx= ny= nz= ngrid;
    xyz= array(0.0, 3, nx,ny,nz);
    xyz(1,..)= span(-1,1,nx);
    xyz(2,..)= span(-1,1,nx)(-,);
    xyz(3,..)= span(-1,1,nx)(-,-,);
    r= abs(xyz(1,..),xyz(2,..),xyz(3,..));
    theta= acos(xyz(3,..)/r);
    phi= atan(xyz(2,..),xyz(1,..)+(!r));
    y32= sin(theta)^2*cos(theta)*cos(2*phi);
    f= r*(1.+y32);
    r= theta= phi= [];

    write, "   test uses "+pr1((nx-1)*(ny-1)*(nz-1))+" cells";

    origin= [-1.0,-1.0,-1.0];
    delta= [2.0/(nx-1.0), 2.0/(nx-1.0), 2.0/(nx-1.0)];
    level= 1.0;
    colr= [0.0,1.0,0.2];
    tris= iso3cenre(origin,delta,f,level,colr);
    if(tris) {
      tris2= CollapseTri3d(tris);
    }

    tnew= slice2only(plane, tris2);
    clear3;
    if(tnew) {
      res= pl3dtrilists(tnew);
    }

    draw3, !making_movie;
    stdview3;
}

local slice2_precision;
/* DOCUMENT slice2_precision= precision
     Controls how slice2 (or slice2x) handles points very close to
     the slicing plane.  PRECISION should be a positive number or zero.
     Zero PRECISION means to clip exactly to the plane, with points
     exactly on the plane acting as if they were slightly on the side
     the normal points toward.  Positive PRECISION means that edges
     are clipped to parallel planes a distance PRECISION on either
     side of the given plane.  (Polygons lying entirely between these
     planes are completely discarded.)

     Default value is 0.0.

   SEE ALSO: slice2, slice2x
 */
slice2_precision= 0.0;

func _nuslice2_part(num,ntri,keep,dp,xyzc,normc,colrc,idsc,xyznew,normnew,colrnew,idsnew)
{
  ncut= keep(sum,);
  one= where(ncut == 1);
  two= where(ncut == 2);
  /* find the points where any edges of the polys cut the plane */
  mask0= (!keep) & keep(next);   /* off-->on */
  list0= where(mask0);
  if (numberof(list0)) {
    list= next(list0);
    dpl= dp(list0)(-,);
    dpu= dp(list)(-,);
    xyz0= (xyzc(,list0)*dpu-xyzc(,list)*dpl)/(dpu-dpl);
  }
  mask1= (!keep) & keep(prev);   /* on-->off */
  list1= where(mask1);
  if (numberof(list1)) {
    list= prev(list1);
    dpl= dp(list1)(-,);
    dpu= dp(list)(-,);
    xyz1= (xyzc(,list1)*dpu-xyzc(,list)*dpl)/(dpu-dpl);
  }

  /* form an index list xold which gives the indices in the original
   * xyzc list at the places in the new, sliced xyzc list */
  mask= keep+mask0+mask1;  /* 0, 1, or 2 */
  list= mask(psum);  /* index values in new list */
  xold= array(0, list(0));
  mlist= where(mask);
  xold(list(mlist))= mlist;
  dups= where(mask==2);
  if (numberof(dups)) xold(list(dups)-1)= dups;

  /* form the new, sliced xyzc vertex list */
  xyzc= xyzc(,xold);
  if (numberof(list0)) xyzc(,list(list0))= xyz0;
  if (numberof(list1)) xyzc(,list(list1)-mask(list1)+1)= xyz1;

  /* get the list of indices into nvertc (or valuec) for each of
   * the points in xyzc */
  ndxs= histogram(1+last)(1:-1);
  ndxs(1)+= 1;
  ndxs= ndxs(psum);
  /* compute the new number of vertices */
  nvertc= histogram(ndxs(xold));
}
