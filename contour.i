/*
 * $Id: contour.i,v 1.1.1.1 2005-09-18 22:08:00 dhmunro Exp $
 * These functions extract iso-surfaces and slicing planes
 * from 3D meshes.
 * The results can be displayed using functions in pl3gl.i.
 * Users ordinarily call functions in contour.i so those 
 * functions are autoloaded by the yorgl plugin.
 * The compiled functions are connected via cntrfunc.i, so 
 * require that cntrfunc.i be "pulled in" any time something in contour.i
 * is called.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

require,"tristruct.i";
require,"cntrfunc.i";

/* by default, make triangle strips when possible */
make_strip= 1;


func is_scalar_color(colr, name)
{
  if(is_void(colr)) {
    set_color= 0;
  } else {
    dimc= dimsof(colr);
    if( !(dimc(1) == 1 && (dimc(2) == 3 || dimc(2) == 4)) )
      error,"colr must be a vector of length 3 or 4 in "+name;
    set_color= 1;
  }
  return set_color;
}

func is_3_vector(var, msg)
{
  if(anyof(dimsof(var) != [1,3])) error(msg);
}

func is_compat(var1, var2, msg)
{
  if(is_void(var2)) return;
  if( anyof( dimsof(var1) != dimsof(var2) ) ) {
    error,msg;
  }
}

func chek_iso_reg(origin, delta, dimsv, var2, name, zcn=)
{
  if(zcn) zoff= 1;
  else zoff= 0;
  if(anyof(dimsof(origin) != [1,3])) error,"origin must be a 3 element vector in "+name;
  if(anyof(dimsof(delta) != [1,3])) error,"delta must be a 3 element vector in "+name;
  if(anyof(dimsv(2:) < 4-zoff)) error,"grid must be at least 4-by-4-by-4 in "+name;
  if(!is_void(var2)) if( anyof( dimsv != dimsof(var2) ) ) error,"var2 must have the same dimensions as var in "+name;
}

func chek_iso_crv(xyz, dimsv, var2, name, zcn=)
{
  dimsx= dimsof(xyz);
  if(is_void(zcn)) zcn= 0;
  if(dimsv(1) != 3) error,"var must be a 3-dimensional array in "+name;
  if(dimsx(1) != 4 || dimsx(2) != 3) error,"xyz must have leading dimension 3 in "+name;
  if(anyof(dimsx(3:5) < 4)) error,"grid must be at least 4-by-4-by-4 in "+name;
  if(zcn) {
    if( anyof(dimsv(2:4) != dimsx(3:5)-1) ) error,"var must have dimensions nx-1, ny-1, and nz-1 in "+name;
  } else {
    if( anyof(dimsv(2:4) != dimsx(3:5)) ) error,"xyz and var must have the same nx, ny, and nz in "+name;
  }
  if(!is_void(var2)) if( anyof(dimsv != dimsof(var2)) ) error,"var2 must have the same dimensions as var in "+name;
}

func get_chunk(sizes)
{
  chunk= CHKSIZ;
  if( chunk(1)*chunk(2)*chunk(3) >= sizes(1)*sizes(2)*sizes(3) ) {
    chunk(1)= sizes(1);
    chunk(2)= sizes(2);
    chunk(3)= sizes(3);
  }
  chunk= min(chunk, sizes);
  return chunk;
}

func iso3cenregngrd(origin, delta, var, level, colr, var2=)
/* DOCUMENT iso3cenregngrd(origin, delta, var, level, colr, var2=)

     Extract an iso-surface from a variable on a rectangular 3D grid.
     Input array var is point centered and NO guard points. 
     var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.

     Uses a 6 tetrahedron decomposition of each hexahedral zone.
     The result is a list of triangle arrays.

   SEE ALSO: iso3cenreg, iso3zcenregngrd, iso3cenregndx, iso3cencrv.
 */
{
  lst= iso3cenregbase(origin, delta, var, level, colr, 0, var2);
  return lst;
}

func iso3(origin, delta, var, level, colr, var2=)
/* DOCUMENT iso3(origin, delta, var, level, colr, var2=)

     Included for compatibility of earlier versions of yorgl.
     Please call iso3cenreg instead.

   SEE ALSO: iso3cenreg.
 */
{
  lst= iso3cenreg(origin, delta, var, level, colr, var2=var2);
  return lst;
}

func iso3cenreg(origin, delta, var, level, colr, var2=)
/* DOCUMENT iso3cenreg(origin, delta, var, level, colr, var2=)

     Extract an iso-surface from a variable on a rectangular 3D grid.
     Input array var is point centered and has one extra point on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.

     Uses a 6 tetrahedron decomposition of each hexahedral zone.
     The result is a list of triangle arrays.

   SEE ALSO: iso3cenregngrd, iso3zcenreg, iso3cenregndx, iso3cencrv.
 */
{
  lst= iso3cenregbase(origin, delta, var, level, colr, 1, var2);
  return lst;
}

func iso3cenregbase(origin, delta, var, level, colr, has_guard, var2)
/* worker routine for point centered data, not called by user */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  // set name of calling routine and the number of guard
  // on each side of the grid
  if(has_guard) {
    name= "iso3cenreg";
    nguard= 1;
  } else {
    name= "iso3cenregngrd";
    nguard= 0;
  }
  make_strip= 0;
  dimsv= dimsof(var);
  sizes= dimsv(2:4);
  chek_iso_reg, origin, delta, dimsv, var2, name, zcn=0;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);

  chunk= get_chunk(sizes);
  cnx= chunk(1);
  cny= chunk(2);
  cnz= chunk(3);
  set_color= is_scalar_color(colr, name);

  nulvar= [];
  /* number of cells in a chunk not counting guard cells */
  numZone= (cnx-3)*(cny-3)*(cnz-3);
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*numZone;
  maxVert= 3*maxTri;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= &array(0.0, 3, maxTri);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriArrayGrp(
          numTri= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, 3, maxTri),
          normals= &array(0.0, 3, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
         );

  xyztmp= array(0.0, 3, cnx, cny, cnz);
  grdtmp= array(0.0, 3, cnx, cny, cnz);
  vartmp= array(0.0, cnx, cny, cnz);
  v2tmp=  array(0.0, cnx, cny, cnz);
  flag=   array(char, cnx, cny, cnz);
  lst= [];

  /* The compiled functions use a chunk (specified by siz) that 
     includes guard points for the chunk. If the input array
     didn't have guard points, klo and khi will point 
     outside the input array when the chunk is at the border. 
     This is handled by the point fetching routine in the compiled code.
     The loops below run over the real cells in each chunk, so 
     nx-1-nguard is the upper limit (i.e. the last real cell in the 
     input array).
     The limits passed in the call includes both real and guard points for  
     the chunk.
   */
  for(k= 1+nguard; k <= nz-1-nguard; k += cnz-3) {
    klo= k-1; /* the lower guard point of the chunk */
    khi= min(nz+1-nguard, klo+cnz-1);
    for(j= 1+nguard; j <= ny-1-nguard; j += cny-3) {
      jlo= j-1;
      jhi= min(ny+1-nguard, jlo+cny-1);
      for(i= 1+nguard; i <= nx-1-nguard; i += cnx-3) {
        ilo= i-1;
        ihi= min(nx+1-nguard, ilo+cnx-1);
        siz= [ihi-ilo+1, jhi-jlo+1, khi-klo+1];
        /* norigin is the location of the first guard point for the 
           chunk, not the first real point. */
        norigin= origin+[ilo-1, jlo-1, klo-1]*delta;
        offsets= [ilo, jlo, klo, nx, ny, nz];
        if(has_guard) {
          res= ContourInitCartGrdPcen(siz,offsets,delta,norigin,var, &var2);
        } else {
          res= ContourInitCartPcen(siz,offsets,delta,norigin,var, &var2);
        }
        if(!res) {
          write,"chunk ",siz(1), siz(2), siz(3), " is too small";
          continue;
        }
        res= ContourTetArray(make_strip,siz,level,
	         vartmp, v2tmp, xyztmp, grdtmp, flag, &arr0);
        if(!res) continue;
        /* have an array of triangles. */
        numTri= arr0.numTri;
        xyzverts= (*(arr0.xyzverts))(,,1:numTri);
        normals= (*(arr0.normals))(,,1:numTri);
        cellIDs= (*(arr0.cellIDs))(1:numTri);
        if(!is_void(var2)) {
          /* need to be able to save one variable value at every
             triangle vertex */
          v2new= (*(arr0.var2))(,1:numTri);
        } else {
          v2new= [];
        }
        lst= TriArrayGrp(
            numTri= numTri,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            /* the next three are not used unless assembling strips */
	    triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
        if(set_color) {
          lst.colors= &float(colr);
        }
      }
    }
  }
  return lst;
}

func iso3zcenregngrd(origin, delta, var, level, colr, var2=)
/* DOCUMENT iso3zcenreg(origin, delta, var, level, colr, var2=)

     Extract an iso-surface from a variable on a rectangular 3D grid.
     Input array var is zone centered and has one extra zone on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.

     Uses a 6 tetrahedron decomposition of each hexahedral zone.
     The result is a list of triangle arrays.

   SEE ALSO: iso3zcenreg, iso3zcenregndx, iso3zcencrv.
 */
{
  lst= iso3zcenregbase(origin, delta, var, level, colr, 0, var2);
  return lst;
}

func iso3zcenreg(origin, delta, var, level, colr, var2=)
/* DOCUMENT iso3zcenreg(origin, delta, var, level, colr, var2=)

     Extract an iso-surface from a variable on a rectangular 3D grid.
     Input array var is zone centered and has one extra zone on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.

     Uses a 6 tetrahedron decomposition of each hexahedral zone.
     The result is a list of triangle arrays.

   SEE ALSO: iso3, iso3crv, iso3strip, iso3stripcrv, slice3stripcrv.
 */
{
  lst= iso3zcenregbase(origin, delta, var, level, colr, 1, var2);
  return lst;
}

func iso3zcenregbase(origin, delta, var, level, colr, has_guard, var2)
/* worker routine for zone centered data, not called by the user */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  // set name of calling routine and the number of guard cells
  // on each side of the grid
  if(has_guard) {
    name= "iso3zcenreg";
    nguard= 1;
  } else {
    name= "iso3zcenregngrd";
    nguard= 0;
  }
  make_strip= 0;
  dimsv= dimsof(var);
  sizes= dimsv(2:4)+1;  /* number of vertices surrounding the cells */
  chek_iso_reg, origin, delta, dimsv, var2, name, zcn=1;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);

  chunk= get_chunk(sizes);
  cnx= chunk(1);
  cny= chunk(2);
  cnz= chunk(3);
  set_color= is_scalar_color(colr, name);

  nulvar= [];
  /* number of cells not counting guard cells */
  numZone= (cnx-3)*(cny-3)*(cnz-3);
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*numZone;
  maxVert= 3*maxTri;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= &array(0.0, 3, maxTri);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriArrayGrp(
          numTri= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, 3, maxTri),
          normals= &array(0.0, 3, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
        );

  xyztmp= array(0.0, 3, cnx, cny, cnz);
  grdtmp= array(0.0, 3, cnx, cny, cnz);
  vartmp= array(0.0, cnx, cny, cnz);
  v2tmp=  array(0.0, cnx, cny, cnz);
  flag=   array(char, cnx, cny, cnz);
  lst= [];

  /* The input array may or may not have guard points on all sides. 
     The compiled functions use a chunk (specified by siz) that 
     includes guard points for the chunk (a chunk has guard points
     even if the input array does not).
     The loops below run over the real cells in each chunk, so nx-1-nguard 
     is the upper limit (i.e. the last real cell in the input array).
     The limits passed in the call include both real and guard points for  
     the chunk. */
  for(k= 1+nguard; k <= nz-1-nguard; k += cnz-3) {
    klo= k-1; /* the lower guard point of the chunk */
    khi= min(nz+1-nguard, klo+cnz-1); // the upper guard point of the chunk
    for(j= 1+nguard; j <= ny-1-nguard; j += cny-3) {
      jlo= j-1;
      jhi= min(ny+1-nguard, jlo+cny-1);
      for(i= 1+nguard; i <= nx-1-nguard; i += cnx-3) {
        ilo= i-1;
        ihi= min(nx+1-nguard, ilo+cnx-1);
        siz= [ihi-ilo+1, jhi-jlo+1, khi-klo+1];
        /* norigin is the location of the first guard point for the 
           chunk, not the first real point. */
        norigin= origin+[ilo-1, jlo-1, klo-1]*delta;
        offsets= [ilo, jlo, klo, nx, ny, nz];
        if(has_guard) {
          res= ContourInitCartGrdZcen(siz,offsets,delta,norigin,var, &var2);
        } else {
          res= ContourInitCartZcen(siz,offsets,delta,norigin,var, &var2);
        }
        if(!res) {
          write,"chunk ",siz(1), siz(2), siz(3), " is too small";
          continue;
        }
        res= ContourTetArray(make_strip,siz,level,
	         vartmp, v2tmp, xyztmp, grdtmp, flag, &arr0);
        if(!res) continue;
        /* have an array of triangles. */
        numTri= arr0.numTri;
        xyzverts= (*(arr0.xyzverts))(,,1:numTri);
        normals= (*(arr0.normals))(,,1:numTri);
        cellIDs= (*(arr0.cellIDs))(1:numTri);
        if(!is_void(var2)) {
          /* need to be able to save one variable value at every
             triangle vertex */
          v2new= (*(arr0.var2))(,1:numTri);
        } else {
          v2new= [];
        }
        lst= TriArrayGrp(
            numTri= numTri,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            /* the next three are not used unless assembling strips */
            triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
        if(set_color) {
          lst.colors= &float(colr);
        }
      }
    }
  }
  return lst;
}

func iso3ndx(origin, delta, var, level, colr, var2=)
/* DOCUMENT iso3ndx(origin, delta, var, level, colr, var2=)

     Included for compatibility with earlier versions of yorgl.
     Please call iso3cenregndx instead. 

     Uses a 6 tetrahedron decomposition of each hexahedral zone.

   SEE ALSO: iso3cenregndx.
 */
{
  lst= iso3cenregndx(origin, delta, var, level, colr, var2=var2);
  return lst;
}

func iso3cenregndx(origin, delta, var, level, colr, var2=)
/* DOCUMENT iso3cenregndx(origin, delta, var, level, colr, var2=)

     Extract an iso-surface from a variable on a rectangular 3D grid.
     Input array var is point centered and has one extra point on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.
     Triangles are specified by indices into vectors of coordinates, 
     gradients, and auxiliary variables. 

     Uses a 6 tetrahedron decomposition of each hexahedral zone.

   SEE ALSO: iso3cenreg, iso3zcenregndx, iso3cencrv, iso3cencrvndx.
 */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  name= "iso3cenregndx";
  make_strip= 0;
  dimsv= dimsof(var);
  sizes= dimsv(2:4);
  chek_iso_reg, origin, delta, dimsv, var2, name, zcn=0;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);

  chunk= get_chunk(sizes);
  cnx= chunk(1);
  cny= chunk(2);
  cnz= chunk(3);
  set_color= is_scalar_color(colr, name);
  nulvar= [];
  /* number of cells not counting guard cells */
  numZone= (cnx-3)*(cny-3)*(cnz-3);
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*numZone;
  maxVert= 3*maxTri;
  maxEdg= 3*cnx*cny*cnz;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       unique vertex */
    v2= &array(0.0, maxEdg);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriVertexGrp(
          numTri= 0,
          numEdg= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, maxEdg),
          normals= &array(0.0, 3, maxEdg),
          ptndx= &array(0, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          next= &nulvar
        );

  xyztmp= array(0.0, 3, cnx, cny, cnz);
  grdtmp= array(0.0, 3, cnx, cny, cnz);
  vartmp= array(0.0, cnx, cny, cnz);
  v2tmp=  array(0.0, cnx, cny, cnz);
  flag=   array(char, cnx, cny, cnz);
  /* there is one index into the vertex array per real edge in the grid */
  ndx= array(long, maxEdg);
  lst= [];

  /* The input array has guard points on all sides. The compiled functions
     use a chunk (specified by siz) that includes guard points for the chunk.
     The loops below run over the real cells in each chunk, so nx-2 is the 
     upper limit (i.e. the last real cell in the input array).
     The limits passed in the call include both real and guard points for  
     the chunk. */
  for(k= 2; k <= nz-2; k += cnz-3) {
    klo= k-1; /* the lower guard point of the chunk */
    khi= min(nz, klo+cnz-1);
    for(j= 2; j <= ny-2; j += cny-3) {
      jlo= j-1;
      jhi= min(ny, jlo+cny-1);
      for(i= 2; i <= nx-2; i += cnx-3) {
        ilo= i-1;
	    ihi= min(nx, ilo+cnx-1);
        siz= [ihi-ilo+1, jhi-jlo+1, khi-klo+1];
        /* norigin is the location of the first guard point for the 
           chunk, not the first real point. */
        norigin= origin+[ilo-1, jlo-1, klo-1]*delta;
        offsets= [ilo, jlo, klo, nx, ny, nz];
        res= ContourInitCartGrdPcenNdx(siz,offsets,delta,norigin,var, &var2);
        if(!res) {
          write,"chunk ",siz(1), siz(2), siz(3), " is too small";
          continue;
        }
        res= ContourTetArrayNdx(make_strip,siz,level,
	         vartmp, v2tmp, xyztmp, grdtmp, flag, ndx, &arr0);
        if(!res) continue;
        /* have an array of triangles specified by indices into vertex arrays. */
        numTri= arr0.numTri;
        numEdg= arr0.numEdg;
        xyzverts= (*(arr0.xyzverts))(,1:numEdg);
        normals= (*(arr0.normals))(,1:numEdg);
        cellIDs= (*(arr0.cellIDs))(1:numTri);
        ptndx= (*(arr0.ptndx))(,1:numTri);
        if(!is_void(var2)) {
          /* need to be able to save one variable value at every
             triangle vertex */
          v2new= (*(arr0.var2))(1:numEdg);
        } else {
          v2new= [];
        }
        lst= TriVertexGrp(
            numTri= numTri,
            numEdg= numEdg,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            ptndx= &ptndx,
            /* the next three are not used unless assembling strips */
            triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
        if(set_color) {
          lst.colors= &float(colr);
        }
      }
    }
  }
  return lst;
}

func iso3zcenregndx(origin, delta, var, level, colr, var2=)
/* DOCUMENT iso3zcenregndx(origin, delta, var, level, colr, var2=)

     Extract an iso-surface from a variable on a rectangular 3D grid.
     Input array var is zone centered and has one extra zone on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.
     Triangles are specified by indices into vectors of coordinates, 
     gradients, and auxiliary variables. 

     Uses a 6 tetrahedron decomposition of each hexahedral zone.

   SEE ALSO: iso3zcenreg, iso3cenregndx, iso3zcencrv, iso3zcencrvndx.
 */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  name= "iso3zcenregndx";
  make_strip= 0;
  dimsv= dimsof(var);
  sizes= dimsv(2:4)+1;  /* number of vertices surrounding the cells */
  chek_iso_reg, origin, delta, dimsv, var2, name, zcn=1;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);

  chunk= get_chunk(sizes);
  cnx= chunk(1);
  cny= chunk(2);
  cnz= chunk(3);
  set_color= is_scalar_color(colr, name);

  nulvar= [];
  /* number of cells not counting guard cells */
  numZone= (cnx-3)*(cny-3)*(cnz-3);
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*numZone;
  maxVert= 3*maxTri;
  maxEdg= 3*cnx*cny*cnz;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       unique vertex */
    v2= &array(0.0, maxEdg);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriVertexGrp(
          numTri= 0,
          numEdg= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, maxEdg),
          normals= &array(0.0, 3, maxEdg),
          ptndx= &array(0, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
        );

  xyztmp= array(0.0, 3, cnx, cny, cnz);
  grdtmp= array(0.0, 3, cnx, cny, cnz);
  vartmp= array(0.0, cnx, cny, cnz);
  v2tmp=  array(0.0, cnx, cny, cnz);
  flag=   array(char, cnx, cny, cnz);
  /* there is one index into the vertex array per real edge in the grid */
  ndx= array(long, maxEdg);
  lst= [];

  /* The input array has guard points on all sides. The compiled functions
     use a chunk (specified by siz) that includes guard points for the chunk.
     The loops below run over the real cells in each chunk, so nx-2 is the 
     upper limit (i.e. the last real cell in the input array).
     The limits passed in the call include both real and guard points for  
     the chunk. */
  for(k= 2; k <= nz-2; k += cnz-3) {
    klo= k-1; /* the lower guard point of the chunk */
    khi= min(nz, klo+cnz-1);
    for(j= 2; j <= ny-2; j += cny-3) {
      jlo= j-1;
      jhi= min(ny, jlo+cny-1);
      for(i= 2; i <= nx-2; i += cnx-3) {
        ilo= i-1;
	    ihi= min(nx, ilo+cnx-1);
        siz= [ihi-ilo+1, jhi-jlo+1, khi-klo+1];
        /* norigin is the location of the first guard point for the 
           chunk, not the first real point. */
        norigin= origin+[ilo-1, jlo-1, klo-1]*delta;
        offsets= [ilo, jlo, klo, nx, ny, nz];
        res= ContourInitCartGrdZcenNdx(siz,offsets,delta,norigin,var, &var2);
        if(!res) {
          write,"chunk ",siz(1), siz(2), siz(3), " is too small";
          continue;
        }
        res= ContourTetArrayNdx(make_strip,siz,level,
	         vartmp, v2tmp, xyztmp, grdtmp, flag, ndx, &arr0);
        if(!res) continue;
        /* have an array of triangles specified by indices into vertex arrays. */
        numTri= arr0.numTri;
        numEdg= arr0.numEdg;
        xyzverts= (*(arr0.xyzverts))(,1:numEdg);
        normals= (*(arr0.normals))(,1:numEdg);
        cellIDs= (*(arr0.cellIDs))(1:numTri);
        ptndx= (*(arr0.ptndx))(,1:numTri);
        if(!is_void(var2)) {
          /* need to be able to save one variable value at every
             triangle vertex */
          v2new= (*(arr0.var2))(1:numEdg);
        } else {
          v2new= [];
        }
        lst= TriVertexGrp(
            numTri= numTri,
            numEdg= numEdg,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            ptndx= &ptndx,
            /* the next three are not used unless assembling strips */
            triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
        if(set_color) {
          lst.colors= &float(colr);
        }
      }
    }
  }
  return lst;
}




func iso3hex(xyz, grad, hexndx, var, level, colr, var2=)
/* DOCUMENT iso3hex(xyz, grad, hexndx, var, level, colr, var2=)
     Extract an iso-surface from a variable on an arbitrarily
     connected 3D grid of hexahedra.
     xyz, grad, and var are coordinates, gradients, and
     variable values at the points of the grid.
     var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     hexndx has 8 indices into the xyz etc. arrays for 
     each zone.

     Uses a 6 tetrahedron decomposition of each hexahedral zone.
     The result is a list of triangle arrays.

   SEE ALSO: iso3cenreg, iso3zcenreg, iso3cencrv, iso3_tree.
 */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  name= "iso3hex";
  make_strip= 0;
  dimsv= dimsof(var);
  ASSERT, (dimsv(1) == 1 && dimsv(2) >= 8), "The variable array must contain at least 8 points in "+name;
  is_compat, var, var2, "var2 must have the same dimensions as var in "+name;
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3 && dimsx(3) == dimsv(2)),
	  "The xyz array must be 3 by N-points in iso3hex";
  ASSERT, allof(dimsx == dimsof(grad)),
	  "The gradient array must be 3 by N-points in "+name;
  dimsh= dimsof(hexndx);
  ASSERT, (dimsh(1) == 2 && dimsh(2) == 8), 
	  "The index array must be 8 by N-zones in "+name;
  nzone= dimsh(3);
  set_color= is_scalar_color(colr, name);

  nulvar= [];
  /* The input array can easily be chunked if the scratch array
     for triangles would be too large */
  chk= CHKSIZ(1)*CHKSIZ(2)*CHKSIZ(3);
  if(chk > nzone) chk= nzone;
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*chk;
  maxVert= 3*maxTri;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= &array(0.0, 3, maxTri);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriArrayGrp(
          numTri= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, 3, maxTri),
          normals= &array(0.0, 3, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
        );

  lst= [];

  /* The input variables include gradients at all points. There is 
     no need for guard cells. */
  for(i= 1; i <= nzone; i += chk) {
    num= min(chk, nzone-i+1);
    res= ContourTetHex(level, i, num, xyz, grad, hexndx, var, &var2, &arr0);
    if(!res) continue;
    /* have an array of triangles. */
    numTri= arr0.numTri;
    xyzverts= (*(arr0.xyzverts))(,,1:numTri);
    normals= (*(arr0.normals))(,,1:numTri);
    cellIDs= (*(arr0.cellIDs))(1:numTri);
    if(!is_void(var2)) {
      /* need to be able to save one variable value at every
         triangle vertex */
      v2new= (*(arr0.var2))(,1:numTri);
    } else {
      v2new= [];
    }
    lst= TriArrayGrp(
            numTri= numTri,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            /* the next three are not used unless assembling strips */
            triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
    if(set_color) {
      lst.colors= &float(colr);
    }
  }
  return lst;
}




func iso3cencrv(xyz, var, level, colr, var2=)
/* DOCUMENT iso3cencrv(xyz
     Extract an iso-surface from a variable on a curvilinear 3D grid.
     Input array var is point centered and has one extra point on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.

     Uses a 6 tetrahedron decomposition of each hexahedral zone.
     The result is a list of triangle arrays.

   SEE ALSO: iso3cenreg, iso3zcencrv, iso3zcencrvndx.
 */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  name= "iso3cencrv";
  make_strip= 0;
  dimsv= dimsof(var);
  sizes= dimsv(2:4);
  chek_iso_crv, xyz, dimsv, var2, name, zcn=0;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);

  chunk= get_chunk(sizes);
  cnx= chunk(1);
  cny= chunk(2);
  cnz= chunk(3);
  set_color= is_scalar_color(colr, name);

  nulvar= [];
  /* number of cells not counting guard cells */
  numZone= (cnx-3)*(cny-3)*(cnz-3);
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*numZone;
  maxVert= 3*maxTri;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= &array(0.0, 3, maxTri);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriArrayGrp(
          numTri= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, 3, maxTri),
          normals= &array(0.0, 3, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
        );

  xyztmp= array(0.0, 3, cnx, cny, cnz);
  grdtmp= array(0.0, 3, cnx, cny, cnz);
  vartmp= array(0.0, cnx, cny, cnz);
  v2tmp=  array(0.0, cnx, cny, cnz);
  flag=   array(char, cnx, cny, cnz);
  lst= [];

  /* The input array has guard points on all sides. The compiled functions
     use a chunk (specified by siz) that includes guard points for the chunk.
     The loops below run over the real cells in each chunk, so nx-2 is the 
     upper limit (i.e. the last real cell in the input array).
     The limits passed in the call include both real and guard points for  
     the chunk. */
  for(k= 2; k <= nz-2; k += cnz-3) {
    klo= k-1; /* the lower guard point of the chunk */
    khi= min(nz, klo+cnz-1);
    for(j= 2; j <= ny-2; j += cny-3) {
      jlo= j-1;
      jhi= min(ny, jlo+cny-1);
      for(i= 2; i <= nx-2; i += cnx-3) {
        ilo= i-1;
	    ihi= min(nx, ilo+cnx-1);
        siz= [ihi-ilo+1, jhi-jlo+1, khi-klo+1];
        /* norigin is the location of the first guard point for the 
           chunk, not the first real point. */
        norigin= origin+[ilo-1, jlo-1, klo-1]*delta;
        offsets= [ilo, jlo, klo, nx, ny, nz];
        res= ContourInitCrvGrdPcen(siz,offsets,xyz,var, &var2);
        if(!res) {
          write,"chunk ",siz(1), siz(2), siz(3), " is too small";
          continue;
        }
        res= ContourTetArray(make_strip,siz,level,
	         vartmp, v2tmp, xyztmp, grdtmp, flag, &arr0);
        if(!res) continue;
        /* have an array of triangles. */
        numTri= arr0.numTri;
        xyzverts= (*(arr0.xyzverts))(,,1:numTri);
        normals= (*(arr0.normals))(,,1:numTri);
        cellIDs= (*(arr0.cellIDs))(1:numTri);
        if(!is_void(var2)) {
          /* need to be able to save one variable value at every
             triangle vertex */
          v2new= (*(arr0.var2))(,1:numTri);
        } else {
          v2new= [];
        }
        lst= TriArrayGrp(
            numTri= numTri,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            /* the next three are not used unless assembling strips */
            triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
        if(set_color) {
          lst.colors= &float(colr);
        }
      }
    }
  }
  return lst;
}

func iso3zcencrv(xyz, var, level, colr, var2=)
/* DOCUMENT iso3zcencrv(xyz, var, level, colr, var2=)

     Extract an iso-surface from a variable on a curvilinear 3D grid.
     Input array var is zone centered and has one extra zone on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.

     Uses a 6 tetrahedron decomposition of each hexahedral zone.
     The result is a list of triangle arrays.

   SEE ALSO: iso3zcenreg, iso3cencrv, iso3zcencrvndx.
 */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  name= "iso3zcencrv";
  make_strip= 0;
  dimsv= dimsof(var);
  sizes= dimsv(2:4)+1;  /* number of vertices surrounding the cells */
  chek_iso_crv, xyz, dimsv, var2, name, zcn=1;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);

  chunk= get_chunk(sizes);
  cnx= chunk(1);
  cny= chunk(2);
  cnz= chunk(3);
  set_color= is_scalar_color(colr, name);

  nulvar= [];
  /* number of cells not counting guard cells */
  numZone= (cnx-3)*(cny-3)*(cnz-3);
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*numZone;
  maxVert= 3*maxTri;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= &array(0.0, 3, maxTri);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriArrayGrp(
          numTri= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, 3, maxTri),
          normals= &array(0.0, 3, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
        );

  xyztmp= array(0.0, 3, cnx, cny, cnz);
  grdtmp= array(0.0, 3, cnx, cny, cnz);
  vartmp= array(0.0, cnx, cny, cnz);
  v2tmp=  array(0.0, cnx, cny, cnz);
  flag=   array(char, cnx, cny, cnz);
  lst= [];

  /* The input array has guard points on all sides. The compiled functions
     use a chunk (specified by siz) that includes guard points for the chunk.
     The loops below run over the real cells in each chunk, so nx-2 is the 
     upper limit (i.e. the last real cell in the input array).
     The limits passed in the call include both real and guard points for  
     the chunk. */
  for(k= 2; k <= nz-2; k += cnz-3) {
    klo= k-1; /* the lower guard point of the chunk */
    khi= min(nz, klo+cnz-1);
    for(j= 2; j <= ny-2; j += cny-3) {
      jlo= j-1;
      jhi= min(ny, jlo+cny-1);
      for(i= 2; i <= nx-2; i += cnx-3) {
        ilo= i-1;
	    ihi= min(nx, ilo+cnx-1);
        siz= [ihi-ilo+1, jhi-jlo+1, khi-klo+1];
        /* norigin is the location of the first guard point for the 
           chunk, not the first real point. */
        norigin= origin+[ilo-1, jlo-1, klo-1]*delta;
        offsets= [ilo, jlo, klo, nx, ny, nz];
        res= ContourInitCrvGrdZcen(siz,offsets,xyz,var, &var2);
        if(!res) {
          write,"chunk ",siz(1), siz(2), siz(3), " is too small";
          continue;
        }
        res= ContourTetArray(make_strip,siz,level,
	         vartmp, v2tmp, xyztmp, grdtmp, flag, &arr0);
        if(!res) continue;
        /* have an array of triangles. */
        numTri= arr0.numTri;
        xyzverts= (*(arr0.xyzverts))(,,1:numTri);
        normals= (*(arr0.normals))(,,1:numTri);
        cellIDs= (*(arr0.cellIDs))(1:numTri);
        if(!is_void(var2)) {
          /* need to be able to save one variable value at every
             triangle vertex */
          v2new= (*(arr0.var2))(,1:numTri);
        } else {
          v2new= [];
        }
        lst= TriArrayGrp(
            numTri= numTri,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            /* the next three are not used unless assembling strips */
            triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
        if(set_color) {
          lst.colors= &float(colr);
        }
      }
    }
  }
  return lst;
}

func iso3cencrvndx(xyz, var, level, colr, var2=)
/* DOCUMENT iso3cencrvndx(xyz, var, level, colr, var2=)

     Extract an iso-surface from a variable on a curvilinear 3D grid.
     Input array var is point centered and has one extra point on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.
     Triangles are specified by indices into vectors of coordinates, 
     gradients, and auxiliary variables. 

     Uses a 6 tetrahedron decomposition of each hexahedral zone.

   SEE ALSO: iso3cenregndx, iso3cencrv, iso3zcencrvndx.
 */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  name= "iso3cencrvndx";
  make_strip= 0;
  dimsv= dimsof(var);
  sizes= dimsv(2:4);
  chek_iso_crv, xyz, dimsv, var2, name, zcn=0;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);

  chunk= get_chunk(sizes);
  cnx= chunk(1);
  cny= chunk(2);
  cnz= chunk(3);
  set_color= is_scalar_color(colr, name);

  nulvar= [];
  /* number of cells not counting guard cells */
  numZone= (cnx-3)*(cny-3)*(cnz-3);
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*numZone;
  maxVert= 3*maxTri;
  maxEdg= 3*cnx*cny*cnz;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       unique vertex */
    v2= &array(0.0, maxEdg);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriVertexGrp(
          numTri= 0,
          numEdg= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, maxEdg),
          normals= &array(0.0, 3, maxEdg),
          ptndx= &array(0, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
        );

  xyztmp= array(0.0, 3, cnx, cny, cnz);
  grdtmp= array(0.0, 3, cnx, cny, cnz);
  vartmp= array(0.0, cnx, cny, cnz);
  v2tmp=  array(0.0, cnx, cny, cnz);
  flag=   array(char, cnx, cny, cnz);
  /* there is one index into the vertex array per real edge in the grid */
  ndx= array(long, maxEdg);
  lst= [];

  /* The input array has guard points on all sides. The compiled functions
     use a chunk (specified by siz) that includes guard points for the chunk.
     The loops below run over the real cells in each chunk, so nx-2 is the 
     upper limit (i.e. the last real cell in the input array).
     The limits passed in the call include both real and guard points for  
     the chunk. */
  for(k= 2; k <= nz-2; k += cnz-3) {
    klo= k-1; /* the lower guard point of the chunk */
    khi= min(nz, klo+cnz-1);
    for(j= 2; j <= ny-2; j += cny-3) {
      jlo= j-1;
      jhi= min(ny, jlo+cny-1);
      for(i= 2; i <= nx-2; i += cnx-3) {
        ilo= i-1;
	    ihi= min(nx, ilo+cnx-1);
        siz= [ihi-ilo+1, jhi-jlo+1, khi-klo+1];
        /* norigin is the location of the first guard point for the 
           chunk, not the first real point. */
        norigin= origin+[ilo-1, jlo-1, klo-1]*delta;
        offsets= [ilo, jlo, klo, nx, ny, nz];
        res= ContourInitCrvGrdPcenNdx(siz,offsets,xyz,var, &var2);
        if(!res) {
          write,"chunk ",siz(1), siz(2), siz(3), " is too small";
          continue;
        }
        res= ContourTetArrayNdx(make_strip,siz,level,
	         vartmp, v2tmp, xyztmp, grdtmp, flag, ndx, &arr0);
        if(!res) continue;
        /* have an array of triangles specified by indices into vertex arrays. */
        numTri= arr0.numTri;
        numEdg= arr0.numEdg;
        xyzverts= (*(arr0.xyzverts))(,1:numEdg);
        normals= (*(arr0.normals))(,1:numEdg);
        cellIDs= (*(arr0.cellIDs))(1:numTri);
        ptndx= (*(arr0.ptndx))(,1:numTri);
        if(!is_void(var2)) {
          /* need to be able to save one variable value at every
             triangle vertex */
          v2new= (*(arr0.var2))(1:numEdg);
        } else {
          v2new= [];
        }
        lst= TriVertexGrp(
            numTri= numTri,
            numEdg= numEdg,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            ptndx= &ptndx,
            /* the next three are not used unless assembling strips */
            triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
        if(set_color) {
          lst.colors= &float(colr);
        }
      }
    }
  }
  return lst;
}

func iso3zcencrvndx(xyz, var, level, colr, var2=)
/* DOCUMENT iso3zcencrvndx(xyz, var, level, colr, var2=)

     Extract an iso-surface from a variable on a curvilinear 3D grid.
     Input array var is zone centered and has one extra zone on
     each side. var2, if suplied, is an auxiliary variable with
     the same centering and size as var.
     The result is a list of groups of triangle arrays. 
     If requested, triangle strips may extend across multiple zones.
     Triangles are specified by indices into vectors of coordinates, 
     gradients, and auxiliary variables. 

     Uses a 6 tetrahedron decomposition of each hexahedral zone.

   SEE ALSO: iso3zcenregndx, iso3zcencrv, iso3cencrvndx.
 */
{
  extern CHKSIZ;
  local make_strip;
  local numTri, nx, ny, nz, lst;

  name= "iso3zcencrvndx";
  make_strip= 0;
  dimsv= dimsof(var);
  sizes= dimsv(2:4)+1;  /* number of vertices surrounding the cells */
  chek_iso_crv, xyz, dimsv, var2, name, zcn=1;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);

  chunk= get_chunk(sizes);
  cnx= chunk(1);
  cny= chunk(2);
  cnz= chunk(3);
  set_color= is_scalar_color(colr, name);

  nulvar= [];
  /* number of cells not counting guard cells */
  numZone= (cnx-3)*(cny-3)*(cnz-3);
  /* note: this triangle count assumes that points on edges
     introduced by the 6 tet decomposition are removed */
  maxTri= 6*numZone;
  maxVert= 3*maxTri;
  maxEdg= 3*cnx*cny*cnz;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       unique vertex */
    v2= &array(0.0, maxEdg);
  } else {
    v2= &[];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriVertexGrp(
          numTri= 0,
          numEdg= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, maxEdg),
          normals= &array(0.0, 3, maxEdg),
          ptndx= &array(0, 3, maxTri),
          var2= v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
        );

  xyztmp= array(0.0, 3, cnx, cny, cnz);
  grdtmp= array(0.0, 3, cnx, cny, cnz);
  vartmp= array(0.0, cnx, cny, cnz);
  v2tmp=  array(0.0, cnx, cny, cnz);
  flag=   array(char, cnx, cny, cnz);
  /* there is one index into the vertex array per real edge in the grid */
  ndx= array(long, maxEdg);
  lst= [];

  /* The input array has guard points on all sides. The compiled functions
     use a chunk (specified by siz) that includes guard points for the chunk.
     The loops below run over the real cells in each chunk, so nx-2 is the 
     upper limit (i.e. the last real cell in the input array).
     The limits passed in the call include both real and guard points for  
     the chunk. */
  for(k= 2; k <= nz-2; k += cnz-3) {
    klo= k-1; /* the lower guard point of the chunk */
    khi= min(nz, klo+cnz-1);
    for(j= 2; j <= ny-2; j += cny-3) {
      jlo= j-1;
      jhi= min(ny, jlo+cny-1);
      for(i= 2; i <= nx-2; i += cnx-3) {
        ilo= i-1;
	    ihi= min(nx, ilo+cnx-1);
        siz= [ihi-ilo+1, jhi-jlo+1, khi-klo+1];
        /* norigin is the location of the first guard point for the 
           chunk, not the first real point. */
        norigin= origin+[ilo-1, jlo-1, klo-1]*delta;
        offsets= [ilo, jlo, klo, nx, ny, nz];
        res= ContourInitCrvGrdZcenNdx(siz,offsets,xyz,var, &var2);
        if(!res) {
          write,"chunk ",siz(1), siz(2), siz(3), " is too small";
          continue;
        }
        res= ContourTetArrayNdx(make_strip,siz,level,
             vartmp, v2tmp, xyztmp, grdtmp, flag, ndx, &arr0);
        if(!res) continue;
        /* have an array of triangles specified by indices into vertex arrays. */
        numTri= arr0.numTri;
        numEdg= arr0.numEdg;
        xyzverts= (*(arr0.xyzverts))(,1:numEdg);
        normals= (*(arr0.normals))(,1:numEdg);
        cellIDs= (*(arr0.cellIDs))(1:numTri);
        ptndx= (*(arr0.ptndx))(,1:numTri);
        if(!is_void(var2)) {
          /* need to be able to save one variable value at every
             triangle vertex */
          v2new= (*(arr0.var2))(1:numEdg);
        } else {
          v2new= [];
        }
        lst= TriVertexGrp(
            numTri= numTri,
            numEdg= numEdg,
            xyzverts= &xyzverts,
            normals= &normals,
            var2= &v2new,
            cellIDs= &cellIDs,
            ptndx= &ptndx,
            /* the next three are not used unless assembling strips */
            triEdg= &nulvar,
            triStart= &nulvar,
            nTris= &nulvar,
            next= &lst);
        if(set_color) {
          lst.colors= &float(colr);
        }
      }
    }
  }
  return lst;
}


func mak_isotree(var)
/* DOCUMENT mak_isotree(var)

     Build an octree for use in making iso-surfaces.

   SEE ALSO: iso3_tree, iso3_treecrv, iso3_treevarr, slice_tree
 */
{
  local nx, ny, nz, i, j, k;

  dimsv= dimsof(var);
  /* NOTE: the variable must have guard cells */
  ASSERT, (dimsv(1) == 3 && allof(dimsv(2:4) > 3)), 
      "The variable array must be nx-by-ny-by-nz in mak_isotree";
  /* NOTE: var is assumed to be point centered */
  sizes= dimsv(2:4);
  rlsiz= sizes-2;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);
  itop= nx-1;
  jtop= ny-1;
  ktop= nz-1;

  /* NOTE: chunk is the number of points in a chunk - the 
     number of cells is one less. A chunk never includes
     either the first or last points in var. */
  chunk= get_chunk(rlsiz);
  tree= [];
  nulvar= [];
  for(k= 2; k < ktop; k += chunk(3)-1) {
    khi= min(ktop, k+chunk(3)-1);
    for(j= 2; j < jtop; j += chunk(2)-1) {
      jhi= min(jtop, j+chunk(2)-1);
      for(i= 2; i < itop; i += chunk(1)-1) {
        ihi= min(itop, i+chunk(1)-1);
        /* chk is the number of vertices in this chunk */
        chk= [ihi-i+1, jhi-j+1, khi-k+1];
        /* determine the number of levels in the octree */
        maxdepth= 1;
        num= chk(max)-1;  /* number of cells */
        while(num > 1) {
          num= (num+1)/2;
          maxdepth++;
        }
        /* Compute the number of "cells" in each direction for all levels
           of the tree. 
        */
        numx= numy= numz= array(0, maxdepth);
        /* nnx, nny, nnz are numbers of cells */
        nnx= chk(1)-1; nny= chk(2)-1; nnz= chk(3)-1;
        for(id= 1; id <= maxdepth; id++) {
          numx(id)= nnx;
          numy(id)= nny;
          numz(id)= nnz;
          nnx= (nnx+1)/2;
          nny= (nny+1)/2;
          nnz= (nnz+1)/2;
        }
        /* Compute C-style offsets into the range array for the start
           of each level and the total size of the range array.
      	   An element of the range is a struct with a min and a max
           value.
           A single range array contains all levels. The offsets
           are counts of range elements, not counts of the doubles
       	   that make up ranges.
        */
        lens= numx*numy*numz;
        /* NOTE: C-style indices into the range array */
        offsets= lens(psum);
        if(numberof(lens) <= 1) {
          offsets= array(0, 1);
        } else {
          offsets= grow(0, offsets(:-1));
        }
        ranges= array(OctRange, lens(sum));
        trsiz= array(0, 3, maxdepth);
        trsiz(1,)= numx;
        trsiz(2,)= numy;
        trsiz(3,)= numz;
        tree= OctTree(
            maxdepth= maxdepth,
            size= &sizes,
            chunk= &chk,
            start= &[i-1,j-1,k-1], /* C-style start in the full array */
            trsiz= &trsiz,
            offsets= &offsets,
            ranges= &ranges,
            next= &tree
        );
        MakeContourTree, var, &tree;
      }
    }
  }
  return tree;
}

func slice_tree(sizes, origin, delta, point, pl_normal, var2, guard=)
/* DOCUMENT slice_tree(sizes, origin, delta, point, pl_normal, var2, guard=)

     Extract a slicing plane on a rectangular 3D grid.
     If var2 is not null, return it's value at each vertex of the
     slice plane.
     sizes= [nx,ny,nz] is the grid size.

   SEE ALSO: slice_treecrv, .
 */
{
  local numTri, nVert, nStrip;

  name= "slice_tree";
  /* assume guard cells by default!!! */
  if(is_void(guard)) guard= 1;
  if(guard) {
    ASSERT, allof(sizes(1:3) >= 4), 
        "The variable array must be nx-by-ny-by-nz in "+name;
  } else {
    ASSERT, allof(sizes(1:3) >= 2), 
        "The variable array must be nx-by-ny-by-nz in "+name;
  }
  is_3_vector, origin, "origin must be a 3 element vector in "+name;
  is_3_vector, delta, "delta must be a 3 element vector in "+name;
  is_3_vector, point, "point must be a 3 element vector in "+name;
  is_3_vector, pl_normal, "pl_normal must be a 3 element vector in "+name;

  /* NOTE: var is assumed to be point centered */
  if(guard) {
    rlsiz= sizes-2;
    ibot= 2; itop= sizes(1)-1;
    jbot= 2; jtop= sizes(2)-1;
    kbot= 2; ktop= sizes(3)-1;
  } else {
    rlsiz= sizes;
    ibot= 1; itop= sizes(1);
    jbot= 1; jtop= sizes(2);
    kbot= 1; ktop= sizes(3);
  }
  /* NOTE: chunk is the number of points in a chunk - the 
     number of cells is one less. A chunk never includes
     either the first or last points in var. */
  chunk= get_chunk(rlsiz);

  /* arr0 must be big enough for the worst possible case - 
     i.e. several triangles in each zone.
	 This can be huge unless the input array is handled in chunks.
	 No chunk can be bigger than the last chunk in the list
	 (the one that starts at [0,0,0]).
  */
//  trptr= &tree;
  numZone= (chunk(1)-1)*(chunk(2)-1)*(chunk(3)-1);
  maxTri= 5*numZone;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= &array(0.0, 3, maxTri);
  } else {
    v2= &[];
  }
  arr0= TriArrayGrp(
          numTri= 0,
	  cellIDs= &array(0, maxTri),
	  xyzverts= &array(0.0, 3, 3, maxTri),
	  normals= &array(0.0, 3, 3, maxTri),
	  var2= v2,
	  triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
         );

  lst= [];
  for(k= kbot; k < ktop; k += chunk(3)-1) {
    khi= min(ktop, k+chunk(3)-1);
    for(j= jbot; j < jtop; j += chunk(2)-1) {
      jhi= min(jtop, j+chunk(2)-1);
      for(i= ibot; i < itop; i += chunk(1)-1) {
        ihi= min(itop, i+chunk(1)-1);
        /* chk is the number of vertices in this chunk */
        chk= [ihi-i+1, jhi-j+1, khi-k+1];
        /* determine the number of levels in the octree */
        maxdepth= 1;
        num= chk(max)-1;  /* number of cells */
        while(num > 1) {
          num= (num+1)/2;
          maxdepth++;
        }
        start= [i-1, j-1, k-1];  /* C-style */
        res= SliceTree(maxdepth,sizes,chk,start,delta,origin,point,
		       pl_normal,&var2,&arr0);
        if(res) {
          /* have an array of triangles. */
          numTri= arr0.numTri;
          xyzverts= (*(arr0.xyzverts))(,,1:numTri);
          normals= array(pl_normal,3,1:numTri);
          cellIDs= (*(arr0.cellIDs))(1:numTri);
          if(!is_void(var2)) {
            /* need to be able to save one variable value at every
               triangle vertex */
            v2new= (*(arr0.var2))(,1:numTri);
          } else {
            v2new= [];
          }
          lst= TriArrayGrp(
             numTri= numTri,
	     xyzverts= &xyzverts,
	     normals= &normals,
	     var2= &v2new,
	     cellIDs= &cellIDs,
             colors= &nulvar,
	     /* the next three are not used unless assembling strips */
	     triEdg= &nulvar,
	     triStart= &nulvar,
	     nTris= &nulvar,
	     next= &lst
          );
    	}
      }
    }
  }
  return lst;
}

func mak_slice_treecrv(xyz,guard=)
/* DOCUMENT mak_slice_treecrv(xyz,guard=)

     Extract a slicing plane on an ijk 3D grid.

   SEE ALSO: slice_tree, slice_treecrv.
 */
{
  local i, j, k;

  name= "mak_slice_treecrv";
  /* sizes is the number of vertices in the grid. It includes
     the points from the guard cells */
  sizes= dimsof(xyz)(3:5);
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);
  /* assume guard cells by default!!! */
  if(is_void(guard)) guard= 1;
  if(guard) {
    ASSERT, allof(sizes >= 4),
        "The variable array must be nx-by-ny-by-nz in "+name;
    rlsiz= sizes-2;
    ibot= 2; itop= nx-1;
    jbot= 2; jtop= ny-1;
    kbot= 2; ktop= nz-1;
  } else {
    ASSERT, allof(sizes >= 2),
        "The variable array must be nx-by-ny-by-nz in "+name;
    rlsiz= sizes;
    ibot= 1; itop= nx;
    jbot= 1; jtop= ny;
    kbot= 1; ktop= nz;
  }

  /* NOTE: chunk is the number of points in a chunk - the 
     number of cells is one less. A chunk never includes
     either the first or last points in the grid if guard is non-zero. */
  chunk= get_chunk(rlsiz);
  tree= [];
  nulvar= [];
  for(k= kbot; k < ktop; k += chunk(3)-1) {
    khi= min(ktop, k+chunk(3)-1);
    for(j= jbot; j < jtop; j += chunk(2)-1) {
      jhi= min(jtop, j+chunk(2)-1);
      for(i= ibot; i < itop; i += chunk(1)-1) {
        ihi= min(itop, i+chunk(1)-1);
        /* chk is the number of vertices in this chunk */
        chk= [ihi-i+1, jhi-j+1, khi-k+1];
        /* determine the number of levels in the octree */
        maxdepth= 1;
        num= chk(max)-1;  /* number of cells */
        while(num > 1) {
          num= (num+1)/2;
          maxdepth++;
        }
        /* Compute the number of "cells" in each direction for all levels
           of the tree. 
        */
        numx= numy= numz= array(0, maxdepth);
        /* nnx, nny, nnz are numbers of cells */
        nnx= chk(1)-1; nny= chk(2)-1; nnz= chk(3)-1;
        for(id= 1; id <= maxdepth; id++) {
          numx(id)= nnx;
          numy(id)= nny;
          numz(id)= nnz;
          nnx= (nnx+1)/2;
          nny= (nny+1)/2;
          nnz= (nnz+1)/2;
        }
        /* Compute C-style offsets into the range array for the start
           of each level and the total size of the range array.
      	   An element of the range is a struct with a min and a max
           value.
           A single range array contains all levels. The offsets
           are counts of range elements, not counts of the doubles
       	   that make up ranges.
        */
        lens= numx*numy*numz;
        /* NOTE: C-style indices into the range array */
        offsets= lens(psum);
        if(numberof(lens) <= 1) {
          offsets= array(0, 1);
        } else {
          offsets= grow(0, offsets(:-1));
        }
        ranges= array(OctSpan, lens(sum));
        trsiz= array(0, 3, maxdepth);
        trsiz(1,)= numx;
        trsiz(2,)= numy;
        trsiz(3,)= numz;
        tree= OctSTree(
                  maxdepth= maxdepth,
		  size= &sizes,
		  chunk= &chk,
		  start= &[i-1,j-1,k-1], /* C-style start in the full array */
                  trsiz= &trsiz,
                  offsets= &offsets,
                  ranges= &ranges,
		  next= &tree
        	 );
        MakeSliceTreeCrv, xyz, &tree;
      }
    }
  }
  return tree;
}

func slice_treecrv(xyz, point, pl_normal, tree, var2=)
/* DOCUMENT slice_treecrv(xyz, point, pl_normal, tree, var2=)

     Extract a slicing plane on an ijk 3D grid using a pre-computed tree.
     var2 is an optional point-centered variable. If present, 
     it will be interpolated to the points on
     the slicing plane. The resulting variable can be 
     turned into a color and the slice plane can then be smoothly
     colored using one color per vertex instead of one
     color per triangle.

   SEE ALSO: slice_tree, mak_slice_treecrv.
 */
{
  name= "slice_treecrv";
  is_3_vector, point, "point must be a 3 element vector in "+name;
  is_3_vector, pl_normal, "pl_normal must be a 3 element vector in "+name;
  if(!is_void(var2)) {
    dimv2= dimsof(var2);
    dimx= dimsof(xyz);
    ASSERT, (dimv2(1) == 3 && allof(dimx(3:5) == dimv2(2:4)) ),
        "var2 must be nx-by-ny-by-nz in "+name;
    /* force var2 to be a double, using a method that is 
       efficient if it already is a double */
    eq_nocopy,var2,double(var2);
  }

  /* arr0 must be big enough for the worst possible case - 
     i.e. several triangles in each zone.
     The tree was pre-computed using reasonable chunks, so
     allocate a triangle array big enough to handle the largest
     chunk (the final one in the linked list is at least as big
     as any other).
  */
  trptr= &tree;
  while(trptr && *trptr) {
    chunk= *(trptr->chunk);
    trptr= trptr->next;
  }
  numZone= (chunk(1)-1)*(chunk(2)-1)*(chunk(3)-1);
  maxTri= 5*numZone;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= &array(0.0, 3, maxTri);
  } else {
    v2= &[];
  }
  arr0= TriArrayGrp(
          numTri= 0,
	  cellIDs= &array(0, maxTri),
	  xyzverts= &array(0.0, 3, 3, maxTri),
	  normals= &array(0.0, 3, 3, maxTri),
	  var2= v2,
	  triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
         );

  lst= [];
  trptr= &tree;
  while(trptr && *trptr) {
    res= SliceTreeCrv(point,pl_normal,xyz,&var2,&arr0,trptr);
    if(res) {
      /* have an array of triangles. */
      numTri= arr0.numTri;
      xyzverts= (*(arr0.xyzverts))(,,1:numTri);
      normals= array(pl_normal,3,1:numTri);
      cellIDs= (*(arr0.cellIDs))(1:numTri);
      if(!is_void(var2)) {
        /* need to be able to save one variable value at every
           triangle vertex */
        v2new= (*(arr0.var2))(,1:numTri);
      } else {
        v2new= [];
      }
      lst= TriArrayGrp(
             numTri= numTri,
	     xyzverts= &xyzverts,
	     normals= &normals,
             var2= &v2new,
             cellIDs= &cellIDs,
             colors= &nulvar,
             /* the next three are not used unless assembling strips */
             triEdg= &nulvar,
             triStart= &nulvar,
             nTris= &nulvar,
	     next= &lst);
    }
    trptr= trptr->next;
  }
  return lst;
}

func iso3_tree(origin, delta, var, level, colr, tree, var2=)
/* DOCUMENT iso3_tree(origin, delta, var, level, colr, tree, var2=)

     Extract an iso-surface from a variable on a rectangular 3D grid.

   SEE ALSO: mak_isotree, iso3cenreg, iso3_treecrv, iso3_treevarr.
 */
{
  local numTri, nVert, nStrip;

  name= "iso3_tree";
  dimsv= dimsof(var);
  sizes= dimsv(2:4);
  chek_iso_reg, origin, delta, dimsv, var2, name, zcn=0;
  set_color= is_scalar_color(colr, name);
 
  /* arr0 must be big enough for the worst possible case - 
     i.e. several triangles in each zone.
     This can be huge unless the input array is handled in chunks.
     No chunk can be bigger than the last chunk in the list
     (the one that starts at [0,0,0]).
  */
  trptr= &tree;
  while(trptr && *trptr) {
    chunk= *(trptr->chunk);
    trptr= trptr->next;
  }
  numZone= (chunk(1)-1)*(chunk(2)-1)*(chunk(3)-1);
  maxTri= 5*numZone;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= array(0.0, 3, maxTri);
  } else {
    v2= [];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriArrayGrp(
          numTri= 0,
	  cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, 3, maxTri),
          normals= &array(0.0, 3, 3, maxTri),
          var2= &v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
         );

  lst= [];
  trptr= &tree;
  while(trptr && *trptr) {
    if(!is_void(var2)) {
      res= ContourTree2(delta,origin,level,var,&var2,&arr0,trptr);
    } else {
      res= ContourTree(delta,origin,level,var,&arr0,trptr);
    }
    if(res) {
      /* have an array of triangles. */
      numTri= arr0.numTri;
      xyzverts= (*(arr0.xyzverts))(,,1:numTri);
      normals= (*(arr0.normals))(,,1:numTri);
      cellIDs= (*(arr0.cellIDs))(1:numTri);
      if(!is_void(var2)) {
        /* need to be able to save one variable value at every
           triangle vertex */
        v2new= (*(arr0.var2))(,1:numTri);
      } else {
        v2new= [];
      }
      lst= TriArrayGrp(
         numTri= numTri,
         xyzverts= &xyzverts,
         normals= &normals,
         var2= &v2new,
         cellIDs= &cellIDs,
         /* the next three are not used unless assembling strips */
         triEdg= &nulvar,
         triStart= &nulvar,
         nTris= &nulvar,
         next= &lst);
      if(set_color) {
        lst.colors= &float(colr);
      }
    }
    trptr= trptr->next;
  }
  return lst;
}

func iso3_treecrv(xyz, var, level, colr, tree, var2=)
/* DOCUMENT iso3_treecrv(xyz, var, level, colr, tree, var2=)

     Extract an iso-surface from a variable on an ijk 3D grid.
     Uses a tree created by mak_isotree.

   SEE ALSO: iso3_tree, mak_isotree, iso3crv.
 */
{
  extern n_tri_3d, CHKSIZ;
  local numTri, nVert, nStrip;

  name= "iso3_treecrv";
  dimsv= dimsof(var);
  ASSERT, (dimsv(1) == 3 && allof(dimsv(2:4) >= 2) ),
	  "The variable array must be nx-by-ny-by-nz in "+name;
  sizes= dimsv(2:4);
  ASSERT, allof(sizes == *(tree.sizes)),
    "The variable size differs from the size used to build the tree in "+name;
  nx= sizes(1);
  ny= sizes(2);
  nz= sizes(3);
  if(is_void(xyz)) {
    xyz= array(0.0, 3, nx,ny,nz);
    xyz(1,..)= span(-1,1,nx);
    xyz(2,..)= span(-1,1,nx)(-,);
    xyz(3,..)= span(-1,1,nx)(-,-,);
  } else {
    dimsx= dimsof(xyz);
    ASSERT, (dimsx(1) == 4 && dimsx(2) == 3 && allof(dimsx(3:5) == sizes) ),
        "dimensions of xyz do not match those of var in "+name;
  }
  set_color= is_scalar_color(colr, name);
  is_compat(var, var2, "var2 must have the same dimensions as var in "+name);

  /* arr0 must be big enough for the worst possible case - 
     i.e. several triangles in each zone.
	 This can be huge unless the input array is handled in chunks.
	 No chunk can be bigger than the last chunk in the list
	 (the one that starts at [0,0,0]).
  */
  trptr= &tree;
  while(trptr && *trptr) {
    chunk= *(trptr->chunk);
    trptr= trptr->next;
  }
  numZone= (chunk(1)-1)*(chunk(2)-1)*(chunk(3)-1);
  maxTri= 5*numZone;
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= array(0.0, 3, maxTri);
  } else {
    v2= [];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  arr0= TriArrayGrp(
          numTri= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, 3, maxTri),
          normals= &array(0.0, 3, 3, maxTri),
          var2= &v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
         );

  lst= [];
  trptr= &tree;
  while(trptr && *trptr) {
    if(!is_void(var2)) {
      res= ContourTreeCrv2(level,xyz,var,&var2,&arr0,trptr);
    } else {
      res= ContourTreeCrv(level,xyz,var,&arr0,trptr);
    }
    if(res) {
      /* have an array of triangles. */
      numTri= arr0.numTri;
      xyzverts= (*(arr0.xyzverts))(,,1:numTri);
      normals= (*(arr0.normals))(,,1:numTri);
      cellIDs= (*(arr0.cellIDs))(1:numTri);
      if(!is_void(var2)) {
        /* need to be able to save one variable value at every
           triangle vertex */
        v2new= (*(arr0.var2))(,1:numTri);
      } else {
        v2new= [];
      }
      lst= TriArrayGrp(
         numTri= numTri,
         xyzverts= &xyzverts,
         normals= &normals,
         var2= &v2new,
         cellIDs= &cellIDs,
         /* the next three are not used unless assembling strips */
         triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &lst);
      if(set_color) {
        lst.colors= &float(colr);
      }
    }
    trptr= trptr->next;
  }
  return lst;
}

func iso3_treevarr(origin, delta, var, level, colr, tree, var2=)
/* DOCUMENT iso3_treevarr(origin, delta, var, level, colr, tree, var2=)

     Extract an iso-surface from a variable on a rectangular 3D grid.
     The triangles are returned in vertex arrays.

   SEE ALSO: iso3_tree, iso3_treecrv, mak_isotree.
 */
{
  local numTri, numEdg, nVert, nStrip;

  name= "iso3_treevarr";
  dimsv= dimsof(var);
  sizes= dimsv(2:4);
  chek_iso_reg, origin, delta, dimsv, var2, name, zcn=0;
  set_color= is_scalar_color(colr, name);
 
  /* arr0 must be big enough for the worst possible case - 
     i.e. several triangles in each zone.
     This can be huge unless the input array is handled in chunks.
     No chunk can be bigger than the last chunk in the list
     (the one that starts at [0,0,0]).
  */
  trptr= &tree;
  while(trptr && *trptr) {
    chunk= *(trptr->chunk);
    trptr= trptr->next;
  }
  numZone= (chunk(1)-1)*(chunk(2)-1)*(chunk(3)-1);
  maxTri= 5*numZone;
  maxEdg= 3*chunk(1)*chunk(2)*chunk(3);
  if(!is_void(var2)) {
    /* need to be able to save one variable value at every
       triangle vertex */
    v2= array(0.0, 3, maxTri);
  } else {
    v2= [];
  }
  /* note: if second input variable is null, the var2 element of the 
     struct will be zero (the address of a null variable) */
  nulvar= [];
  arr0= TriVertexGrp(
          numTri= 0,
          numEdg= 0,
          cellIDs= &array(0, maxTri),
          xyzverts= &array(0.0, 3, maxEdg),
          normals= &array(0.0, 3, maxEdg),
          ptndx= &array(long, 3*maxTri),
          var2= &v2,
          triEdg= &nulvar,
          triStart= &nulvar,
          nTris= &nulvar,
          next= &nulvar
         );
  /* this scratch array holds the index into the list of vertices
     associated with the point (if any) on each edge in the chunk */
  edgndx= array(-1, 3*chunk(1)*chunk(2)*chunk(3));

  lst= [];
  /* This is a loop over chunks. Each chunk has its own tree. */
  trptr= &tree;
  while(trptr && *trptr) {
    edgndx(*)= -1;
    if(!is_void(var2)) {
      res= ContourTreeVarr2(delta,origin,level,var,&var2,&arr0,trptr,edgndx);
    } else {
      res= ContourTreeVarr(delta,origin,level,var,&arr0,trptr,edgndx);
    }
    if(res) {
      /* have an array of triangles. */
      numTri= arr0.numTri;
      numEdg= arr0.numEdg;
      if(!is_void(var2)) {
        /* need to be able to save one variable value at every
           triangle vertex */
        v2new= (*(arr0.var2))(1:numEdg);
      } else {
        v2new= [];
      }
      xyzverts= (*(arr0.xyzverts))(,1:numEdg);
      normals= (*(arr0.normals))(,1:numEdg);
      cellIDs= (*(arr0.cellIDs))(1:numTri);
      ptndx= (*(arr0.ptndx))(1:3*numTri);
      lst= TriVertexGrp(
        numTri= numTri,
        numEdg= numEdg,
        xyzverts= &xyzverts,
        normals= &normals,
        var2= &v2new,
        cellIDs= &cellIDs,
        ptndx= &ptndx,
        /* the next three are not used unless assembling strips */
        triEdg= &nulvar,
        triStart= &nulvar,
        nTris= &nulvar,
        next= &lst);
      if(set_color) {
        lst.colors= &float(colr);
      }
    }
    trptr= trptr->next;
  }
  return lst;
}
