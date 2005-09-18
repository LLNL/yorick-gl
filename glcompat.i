/*
 * $Id: glcompat.i,v 1.1 2005-09-18 22:07:49 dhmunro Exp $
 * This file contains functions from earlier versions of
 * yorgl that have been replaced or renamed. They are provided
 * to allow earlier decks to run with newer versions of yorgl.
 * Please do not rely on these functions for the long term -
 * they may break at some point in the future.  
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

func draw3(called_as_idler)
{
  draw3d,called_as_idler;
}

func draw3_trigger
{
  draw3d_trigger;
}

func clear3(void)
{
  clear3d;
}

func clear3_direct(void)
{
  clear3d_direct;
}

func clear3_cache(void)
{
  clear3d_cache;
}

func gl_getbackrgb(void)
{
  return getbackrgb3d();
}

func gl_getcage_rgb(void)
{
  return getcage_rgb3d();
}

func gl_getgrid_rgb(void)
{
  return getgrid_rgb3d();
}

func light3(ambient=,diffuse=,specular=,spower=,sdir=)
{
  light3(ambient=ambient,diffuse=diffuse,specular=specular,
         spower=spower,sdir=sdir);
}

func stdview3(dummy)
{
  stdview3;
}

func prtview3
{
  prtview3;
}

func lookat(eye,center,up)
{
  return lookat3d(eye,center,up);
}

func get3lims(dummy)
{
  return get3lims(dummy);
}

func get3_normal(xyz, nxyz)
{
  return get_normal3d(xyz, nxyz);
}

func get3_centroid(xyz, nxyz)
{
  return get_centroid3d(xyz, nxyz);
}

func pl3dpoly_gl(nv, xyz, color, norm, draw_edge=)
{
  return plpoly3d(nv, xyz, color, norm, draw_edge=draw_edge);
}

func pl3dcell_gl(corners, color)
{
  return plcell3d(corners, color);
}

func plm3d_gl(xyz, color)
{
  return plm3d(xyz, color);
}

func plf3d_gl(xyz, color)
{
  return plf3d(xyz, color);
}

func surf3d_gl(xyz, norm, color, flip=)
{
  return plsurf3d(xyz, norm, color, flip=flip);
}

func colrsurf3d_gl(xyz, norm, color, flip=)
{
  return plcolrsurf3d(xyz, norm, color, flip=flip);
}

func pl3dlines_gl(xyz, color)
{
  return pllines3d(xyz, color);
}

func pl3dpoints_gl(xyz, color)
{
  return plpoints3d(xyz, color);
}

func pl3dglyphs_gl(xyz, scal, theta, phi, color)
{
  return plglyphs3d(xyz, scal, theta, phi, color);
}

func pl3dtrilists(tris,flip=,offset=,cubemap=,emit=)
{
  return pltrilists3d(tris,flip=flip,offset=offset,cubemap=cubemap,emit=emit);
}

func pl3dtstrips_gl(nv, xyz, color, norm, draw_edge=)
{
  return pltstrips3d(nv, xyz, color, norm, draw_edge=draw_edge);
}

func pl3dqstrips_gl(nv, xyz, color, norm, draw_edge=)
{
  return plqstrips3d(nv, xyz, color, norm, draw_edge=draw_edge);
}

func pl3dtivstrips_gl(nv, ndx, xyz, norm, color, draw_edge=)
{
  return pltivstrips3d(nv, ndx, xyz, norm, color, draw_edge=draw_edge);
}

func pl3dtarray_gl(xyz, norm, color, ntri, cubemap, emit, draw_edge=)
{
  return pltarray3d(xyz, norm, color, ntri, cubemap, emit, draw_edge=draw_edge);
}

func pl3dqarray_gl(xyz, norm, color, nquad, draw_edge=)
{
  return plqarray3d(xyz, norm, color, nquad, draw_edge=);
}

func pl3dtivarray_gl(ptndx, xyz, norm, color, ntri, nvert)
{
  return pltivarray3d(ptndx, xyz, norm, color, ntri, nvert);
}

func putpix_gl(arr)
{
  return putpix3d(arr);
}

func pl3dtexcell2_gl(delta, texval)
{
  return pltex2dvol(delta, texval);
}

func pltex3d_gl(nslab, boxsiz, texval, origin=)
{
  return pltex3dvol(nslab, boxsiz, texval, origin=origin);
}
