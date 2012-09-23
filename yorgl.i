/*
 * $Id: yorgl.i,v 1.1 2005-09-18 22:08:03 dhmunro Exp $
 * autoloads for yorgl package
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

autoload, "dlist3d.i", CollapseTri, SortTri, clear3d, clear3d_cache;
autoload, "dlist3d.i", clear3d_direct, draw3d, get_centroid3d, get_lims3d;
autoload, "dlist3d.i", get_normal3d, getpix3d, gl_rr, light3d, lookat3d;
autoload, "dlist3d.i", palette3d, plcell3d, plcolrsurf3d, plf3d, plglyphs3d;
autoload, "dlist3d.i", pllines3d, plm3d, plpoints3d, plpoly3d, plqarray3d;
autoload, "dlist3d.i", plqstrips3d, plsurf3d, pltarray3d, pltex2dvol;
autoload, "dlist3d.i", pltex3dvol, pltivarray3d, pltivstrips3d, pltrilists3d;
autoload, "dlist3d.i", pltstrips3d, pltvarray3d, prtview3d, putpix3d;
autoload, "dlist3d.i", stdview3d;

autoload, "slicenew.i", nuslice2, nuslice2x, slice2_precision, slice2only;
