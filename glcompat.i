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

draw3 = draw3d;
draw3_trigger = draw3d_trigger;
clear3 = clear3d;
clear3_direct = clear3d_direct;
clear3_cache = clear3d_cache;
gl_getbackrgb = getbackrgb3d;
gl_getcage_rgb = getcage_rgb3d;
gl_getgrid_rgb = getgrid_rgb3d;
light3 = light3d;
stdview3 = stdview3d;
prtview3 = prtview3d;
lookat = lookat3d;
get3lims = get_lims3d;
get3_normal = get_normal3d;
get3_centroid = get_centroid3d;
pl3dpoly_gl = plpoly3d;
pl3dcell_gl = plcell3d;
plm3d_gl = plm3d;
plf3d_gl = plf3d;
surf3d_gl = plsurf3d;
colrsurf3d_gl = plcolrsurf3d;
pl3dlines_gl = pllines3d;
pl3dpoints_gl = plpoints3d;
pl3dglyphs_gl = plglyphs3d;
pl3dtrilists = pltrilists3d;
pl3dtstrips_gl = pltstrips3d;
pl3dqstrips_gl = plqstrips3d;
pl3dtivstrips_gl = pltivstrips3d;
pl3dtarray_gl = pltarray3d;
pl3dqarray_gl = plqarray3d;
pl3dtivarray_gl = pltivarray3d;
putpix_gl = putpix3d;
pl3dtexcell2_gl = pltex2dvol;
pltex3d_gl = pltex3dvol;
