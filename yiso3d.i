/*
 * $Id: yiso3d.i,v 1.1 2005-09-18 22:07:52 dhmunro Exp $
 * autoloads for yiso3d subpackage of yorgl package
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

autoload, "contour.i", iso3, iso3_tree, iso3_treecrv, iso3_treevarr;
autoload, "contour.i", iso3cencrv, iso3cencrvndx, iso3cenreg, iso3cenregndx;
autoload, "contour.i", iso3cenregngrd, iso3hex, iso3ndx, iso3zcencrv;
autoload, "contour.i", iso3zcencrvndx, iso3zcenreg, iso3zcenregndx;
autoload, "contour.i", iso3zcenregngrd, mak_isotree, mak_slice_treecrv;
autoload, "contour.i", slice_tree, slice_treecrv;
