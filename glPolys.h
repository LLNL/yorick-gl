/*
 * $Id: glPolys.h,v 1.1.1.1 2005-09-18 22:07:46 dhmunro Exp $
 * Header file for the routines that draw lists of polygons.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __GLPOLYS__
#define __GLPOLYS__

#ifdef __cplusplus

extern "C" {

#endif

void gl_polys_arr(long npoly, long *nverts, float *xyz, float *norm, 
			  float *colr);
void gl_polys_no_arr(long npoly, long *nverts, float *xyz, float *norm, 
			  float *colr);
void gl_polys_arr_no_lite(long npoly, long *nverts, float *xyz, 
			  float *colr);
void gl_polys_no_arr_no_lite(long npoly, long *nverts, float *xyz, 
			  float *colr);
void gl_polys_sm_arr(long npoly, long *nverts, float *xyz, float *norm, 
			  float *colr);
void gl_polys_sm_no_arr(long npoly, long *nverts, float *xyz, float *norm, 
			  float *colr);

#ifdef __cplusplus
	}
#endif

#endif /* Include/Define */
