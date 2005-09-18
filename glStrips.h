/*
 * $Id: glStrips.h,v 1.1.1.1 2005-09-18 22:07:47 dhmunro Exp $
 * Header file for the routines that draw lists of polygons.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __GLSTRIPS__
#define __GLSTRIPS__

#ifdef __cplusplus

extern "C" {

#endif

extern void yglTstrips(long nstrip, long *len, float *xyz, float *norm, 
			   float *colr, long edge, long smooth, long do_light);
extern void yglTstripsAlpha(long nstrip, long *len, float *xyz, float *norm, 
			   float *colr, long edge, long smooth, long do_light);
extern void yglQstrips(long nstrip, long *len, float *xyz, float *norm, 
			   float *colr, long edge, long smooth, long do_light);
extern void yglQstripsAlpha(long nstrip, long *len, float *xyz, float *norm, 
			   float *colr, long edge, long smooth, long do_light);
extern void yglTstripsNdx(long nstrip, long numedg, long ntri, long *len, long *ndx,
				   float *xyz, float *norm, float *colr, long edge);
extern void yglTstripsAlphaNdx(long nstrip, long numedg, long ntri, long *len, 
						long *ndx, float *xyz, float *norm, float *colr, long edge);

void yglTstrips_sm_no_arr(long nstrip, long *len, float *xyz, float *norm, float *colr);
void yglQstrips_sm_no_arr(long nstrip, long *len, float *xyz, float *norm, float *colr);

void yglTstrip_sm_no_arr(long nobj, float *xyz, float *norm, float *colr);
void yglQstrip_sm_no_arr(long nobj, float *xyz, float *norm, float *colr);
void yglTstrip_sm_arr(long nobj, float *xyz, float *norm, float *colr);
void yglQstrip_sm_arr(long nobj, float *xyz, float *norm, float *colr);

void yglTstrip_arr(long nobj, float *xyz, float *norm, float *colr);
void yglQstrip_arr(long nobj, float *xyz, float *norm, float *colr);
void yglTstrip_no_arr(long nobj, float *xyz, float *norm, float *colr);
void yglQstrip_no_arr(long nobj, float *xyz, float *norm, float *colr);

void yglTstrip_arr_no_lite(long nvert, float *xyz, float *colr);
void yglQstrip_arr_no_lite(long nvert, float *xyz, float *colr);
void yglTstrip_no_arr_no_lite(long nvert, float *xyz, float *colr);
void yglQstrip_no_arr_no_lite(long nvert, float *xyz, float *colr);

extern void yglTarrayCubeMapAlpha(long ntri, float *xyz, float *norm, float *colr, long cpervrt);
extern void yglTarrayCubeMap(long ntri, float *xyz, float *norm, float *colr, long cpervrt);

void yglTarray(long smooth, long ntri, float *xyz, float *norm, float *colr, 
			   long edge, long cpervrt, long emit);
void yglTarrayAlpha(long smooth, long ntri, float *xyz, float *norm, 
					float *colr, long edge, long cpervrt, long emit);
void yglTarrayEmit(long do_alpha, long ntri, float *xyz, float *colr, long cpervrt);

void yglQarray(long smooth, long nquad, float *xyz, float *norm, float *colr, 
			   long edge, long cpervrt);
void yglQarrayAlpha(long smooth, long nquad, float *xyz, float *norm, float *colr, 
			   long edge, long cpervrt);

#ifdef __cplusplus
	}
#endif

#endif /* Include/Define */
