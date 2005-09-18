/*
 * $Id: glviewpoint.h,v 1.1.1.1 2005-09-18 22:07:53 dhmunro Exp $
 * Declarations for the viewpoint transformation..
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __GLVIEWPOINT__
#define __GLVIEWPOINT__

#ifdef __cplusplus

extern "C" {

#endif

extern float ygl_fov;
extern double ygl_eye[3], ygl_center[3], ygl_up[3], ygl_view[3], ygl_viewdist;

#ifdef __cplusplus
	}
#endif

#endif /* Include/Define */
