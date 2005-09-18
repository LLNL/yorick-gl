/*
 * $Id: glMouse.h,v 1.1 2005-09-18 22:07:46 dhmunro Exp $
 * Header file for the routines that handle mouse motion
 * in OpenGL windows.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __GLMOUSE__
#define __GLMOUSE__

#include "glcode.h"

#ifdef __cplusplus

extern "C" {

#endif

typedef struct mspoint {
  int x;
  int y;
} MsPoint;

#define ACTION_ROTATE 1
#define ACTION_ZOOM 2
#define ACTION_PAN 3

extern void DoLButtonDown(int flags, int x, int y);
extern void DoLButtonUp(int flags, int x, int y, glWinProp *theWin3d);
extern void new_mouse_pos(int flags, int x, int y, int force, glWinProp *theWin3d);
extern void yglMouseRot(int x, int y);
extern void yglMousePan(int x, int y);
extern void yglMouseZoom(int x, int y);
extern void yglMouseMove(int flags, int x, int y, glWinProp *theWin3d); 

extern int ygl_use_mouse_move;
extern double ygl_ms_mov_val;

#ifdef __cplusplus
	}
#endif

#endif /* Include/Define */
