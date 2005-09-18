/*
 * $Id: glBasic.h,v 1.1.1.1 2005-09-18 22:07:45 dhmunro Exp $
 * Header file for window creation and drawing related functions
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef __GLBASIC__
#define __GLBASIC__

#include "glcode.h"

#ifdef __cplusplus

extern "C" {

#endif

extern void yglDoGetPixels(long nx, long ny, unsigned char *pix);
extern glWinProp *yglMakWin(char *displayName, int width, int height, char *title);
extern int  yglWin3d(int num, int w, int h);
extern int  yglWinKill3d(int num);
extern int shutdown3d(glWinProp *win3d);
extern void resetcurrwin3d(void);
extern int winnum3d(glWinProp *theWin3d);
extern void yglDoDraw(glWinProp *theWin3d);
extern int isWin3d(glWinProp *win);
extern void yglForceDraw(void);
extern void yglPrepDraw(glWinProp *theWin3d);
extern void yglFinFrame(void);
extern void yglFinDirect(void);
extern void yglFinCache(void);
extern void yglMakeCurrent(glWinProp *theWin3d);
extern void yglPrepList(void);

extern glWinProp *glWin3dList[8];

#ifdef __cplusplus
	}
#endif

#endif /* Include/Define */
