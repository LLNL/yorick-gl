/*
 * $Id: playgl.h,v 1.1.1.1 2005-09-18 22:08:03 dhmunro Exp $
 * OpenGL portability layer interface declarations
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

/*************************************************************************/
/*                  WARNING
 * PLUG_EXPORT below is correct only for the package which is
 * building oglx.c/oglw.c, which defines these functions.
 * Another package should use PLUG_API.  Since there are no
 * other packages currently, I hardwired the changes here.
 * The correct procedure would be to define a PLAY_GL_INTERNAL
 * macro in all the source files which include playgl.h in this
 * directory.  Since playgl.h is not installed, the chance of
 * misuse for the time being is nil.
 */
/*************************************************************************/

/* include play.h, then playgl.h, then GL/gl.h */

#include "plugin.h"

typedef struct p_glwin p_glwin;   /* opaque to platform independent code */

BEGIN_EXTERN_C

PLUG_EXPORT p_glwin *p_glcreate(p_win *parent, int width, int height,
                             int x, int y, void *ctx);
PLUG_EXPORT void p_gldestroy(p_glwin *w);
PLUG_EXPORT void p_glresize(p_glwin *w, int width, int height, int x, int y);
PLUG_EXPORT void p_glswap(p_glwin *w);
PLUG_EXPORT void p_glcurrent(p_glwin *w);

END_EXTERN_C

/* Mr. Bill buggered <GL/gl.h> so it won't work without this
 * -- better to have it confined to this one place */
#ifdef _WIN32
#include <windows.h>
#endif
