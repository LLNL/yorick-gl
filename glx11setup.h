/*
 * $Id: glx11setup.h,v 1.1 2005-09-18 22:07:53 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifndef GL_X11SETUP_H
#define GL_X11SETUP_H

#include <X11/Xlib.h>
/* XVisualInfo declared in Xutil.h */
#include <X11/Xutil.h>

typedef void (*GxHandler)(Engine *, Drawing *, XEvent *);

#endif
