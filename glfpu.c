/* glfpu.c
 * Functions to turn on/off FPE interrupt masks before and after OpenGL
 * calls.  The OpenGL standard does not permit FPE trapping required by
 * yorick (yet doesn't bother turning it on and off itself like libm).
 * This hasn't been a problem on any platform until Mac OSX 10.7, which
 * somehow manages to generate SIGFPE while attempting to clear the newly
 * created GLX window.  It also apparently has installed its own SIGFPE
 * handler, since yorick's is never invoked, and drops into an infinite
 * loop.  This code uses fenv.h to restore the default FPU mode on entry
 * to any routine which calls an OpenGL function, and put back yorick's
 * FPU mode upon return.  The code is fragile, since it may not restore
 * the original mode if the code is interrupted (for example by C-c).
 * A separate interpreted API has been added to restore the yorick FPU
 * environment in case this happens.
 */
/* Copyright (c) 2012, David H. Munro.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#ifdef YGL_CONFIG_MAIN
extern void ygl_fpemask(int on);
int main(int argc, char *argv[])
{
  ygl_fpemask(1);
  return 0;
}

#else
#include "glfunc.h"
#endif

#ifndef MISSING_FENV_H

#include "fenv.h"

static fenv_t ygl_fenv;
static int ygl_valid_fenv = 0;

void
ygl_fpemask(int on)
{
  ygl_valid_fenv = (ygl_valid_fenv ||
                    !fegetenv(&ygl_fenv));
  if (ygl_valid_fenv) {
    if (on)
      fesetenv(&ygl_fenv);
    else
      fesetenv(FE_DFL_ENV);
  }
}

#else

void
ygl_fpemask(int on)
{
}

#endif
