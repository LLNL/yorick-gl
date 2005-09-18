/*
 * $Id: gltexsubs.c,v 1.1.1.1 2005-09-18 22:07:52 dhmunro Exp $
 * This file contains functions that deal with version and OS
 * specific features of 3D textures.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include "glcode.h"
#include "glfunc.h"
#include "glBasic.h"
#include "gl3dtex.h"
#include "glcubetex.h"
#include "yio.h"
#include <errno.h>
#include <stdlib.h>

extern void texpal_make_current(void);
extern void yglLdTexPal(long nx, long ny, long nz, unsigned char *tex);

extern int isExtensionSupported(const char *extension);
extern void *LookupFunction(const char *funcName);
extern int TexExtSetup(void);
extern int yglTexExtSetup(void);
static void yglGenCubeTex(void);

/* It is not possible to infer anything about the openGL features
   that will be available at run time from the OpenGL features
   present on the computer where the code is compiled.
   The first reason is the obvious point that the code might
   be run on a different machine than the one that compiled it.
   The other point is that the program might display on a
   remote computer using glx (the features then depend on the
   remote computer's openGL).
   Another problem is that Microsoft has no plans to release a
   version of OpenGL newer than version 1.1. Many Window's OpenGL
   drivers report that they are version 1.2, 1.3, or 1.4,
   so the features of newer versions of OpenGL are often available
   on Windows, they just can't be accessed through Microsoft's
   OpenGL library.
   On Windows, OpenGL extensions are accessed via function pointers
   acquired at run time. There is thus no point trying to decide at
   compile time what features are present. The best approach seems to
   be to define constants for any interesting extensions that aren't
   present in <GL/gl.h> and query at runtime.
   In the case of Unix systems, extensions are accessed through entries
   must be present in the OpenGL library or the machine where the code
   is running. On Unix systems it is not possible to call an
   extension unless it is available at compile time. However, it is
   still necessary to check for the presence of the extension at run time.
   
   The most reasonable approach to using funcitonality newer than OpenGL 1.1
   seems to be to only use extensions. On Windows machines, all attempts to
   check for the extension are deferred until runtime. On Unix systems,
   the extension must be present on the machine that is compiling the code
   and on the machine displaying the output during the run.
   The only reasonably portable approach is to access all OpenGL features
   added after version 1.1 as extensions. One problem with this approach
   is that the headers on a Windows computer will not define any
   useful extensions. The solution is to hand code all the constants
   needed to query for the extension.

   As near as I can tell, it is not possible to build a code on Linux
   that can use both an ATI and an Nvidia extension (the code will
   link against either a library supplied by ATI or a library supplied
   by Nvidia, but not both). This argues pretty strongly against
   using vendor specific extensions.

   To keep the messing about with extensions localized, I will define
   my own function for each type of functionality I provided by
   extensions that I want to use. If the functionality is not present,
   the function will simply return.
*/

/* Set a flag indicating whether the computer on which yorgl was compiled
   has 3D textures available, which means OpenGL 1.2 or later.
   Whether this function should be called at run time depends 
   on whether the OpenGL on the "desktop" computer is version 1.2
   or newer. An example of where this might fail is a program running
   on a Linux computer and displaying back to a Windows computer
   via glx. Windows only has 3D textures as an extension. */

#ifdef GL_TEXTURE_3D
int host_has_3dtex= 1;
#else
int host_has_3dtex= 0;
#endif


void myglTexImage3D(GLenum target, GLint level, GLint internalformat,
                    GLsizei width, GLsizei height, GLsizei depth,
                    GLint border, GLenum format, GLenum type,
                    const GLvoid* pixels)
{
  /* If earlier tests indicate that both client and server have 3D textures
     built-in, call that function. This test ensures this branch will
     not be taken on any machine where the ifdef causes the function
     call to be absent */
  if(glCurrWin3d->hasTex3d) {
#ifdef GL_TEXTURE_3D
    glTexImage3D(target, level, internalformat, width, height, depth,
               border, format, type, pixels);
#endif
    return;
  }
  /* This code the 3D texture extension on either Windows or a Unix system.
     It returns doing nothing if the 3D texture extension isn't available.
  */
  if(!glCurrWin3d->hasTex3dExt) return;
#ifdef WIN32
  (*(glCurrWin3d->myglTexImage3D_ptr)) (target, level, internalformat, width, height, depth,
           border, format, type, pixels);
#else
  glTexImage3DEXT(target, level, internalformat, width, height, depth,
               border, format, type, pixels);
#endif
}

int yglQueryTex3d(glWinProp *win)
{
  char *version, msg[100];
  double vernum;
  int res;
  GLenum err;
  const GLubyte *errstr;
  char *rest;

  /* clear earlier GL errors */
  err= glGetError();
  errstr= gluErrorString(err);
  /* the 3D window may not have been created yet */
  if(!win) {
    glWinProp *newWin= 0;  /* zero forces window creation */
    yglPrepDraw(newWin);
    win= glCurrWin3d;
    if(!win) return 0;
  }
  if(win->hasTex3d || win->hasTex3dExt) return 1;
  if(win->tex3dChecked) return 0;  /* previous check showed no 3D text. */
  /* Make sure the OpenGL context is current */
  yglMakeCurrent(win);
  version = (char*) glGetString(GL_VERSION);
  DEMAND(version, "Failed to get OpenGL version number")
  sprintf(msg, "OpenGL version number is %s\n", version);
  YputsOut(msg);
  vernum= strtod(version, &rest);
  glCurrWin3d->tex3dChecked= 1;  /* note check for 3D texture occurred */

  if(host_has_3dtex && vernum > 1.199) {
    /* This code was compiled on a computer the OpenGL 1.2 or later
       and the display is on a computer with OpenGL 1.2 or later.
    */
    win->hasTex3d= 1;   /* has 3D texture support */
#ifdef GL_TEXTURE_3D
    win->myGL_TEXTURE_3D= GL_TEXTURE_3D;
    win->myGL_PROXY_TEXTURE_3D= GL_PROXY_TEXTURE_3D;
#endif
    return 1;
  }
  
  res = isExtensionSupported("GL_EXT_texture3D");
  if(res) {
    win->hasTex3dExt= 1;   /* has 3D texture extension */
    win->myGL_TEXTURE_3D= GL_TEXTURE_3D_EXT;
    win->myGL_PROXY_TEXTURE_3D= GL_PROXY_TEXTURE_3D_EXT;
#ifdef WIN32
    win->myglTexImage3D_ptr = (MYTEX3DFUNC) LookupFunction("glTexImage3DEXT");
#endif
    return 1;
  } else {
    return 0;   /* no 3D texture support */
  }
}

void *LookupFunction(const char *funcName)
{
  /* The LookupFunction function finds function pointers 
     in dynamic libraries. */
#ifdef WIN32
  return wglGetProcAddress(funcName);
#elif defined(IRIX)
  void *libHandle = dlopen("libgl.so", RTLD_LAZY);
  void *func = dlsym(libHandle, funcName);
  dlclose(libHandle);
  return func;
#else
  /* other OSes... */
  return 0;
#endif
}
