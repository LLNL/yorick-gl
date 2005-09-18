/*
 * $Id: gl3dtex.h,v 1.1 2005-09-18 22:07:45 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#ifndef __H_3DTEX__
#define __H_3DTEX__

#ifdef __cplusplus
extern "C" {
#endif


/*  Extension: GL_EXT_texture3D */

/* OpenGl 1.2 and greater has 3D textures standard. Version 1.1 may
   have it as an extension. There may be no 3D textures.
   It is hard to write code that covers all cases.
   The strategy is to include the GL headers in glcode.h and assume
   it will define the appropriate 3D texture constants.
   Macros are defined based on which, if any, 3D texture constants 
   are defined. All code accesses 3D texture features through
   macros I define.
   WARNING: Linux with an Nvidia driver does not currently replace
   the OpenGL headers - it puts them in /usr/share/doc/NVIDIA or some 
   such place. To use 3D textures, modify the include path to pick
   up the NVIDIA version.
   WARNING: Windows apparently will never have an OpenGL library or
   headers beyond version 1.1. many OpenGL drivers on Windows are
   version 1.2 or above. The proper appraoch is to try and access
   advanced GL features as extensions. To do this, we need to
   augment gl.h with additional constants.
*/
#ifdef WIN32
#define GL_TEXTURE_BINDING_3D_EXT       0x806A
#define GL_PACK_SKIP_IMAGES_EXT         0x806B
#define GL_PACK_IMAGE_HEIGHT_EXT        0x806C
#define GL_UNPACK_SKIP_IMAGES_EXT       0x806D
#define GL_UNPACK_IMAGE_HEIGHT_EXT      0x806E
#define GL_TEXTURE_3D_EXT               0x806F
/* #define GL_TEXTURE_3D                   0x806F */
#define GL_PROXY_TEXTURE_3D_EXT         0x8070
/* #define GL_PROXY_TEXTURE_3D             0x8070 */
#define GL_TEXTURE_DEPTH_EXT            0x8071
#define GL_TEXTURE_WRAP_R_EXT           0x8072
#define GL_MAX_3D_TEXTURE_SIZE_EXT      0x8073
#endif

#ifndef WIN32
#ifndef GL_TEXTURE_3D
#ifndef GL_TEXTURE_3D_EXT
#define GL_TEXTURE_DEPTH_EXT  0x8071
#define GL_TEXTURE_WRAP_R_EXT 0x8072
#endif
#endif
#endif

extern void myglTexImage3D(GLenum target, GLint level, GLint internalformat,
                           GLsizei width, GLsizei height, GLsizei depth,
                           GLint border, GLenum format, GLenum type,
                           const GLvoid* pixels);

#ifdef __cplusplus
}
#endif

#endif /* __H_3DTEX__ */
