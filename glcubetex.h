/*
 * $Id: glcubetex.h,v 1.1 2005-09-18 22:07:49 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#ifndef __H_3DCUBETEX__
#define __H_3DCUBETEX__

#ifdef __cplusplus
extern "C" {
#endif

/* In case your <GL/gl.h> does not advertise EXT_texture_cube_map... 
   NOTE: there are no new functions needed for cube maps */

#ifndef GL_EXT_texture_cube_map

# define GL_NORMAL_MAP_EXT                   0x8511
# define GL_REFLECTION_MAP_EXT               0x8512
# define GL_TEXTURE_CUBE_MAP_EXT             0x8513
# define GL_TEXTURE_BINDING_CUBE_MAP_EXT     0x8514
# define GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT  0x8515
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_X_EXT  0x8516
# define GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT  0x8517
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT  0x8518
# define GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT  0x8519
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT  0x851A
# define GL_PROXY_TEXTURE_CUBE_MAP_EXT       0x851B
# define GL_MAX_CUBE_MAP_TEXTURE_SIZE_EXT    0x851C

#endif

#ifdef __cplusplus
}
#endif

#endif /* __H_3DCUBETEX__ */
