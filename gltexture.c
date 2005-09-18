/*
 * $Id: gltexture.c,v 1.1.1.1 2005-09-18 22:07:53 dhmunro Exp $
 * This file contains functions that use OpenGL texture capabilities.
 * A significant application is interactive volume visualization
 * using either 2D or 3D textures (prefereable 3D).
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include "glcode.h"
#include "glfunc.h"
#include "glMouse.h"
#include "glBasic.h"
#include "gl3dtex.h"
#include "glcubetex.h"
#include "Contour3D.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef DEBUG_PRINT
#include "play.h"
#include <stdio.h>
#endif

extern int ycInitCartPcen(long sizes[3], long offsets[6], 
            double deltas[3], double origin[3], 
            double *var, double *var2);
extern int ycContourTet_array(long make_strip, long sizes[3], 
            double level, double *var, double *var2, yPoint3D *xyz,
            yPoint3D *grd, unsigned char *flag,
            TriArrayGrp *triangles);

extern int isExtensionSupported(const char *extension);
extern void *LookupFunction(const char *funcName);
extern int TexExtSetup(void);
static void slice_box(double s, double sbox[8], double *origin, 
                      double *len, TriArrayGrp *tris);
extern int yglTexExtSetup(void);
static void yglGenCubeTex(void);

static GLuint texName= 0;
static GLuint texName3d= 0;
static unsigned char *texImage3d= 0;
static TriArrayGrp *tris_tex3d= 0;
static float fracx3d, fracy3d, fracz3d;
static long  nx3d, ny3d, nz3d;

static GLuint texNameCube= 0;
#define ygl_cubemap_size 64
/* 6 square faces of RGBA texels */
static unsigned char cubeMaps[6][ygl_cubemap_size][ygl_cubemap_size][4];

static GLenum faceTarget[6] = {
  GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT,
  GL_TEXTURE_CUBE_MAP_NEGATIVE_X_EXT,
  GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT,
  GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT,
  GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT,
  GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT
};

/* J. F. Blinn and M. E. Newell. Texture and reflection in computer generated images. 
     Communications of the ACM, 19(10):542-547, 1976. 
*/

int yglQueryTexCube(void)
{
  int res;

  if(glCurrWin3d->hascubetex >= 0) return glCurrWin3d->hascubetex;
  /* Make sure the OpenGL context is current */
  yglMakeCurrent(glCurrWin3d);
  res = isExtensionSupported("GL_EXT_texture_cube_map");
  if(res) {
    glCurrWin3d->hascubetex= 1;   /* has cube-map texture support */
  } else {
    glCurrWin3d->hascubetex= 0;   /* not supported */
  }
  return glCurrWin3d->hascubetex;
}

double yglGetVers3d(void)
{
  /* Return the OpenGL version number. */
  char *version;
  double glver;

  yglMakeCurrent(glCurrWin3d);
  version = (char*) glGetString(GL_VERSION);
  glver= atof(version);
  return glver;
}

int yglTexExtSetup(void)
{
  int res;

  yglMakeCurrent(glCurrWin3d);

  if(glCurrWin3d->hasTexExt >= 0) return glCurrWin3d->hasTexExt;
  res = isExtensionSupported("GL_EXT_texture");
  res= 1;
  if(res) {
    glCurrWin3d->hasTexExt= 1;
    glCurrWin3d->myglBindTexture3D = (glBindTextureProc) LookupFunction("glBindTexture3DEXT");
  } else {
    glCurrWin3d->hasTexExt= 0;
    glCurrWin3d->myglBindTexture3D = NULL;
  }
  return glCurrWin3d->hasTexExt;
}

int isExtensionSupported(const char *extension)
{
  /* This function is based on the example in Mark Kilgard's
     article about openGL extensions.
     WARNING!!! an OpenGL context must be established before
     calling this function!!! */
  const GLubyte *extensions= NULL;
  const GLubyte *start;
  GLubyte *where, *terminator;
  const GLubyte *errstr;
  GLenum err;

  /* Extension names should not have spaces,
	   so return "not supported" if the name is malformed. */
  where= (GLubyte *) strchr(extension, ' ');
  if(where || *extension == '\0') return 0;

  extensions= glGetString(GL_EXTENSIONS);
  err= glGetError();
  errstr= gluErrorString(err);
  if(!extensions) return 0;

	/* Parsing OpenGL extension strings is tricky. 
	   handle sub-strings etc. */
  start= extensions;
  for( ; ; ) {
    where= (GLubyte *) strstr((const char *) start, extension);
    if(!where) break;
    terminator= where+strlen(extension);
    if(where == start || *(where-1) == ' ') {
      if(*terminator == ' ' || *terminator == '\0') {
        return 1;
      }
    }
    start= terminator;
  }
  return 0;
}

#define USE_DATA

void yglTexcells(long nx, long ny, long nz, double delta[3], 
                 unsigned char *data, unsigned char *ctab)
{
  /* The input is a 3D array of characters, two corners,
     and a table that maps a character to RGBA.
  */
  long i, j, k, n, nnx, nny, nnz, nh, nv, nhs, nvs, indt, dt, dir;
  float fracx, fracy, fracz;
  float x0, y0, z0, x1, y1, z1, dx, dy, dz;
  float x, y, z;
  double absx, absy, absz;
  unsigned char *newTexImage, *subImage;
#ifdef DO_TEX_TEST
  GLint format, redSize, greenSize, blueSize, alphaSize;
  GLint texWidth, texHeight;
#endif

  yglPrepTex2d();

  /* Textures must be a power of 2 in size. Find the smallest 
     power of 2 that will hold the input array */
  n= nx-1;
  nnx= 1;
  if(n) {
    do {
      nnx *= 2;
    } while(n= n/2);
  }
  fracx= (float)nx/(float)nnx;
  n= ny-1;
  nny= 1;
  if(n) {
    do {
      nny *= 2;
    } while(n= n/2);
  }
  fracy= (float)ny/(float)nny;
  n= nz-1;
  nnz= 1;
  if(n) {
    do {
      nnz *= 2;
    } while(n= n/2);
  }
  fracz= (float)nz/(float)nnz;

  /* There are nx by ny by nz cells to be volume rendered.
     FOR NOW, assume that the cells are aligned with the
     axes. */
  dx= (float) delta[0];
  dy= (float) delta[1];
  dz= (float) delta[2];

#ifdef USE_DATA
  /* find the axis most closely matching the viewing direction */
  absx= glCurrWin3d->view[0] < 0 ? -glCurrWin3d->view[0] : glCurrWin3d->view[0];
  absy= glCurrWin3d->view[1] < 0 ? -glCurrWin3d->view[1] : glCurrWin3d->view[1];
  absz= glCurrWin3d->view[2] < 0 ? -glCurrWin3d->view[2] : glCurrWin3d->view[2];
  if(absx < absy) {
    if(absy < absz) {
      /* the view direction is most nearly along z */
      dir= 3;
      if(glCurrWin3d->view[2] < 0) dir= -3;
      nh= nnx;
      nv= nny;
      nhs= nx;
      nvs= ny;
    } else {
      /* the view direction is most nearly along y */
      dir= 2;
      if(glCurrWin3d->view[1] < 0) dir= -2;
      nh= nnz;
      nv= nnx;
      nhs= nz;
      nvs= nx;
    }
  } else {
    if(absx < absz) {
      /* the view direction is most nearly along z */
      dir= 3;
      if(glCurrWin3d->view[2] < 0) dir= -3;
      nh= nnx;
      nv= nny;
      nhs= nx;
      nvs= ny;
    } else {
      /* the view direction is most nearly along x */
      dir= 1;
      if(glCurrWin3d->view[0] < 0) dir= -1;
      nh= nny;
      nv= nnz;
      nhs= ny;
      nvs= nz;
    }
  }
#else
  nh= nnx;
  nv= nny;
#endif
  newTexImage= (unsigned char *)malloc(sizeof(unsigned char)*4*nh*nv);
  subImage= (unsigned char *)malloc(sizeof(unsigned char)*4*nhs*nvs);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, nh, nv, 0, GL_RGBA, 
               GL_UNSIGNED_BYTE, newTexImage);

#ifdef DO_TEX_TEST
  glTexImage2D(GL_PROXY_TEXTURE_2D, 0, GL_RGBA, nnx, nny, 0, GL_RGBA, 
               GL_UNSIGNED_BYTE, NULL);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_RED_SIZE, &redSize);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_GREEN_SIZE, &greenSize);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_BLUE_SIZE, &blueSize);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_ALPHA_SIZE, &alphaSize);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_WIDTH, &texWidth);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_HEIGHT, &texHeight);
#endif

  glGenTextures(1, &texName);
  glBindTexture(GL_TEXTURE_2D, texName);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  /*  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); */

  glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

#ifdef USE_DATA
  x0= 0.0;
  x1= x0+(nx-1)*dx;
  y0= 0.0;
  y1= y0+(ny-1)*dy;
  z0= 0.0;
  z1= z0+(nz-1)*dz;

  switch(dir) {
  case 1:
    /* the view direction is most nearly along x */
    x= x0;
    /* draw planes in back to front order */
    for(i= 0; i < nx; i++) {
      /* load the proper texture plane */
      for(k= 0; k < nz; k++) {
        for(j= 0; j < ny; j++) {
          indt= 4*(j+k*ny);
          dt= 4*data[i+j*nx+k*nx*ny];
          subImage[indt  ]= ctab[dt];
          subImage[indt+1]= ctab[dt+1];
          subImage[indt+2]= ctab[dt+2];
          subImage[indt+3]= ctab[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ny, nz, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x, y0, z0);
      glTexCoord2f(fracy, 0.0f);
      glVertex3f(x, y1, z0);
      glTexCoord2f(fracy, fracz);
      glVertex3f(x, y1, z1);
      glTexCoord2f(0.0f, fracz);
      glVertex3f(x, y0, z1);
      glEnd();
      x += dx;
    }
    break;
  case -1:
    /* the view direction is most nearly along negative x */
    x= x1;
    /* draw planes in back to front order */
    for(i= nx-1; i >= 0; i--) {
      /* load the proper texture plane */
      for(k= 0; k < nz; k++) {
        for(j= 0; j < ny; j++) {
          indt= 4*(j+k*ny);
          dt= 4*data[i+j*nx+k*nx*ny];
          subImage[indt  ]= ctab[dt];
          subImage[indt+1]= ctab[dt+1];
          subImage[indt+2]= ctab[dt+2];
          subImage[indt+3]= ctab[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ny, nz, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x, y0, z0);
      glTexCoord2f(fracy, 0.0f);
      glVertex3f(x, y1, z0);
      glTexCoord2f(fracy, fracz);
      glVertex3f(x, y1, z1);
      glTexCoord2f(0.0f, fracz);
      glVertex3f(x, y0, z1);
      glEnd();
      x -= dx;
    }
    break;
  case 2:
    /* the view direction is most nearly along y */
    y= y0;
    /* draw planes in back to front order */
    for(j= 0; j < ny; j++) {
      /* load the poper texture plane */
      for(i= 0; i < nx; i++) {
        for(k= 0; k < nz; k++) {
          indt= 4*(k+i*nz);
          dt= 4*data[i+j*nx+k*nx*ny];
          subImage[indt  ]= ctab[dt];
          subImage[indt+1]= ctab[dt+1];
          subImage[indt+2]= ctab[dt+2];
          subImage[indt+3]= ctab[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nz, nx, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x0, y, z0);
      glTexCoord2f(0.0f, fracz);
      glVertex3f(x0, y, z1);
      glTexCoord2f(fracx, fracz);
      glVertex3f(x1, y, z1);
      glTexCoord2f(fracx, 0.0f);
      glVertex3f(x1, y, z0);
      glEnd();
      y += dy;
    }
    break;
  case -2:
    /* the view direction is most nearly along negative y */
    y= y1;
    /* draw planes in back to front order */
    for(j= ny-1; j >=0; j--) {
      /* load the poper texture plane */
      for(i= 0; i < nx; i++) {
        for(k= 0; k < nz; k++) {
          indt= 4*(k+i*nz);
          dt= 4*data[i+j*nx+k*nx*ny];
          subImage[indt  ]= ctab[dt];
          subImage[indt+1]= ctab[dt+1];
          subImage[indt+2]= ctab[dt+2];
          subImage[indt+3]= ctab[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nz, nx, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x0, y, z0);
      glTexCoord2f(0.0f, fracz);
      glVertex3f(x0, y, z1);
      glTexCoord2f(fracx, fracz);
      glVertex3f(x1, y, z1);
      glTexCoord2f(fracx, 0.0f);
      glVertex3f(x1, y, z0);
      glEnd();
      y -= dy;
    }
    break;
  case 3:
    /* the view direction is most nearly along z */
    z= z0;
    /* draw planes in back to front order */
    for(k= 0; k < nz; k++) {
      /* load the poper texture plane */
      for(j= 0; j < ny; j++) {
        for(i= 0; i < nx; i++) {
          indt= 4*(i+j*nx);
          dt= 4*data[i+j*nx+k*nx*ny];
          subImage[indt  ]= ctab[dt];
          subImage[indt+1]= ctab[dt+1];
          subImage[indt+2]= ctab[dt+2];
          subImage[indt+3]= ctab[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nx, ny, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBindTexture(GL_TEXTURE_2D, texName);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x0, y0, z);
      glTexCoord2f(0.0f, fracy);
      glVertex3f(x0, y1, z);
      glTexCoord2f(fracx, fracy);
      glVertex3f(x1, y1, z);
      glTexCoord2f(fracx, 0.0f);
      glVertex3f(x1, y0, z);
      glEnd();
      z += dz;
    }
    break;
  case -3:
    /* the view direction is most nearly along negative z */
    z= z1;
    /* draw planes in back to front order */
    for(k= nz-1; k >= 0; k--) {
      /* load the poper texture plane */
      for(j= 0; j < ny; j++) {
        for(i= 0; i < nx; i++) {
          indt= 4*(i+j*nx);
          dt= 4*data[i+j*nx+k*nx*ny];
          subImage[indt  ]= ctab[dt];
          subImage[indt+1]= ctab[dt+1];
          subImage[indt+2]= ctab[dt+2];
          subImage[indt+3]= ctab[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nx, ny, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBindTexture(GL_TEXTURE_2D, texName);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x0, y0, z);
      glTexCoord2f(0.0f, fracy);
      glVertex3f(x0, y1, z);
      glTexCoord2f(fracx, fracy);
      glVertex3f(x1, y1, z);
      glTexCoord2f(fracx, 0.0f);
      glVertex3f(x1, y0, z);
      glEnd();
      z -= dz;
    }
    break;
  }
#else
  /* the view direction is most nearly along z */
  z= 0.0;
  y0= 0.0;
  x0= 0.0;
  y1= y0+(ny-1)*dy;
  x1= x0+(nx-1)*dx;
  /* draw planes in back to front order */
  for(k= 0; k < nz; k++) {
    /* load the poper texture plane */
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        subImage[4*(i+j*nnx)]= (int)((i+j+k)*255.0/((nx-1.0)+(ny-1.0))) %255;
        subImage[4*(i+j*nnx)+1]= (int)((i+j+k)*255.0/((nx-1.0)+(ny-1.0))+80) %255;
        subImage[4*(i+j*nnx)+2]= (int)((i+j+k)*255.0/((nx-1.0)+(ny-1.0))+160) %255;
        subImage[4*(i+j*nnx)+3]= 50;
      }
    }
    /* replace the portion of the texture for which there is data */
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nnx, nny, GL_RGBA, 
                    GL_UNSIGNED_BYTE, subImage);
    glBindTexture(GL_TEXTURE_2D, texName);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f);
    glVertex3f(x0, y0, z);
    glTexCoord2f(0.0f, fracy);
    glVertex3f(x0, y1, z);
    glTexCoord2f(fracx, fracy);
    glVertex3f(x1, y1, z);
    glTexCoord2f(fracx, 0.0f);
    glVertex3f(x1, y0, z);
    glEnd();
    z += dz;
  }
#endif

  glDeleteTextures(1, &texName);
  free(newTexImage);
  free(subImage);
  yglEndTex2d();
}

void yglTexcell2(long nx, long ny, long nz, double delta[3], 
                 unsigned char *tex)
{
  /* The input is a 3D array of RGBA's and a vector of
     zone sizes.
  */
  long i, j, k, n, nnx, nny, nnz, nh, nv, nhs, nvs, indt, dt, dir;
  float fracx, fracy, fracz;
  float x0, y0, z0, x1, y1, z1, dx, dy, dz;
  float x, y, z;
  double absx, absy, absz;
  unsigned char *newTexImage, *subImage;
#ifdef DO_TEX_TEST
  GLint format, redSize, greenSize, blueSize, alphaSize;
  GLint texWidth, texHeight;
#endif

  yglPrepTex2d();

  /* Textures must be a power of 2 in size. Find the smallest 
     power of 2 that will hold the input array */
  n= nx-1;
  nnx= 1;
  if(n) {
    do {
      nnx *= 2;
    } while(n= n/2);
  }
  fracx= (float)nx/(float)nnx;
  n= ny-1;
  nny= 1;
  if(n) {
    do {
      nny *= 2;
    } while(n= n/2);
  }
  fracy= (float)ny/(float)nny;
  n= nz-1;
  nnz= 1;
  if(n) {
    do {
      nnz *= 2;
    } while(n= n/2);
  }
  fracz= (float)nz/(float)nnz;

  /* There are nx by ny by nz cells to be volume rendered.
     FOR NOW, assume that the cells are aligned with the
     axes. */
  dx= (float) delta[0];
  dy= (float) delta[1];
  dz= (float) delta[2];

#ifdef USE_DATA
  /* find the axis most closely matching the viewing direction */
  absx= glCurrWin3d->view[0] < 0 ? -glCurrWin3d->view[0] : glCurrWin3d->view[0];
  absy= glCurrWin3d->view[1] < 0 ? -glCurrWin3d->view[1] : glCurrWin3d->view[1];
  absz= glCurrWin3d->view[2] < 0 ? -glCurrWin3d->view[2] : glCurrWin3d->view[2];
  if(absx < absy) {
    if(absy < absz) {
      /* the view direction is most nearly along z */
      dir= 3;
      if(glCurrWin3d->view[2] < 0) dir= -3;
      nh= nnx;
      nv= nny;
      nhs= nx;
      nvs= ny;
    } else {
      /* the view direction is most nearly along y */
      dir= 2;
      if(glCurrWin3d->view[1] < 0) dir= -2;
      nh= nnz;
      nv= nnx;
      nhs= nz;
      nvs= nx;
    }
  } else {
    if(absx < absz) {
      /* the view direction is most nearly along z */
      dir= 3;
      if(glCurrWin3d->view[2] < 0) dir= -3;
      nh= nnx;
      nv= nny;
      nhs= nx;
      nvs= ny;
    } else {
      /* the view direction is most nearly along x */
      dir= 1;
      if(glCurrWin3d->view[0] < 0) dir= -1;
      nh= nny;
      nv= nnz;
      nhs= ny;
      nvs= nz;
    }
  }
#else
  nh= nnx;
  nv= nny;
#endif
  glGenTextures(1, &texName);
  glBindTexture(GL_TEXTURE_2D, texName);
  newTexImage= (unsigned char *)malloc(sizeof(unsigned char)*4*nh*nv);
  subImage= (unsigned char *)malloc(sizeof(unsigned char)*4*nhs*nvs);
  {
    int err;
    const GLubyte *errstr;
    err= glGetError();
    errstr= gluErrorString(err);
    err++;
  }

#ifdef DO_TEX_TEST
  glTexImage2D(GL_PROXY_TEXTURE_2D, 0, GL_RGBA, nnx, nny, 0, GL_RGBA, 
               GL_UNSIGNED_BYTE, NULL);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_RED_SIZE, &redSize);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_GREEN_SIZE, &greenSize);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_BLUE_SIZE, &blueSize);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_ALPHA_SIZE, &alphaSize);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_WIDTH, &texWidth);
  glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, 
                           GL_TEXTURE_HEIGHT, &texHeight);
#endif

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  /*  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); */
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, nh, nv, 0, GL_RGBA, 
               GL_UNSIGNED_BYTE, newTexImage);

  glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  {
    int err;
    const GLubyte *errstr;
    err= glGetError();
    errstr= gluErrorString(err);
    err++;
  }

#ifdef USE_DATA
  x0= 0.0;
  x1= x0+(nx-1)*dx;
  y0= 0.0;
  y1= y0+(ny-1)*dy;
  z0= 0.0;
  z1= z0+(nz-1)*dz;

#ifdef DEBUG_PRINT
  {
    char msg[100];
    sprintf(msg,"direction is %d\n",dir);
    p_stderr(msg);
  }
#endif
  switch(dir) {
  case 1:
    /* the view direction is most nearly along x */
    x= x0;
    /* draw planes in back to front order */
    for(i= 0; i < nx; i++) {
      /* load the proper texture plane */
      for(k= 0; k < nz; k++) {
        for(j= 0; j < ny; j++) {
          indt= 4*(j+k*ny);
          dt= 4*(i+j*nx+k*nx*ny);
          subImage[indt  ]= tex[dt];
          subImage[indt+1]= tex[dt+1];
          subImage[indt+2]= tex[dt+2];
          subImage[indt+3]= tex[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ny, nz, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x, y0, z0);
      glTexCoord2f(fracy, 0.0f);
      glVertex3f(x, y1, z0);
      glTexCoord2f(fracy, fracz);
      glVertex3f(x, y1, z1);
      glTexCoord2f(0.0f, fracz);
      glVertex3f(x, y0, z1);
      glEnd();
      x += dx;
    }
    break;
  case -1:
    /* the view direction is most nearly along negative x */
    x= x1;
    /* draw planes in back to front order */
    for(i= nx-1; i >= 0; i--) {
      /* load the proper texture plane */
      for(k= 0; k < nz; k++) {
        for(j= 0; j < ny; j++) {
          indt= 4*(j+k*ny);
          dt= 4*(i+j*nx+k*nx*ny);
          subImage[indt  ]= tex[dt];
          subImage[indt+1]= tex[dt+1];
          subImage[indt+2]= tex[dt+2];
          subImage[indt+3]= tex[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ny, nz, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x, y0, z0);
      glTexCoord2f(0.0f, fracz);
      glVertex3f(x, y0, z1);
      glTexCoord2f(fracy, fracz);
      glVertex3f(x, y1, z1);
      glTexCoord2f(fracy, 0.0f);
      glVertex3f(x, y1, z0);
      glEnd();
      x -= dx;
    }
    break;
  case 2:
    /* the view direction is most nearly along y */
    y= y0;
    /* draw planes in back to front order */
    for(j= 0; j < ny; j++) {
      /* load the poper texture plane */
      for(k= 0; k < nz; k++) {
        for(i= 0; i < nx; i++) {
          indt= 4*(i+k*nx);
          dt= 4*(i+j*nx+k*nx*ny);
          subImage[indt  ]= tex[dt];
          subImage[indt+1]= tex[dt+1];
          subImage[indt+2]= tex[dt+2];
          subImage[indt+3]= tex[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nx, nz, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x0, y, z0);
      glTexCoord2f(0.0f, fracz);
      glVertex3f(x0, y, z1);
      glTexCoord2f(fracx, fracz);
      glVertex3f(x1, y, z1);
      glTexCoord2f(fracx, 0.0f);
      glVertex3f(x1, y, z0);
      glEnd();
      y += dy;
    }
    break;
  case -2:
    /* the view direction is most nearly along negative y */
    y= y1;
    /* draw planes in back to front order */
    for(j= ny-1; j >=0; j--) {
      /* load the poper texture plane */
      for(k= 0; k < nz; k++) {
        for(i= 0; i < nx; i++) {
          indt= 4*(i+k*nx);
          dt= 4*(i+j*nx+k*nx*ny);
          subImage[indt  ]= tex[dt];
          subImage[indt+1]= tex[dt+1];
          subImage[indt+2]= tex[dt+2];
          subImage[indt+3]= tex[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nx, nz, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x0, y, z0);
      glTexCoord2f(fracx, 0.0f);
      glVertex3f(x1, y, z0);
      glTexCoord2f(fracx, fracz);
      glVertex3f(x1, y, z1);
      glTexCoord2f(0.0f, fracz);
      glVertex3f(x0, y, z1);
      glEnd();
      y -= dy;
    }
    break;
  case 3:
    /* the view direction is most nearly along z */
    z= z0;
    /* draw planes in back to front order */
    for(k= 0; k < nz; k++) {
      /* load the poper texture plane */
      for(j= 0; j < ny; j++) {
        for(i= 0; i < nx; i++) {
          indt= 4*(i+j*nx);
          dt= 4*(i+j*nx+k*nx*ny);
          subImage[indt  ]= tex[dt];
          subImage[indt+1]= tex[dt+1];
          subImage[indt+2]= tex[dt+2];
          subImage[indt+3]= tex[dt+3];
        }
      }
      glBindTexture(GL_TEXTURE_2D, texName);
      {
        int err;
        const GLubyte *errstr;
        err= glGetError();
        errstr= gluErrorString(err);
        err++;
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nx, ny, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      {
        int err;
        const GLubyte *errstr;
        err= glGetError();
        errstr= gluErrorString(err);
        err++;
      }
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x0, y0, z);
      glTexCoord2f(fracx, 0.0f);
      glVertex3f(x1, y0, z);
      glTexCoord2f(fracx, fracy);
      glVertex3f(x1, y1, z);
      glTexCoord2f(0.0f, fracy);
      glVertex3f(x0, y1, z);
      glEnd();
      z += dz;
    }
    break;
  case -3:
    /* the view direction is most nearly along negative z */
    z= z1;
    /* draw planes in back to front order */
    for(k= nz-1; k >= 0; k--) {
      /* load the poper texture plane */
      for(j= 0; j < ny; j++) {
        for(i= 0; i < nx; i++) {
          indt= 4*(i+j*nx);
          dt= 4*(i+j*nx+k*nx*ny);
          subImage[indt  ]= tex[dt];
          subImage[indt+1]= tex[dt+1];
          subImage[indt+2]= tex[dt+2];
          subImage[indt+3]= tex[dt+3];
        }
      }
      /* replace the portion of the texture for which there is data */
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nx, ny, GL_RGBA, 
                      GL_UNSIGNED_BYTE, subImage);
      glBindTexture(GL_TEXTURE_2D, texName);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0f, 0.0f);
      glVertex3f(x0, y0, z);
      glTexCoord2f(0.0f, fracy);
      glVertex3f(x0, y1, z);
      glTexCoord2f(fracx, fracy);
      glVertex3f(x1, y1, z);
      glTexCoord2f(fracx, 0.0f);
      glVertex3f(x1, y0, z);
      glEnd();
      z -= dz;
    }
    break;
  }
#else
  /* the view direction is most nearly along z */
  z= 0.0;
  y0= 0.0;
  x0= 0.0;
  y1= y0+(ny-1)*dy;
  x1= x0+(nx-1)*dx;
  /* draw planes in back to front order */
  for(k= 0; k < nz; k++) {
    /* load the poper texture plane */
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        subImage[4*(i+j*nnx)]= (int)((i+j+k)*255.0/((nx-1.0)+(ny-1.0))) %255;
        subImage[4*(i+j*nnx)+1]= (int)((i+j+k)*255.0/((nx-1.0)+(ny-1.0))+80) %255;
        subImage[4*(i+j*nnx)+2]= (int)((i+j+k)*255.0/((nx-1.0)+(ny-1.0))+160) %255;
        subImage[4*(i+j*nnx)+3]= 50;
      }
    }
    /* replace the portion of the texture for which there is data */
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nnx, nny, GL_RGBA, 
                    GL_UNSIGNED_BYTE, subImage);
    glBindTexture(GL_TEXTURE_2D, texName);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f);
    glVertex3f(x0, y0, z);
    glTexCoord2f(0.0f, fracy);
    glVertex3f(x0, y1, z);
    glTexCoord2f(fracx, fracy);
    glVertex3f(x1, y1, z);
    glTexCoord2f(fracx, 0.0f);
    glVertex3f(x1, y0, z);
    glEnd();
    z += dz;
  }
#endif

  glDeleteTextures(1, &texName);
  free(newTexImage);
  free(subImage);
  yglEndTex2d();
}

void yglPrepCubeTex(void)
{
  /*  glDisable(GL_LIGHT0); */
  glEnable(GL_TEXTURE_CUBE_MAP_EXT);
  /*  glEnable(GL_BLEND); */
  /*  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); */
  /*  glDisable(GL_LIGHT0); */
  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
}

void yglEndCubeTex(void)
{
  glDisable(GL_TEXTURE_GEN_S);
  glDisable(GL_TEXTURE_GEN_T);
  glDisable(GL_TEXTURE_GEN_R);
  /* Other texture modes (1D, 2D, or 3D) are overridden by cube maps */
  glDisable(GL_TEXTURE_CUBE_MAP_EXT);
  /*  glEnable(GL_LIGHT0); */
  /*  glDisable(GL_BLEND); */
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
}

static void yglGenCubeTex(void)
{
  long i, j, icen, iwid;
  unsigned char alpha;
  double wid2;
 
  /* These textures will be used for lighting using cube maps.
     The light is white, brightest at the center of the face,
     and has a width of iwid texels.
     It is only on the front and back faces.
     Lighting is handled by making all texels white and varying
     the opacity (high in the specular highlight, low
     elsewhere) */
  memset(cubeMaps, 127, 6*ygl_cubemap_size*ygl_cubemap_size*4);
  icen= ygl_cubemap_size/2;
  iwid= (long) (ygl_cubemap_size/1.5);
  wid2= iwid*iwid;
  for(i= 0; i < ygl_cubemap_size; i++) {
    for(j= 0; j < ygl_cubemap_size; j++) {
      alpha= (unsigned char)(255.0*exp( -((i-icen)*(i-icen)+(j-icen)*(j-icen))/wid2 ) );
      cubeMaps[4][i][j][3]= alpha;
      cubeMaps[5][i][j][3]= alpha;
    }
  }
}

void yglLdCubeTex(void)
{
  int i;

  yglMakeCurrent(glCurrWin3d);
  if( !yglQueryTexCube() ) return;

  /* generate the cube map textures, if necessary */
  if(texNameCube) {
    glBindTexture(GL_TEXTURE_CUBE_MAP_EXT, texNameCube);
  } else {
    yglGenCubeTex();  /* generate the texture data */
    glGenTextures(1, &texNameCube);  /* make the texture object */
    glBindTexture(GL_TEXTURE_CUBE_MAP_EXT, texNameCube);
    glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    for(i= 0; i < 6; i++) {
      glTexImage2D(faceTarget[i], 0, GL_RGBA8,
                   ygl_cubemap_size, ygl_cubemap_size, 0,
                   GL_RGBA, GL_UNSIGNED_BYTE, &(cubeMaps[i][0][0][0]) );
    }
  }
  glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP_EXT);
  glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP_EXT);
  glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP_EXT);
  glEnable(GL_TEXTURE_CUBE_MAP_EXT);
  glEnable(GL_TEXTURE_GEN_S);
  glEnable(GL_TEXTURE_GEN_T);
  glEnable(GL_TEXTURE_GEN_R);
  glEnable(GL_NORMALIZE);
  /*
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  */
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
}

void yglPrepTex2d(void)
{
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);	/* No hidden surface removal for translucent */
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_TEXTURE_2D);
}

void yglEndTex2d(void)
{
  /* 1D textures are not available while 2D textures are enabled */
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);	/* Hidden surface removal */
}

void yglPrepTex3d(void)
{
  glDisable(GL_LIGHTING);
  /*  glDisable(GL_DEPTH_TEST); */	/* No hidden surface removal for translucent */
  glDepthMask(GL_FALSE);	/* do not update depth buffer, but leave enabled */
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(glCurrWin3d->myGL_TEXTURE_3D);
}

void yglEndTex3d(void)
{
  /* 2D textures are not available while 3D textures are enabled */
  glDisable(glCurrWin3d->myGL_TEXTURE_3D);
  glDisable(GL_BLEND);
  glEnable(GL_LIGHTING);
  glDepthMask(GL_TRUE);	/* resume updates of update depth buffer */
  /*  glEnable(GL_DEPTH_TEST); */	/* Hidden surface removal */
}

/* WARNING: Under Visual C++ 6.0 this function seg faults if
   optimized.
*/
#pragma optimize( "", off )
void yglLdTex3d(long nx, long ny, long nz, unsigned char *tex)
{
  /* The input is a 3D array of RGBA's. */
  long i, j, k, n, ind, indt, indup, n2x, n2y, n2z, nxh, nyh;
  long tex_siz;

  if(!yglQueryTex3d(glCurrWin3d)) {
    YError("This computer does not have 3D textures");
  }
  /* Textures must be a power of 2 in size. Find the smallest 
     power of 2 that will hold the input array. */
  n= nx-1;
  n2x= 1;
  if(n) {
    do {
      n2x *= 2;
    } while(n= n/2);
  }
  fracx3d= (float)nx/(float)n2x;
  n= ny-1;
  n2y= 1;
  if(n) {
    do {
      n2y *= 2;
    } while(n= n/2);
  }
  fracy3d= (float)ny/(float)n2y;
  n= nz-1;
  n2z= 1;
  if(n) {
    do {
      n2z *= 2;
    } while(n= n/2);
  }
  fracz3d= (float)nz/(float)n2z;

  if(texImage3d) {
    /* use the previous array if it was the same size */
    if(nx3d != n2x || ny3d != n2y || nz3d != n2z) {
      free(texImage3d);
      nx3d= n2x;
      ny3d= n2y;
      nz3d= n2z;
      tex_siz= sizeof(unsigned char)*4*nx3d*ny3d*nz3d;
      texImage3d= (unsigned char *)malloc(tex_siz);
    }
  } else {
    /* no existing texture, so make a new one */
    nx3d= n2x;
    ny3d= n2y;
    nz3d= n2z;
    tex_siz= sizeof(unsigned char)*4*nx3d*ny3d*nz3d;
    texImage3d= (unsigned char *)malloc(tex_siz);
  }
  for(k= 0; k < nz; k++) {
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        ind= 4*(i+j*nx+k*nx*ny);
        indt= 4*(i+j*nx3d+k*nx3d*ny3d);
        texImage3d[indt]= tex[ind];
        texImage3d[indt+1]= tex[ind+1];
        texImage3d[indt+2]= tex[ind+2];
        texImage3d[indt+3]= tex[ind+3];
      }
    }
  }
  /* If linear interpolation of texture values is turned on, it will
     interpolate off the edge of the valid texture data unless the input
	 array is a power of two in size. Deal with this by duplicating
	 boundary values. */
  if(nx < nx3d) {
    i= nx-1;
    for(k= 0; k < nz; k++) {
      for(j= 0; j < ny; j++) {
        indt= 4*(i+j*nx3d+k*nx3d*ny3d);
        indup= indt+4;
        texImage3d[indup]= texImage3d[indt];
        texImage3d[indup+1]= texImage3d[indt+1];
        texImage3d[indup+2]= texImage3d[indt+2];
        texImage3d[indup+3]= texImage3d[indt+3];
      }
    }
    nxh= nx+1;
  } else {
    nxh= nx;
  }
  if(ny < ny3d) {
    j= ny-1;
    for(k= 0; k < nz; k++) {
      for(i= 0; i < nxh; i++) {
        indt= 4*(i+j*nx3d+k*nx3d*ny3d);
        indup= indt+4*nx3d;
        texImage3d[indup]= texImage3d[indt];
        texImage3d[indup+1]= texImage3d[indt+1];
        texImage3d[indup+2]= texImage3d[indt+2];
        texImage3d[indup+3]= texImage3d[indt+3];
      }
    }
    nyh= ny+1;
  } else {
    nyh= ny;
  }
  if(nz < nz3d) {
    k= nz-1;
    for(j= 0; j < nyh; j++) {
      for(i= 0; i < nxh; i++) {
        indt= 4*(i+j*nx3d+k*nx3d*ny3d);
        indup= indt+4*nx3d*ny3d;
        texImage3d[indup]= texImage3d[indt];
        texImage3d[indup+1]= texImage3d[indt+1];
        texImage3d[indup+2]= texImage3d[indt+2];
        texImage3d[indup+3]= texImage3d[indt+3];
      }
    }
  }
  yglMakeCurrent(glCurrWin3d);

  if(!texName3d) {
    glGenTextures(1, &texName3d);
  }
  glBindTexture(glCurrWin3d->myGL_TEXTURE_3D, texName3d);
  glTexParameteri(glCurrWin3d->myGL_TEXTURE_3D, GL_TEXTURE_WRAP_R_EXT, GL_CLAMP);
  glTexParameteri(glCurrWin3d->myGL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(glCurrWin3d->myGL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(glCurrWin3d->myGL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  /*  glTexParameteri(glCurrWin3d->myGL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); */
  glTexParameteri(glCurrWin3d->myGL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  myglTexImage3D(glCurrWin3d->myGL_TEXTURE_3D, 0, GL_RGBA, nx3d, ny3d, nz3d, 0, GL_RGBA, 
                    GL_UNSIGNED_BYTE, texImage3d);
#undef DO_TEX3_TEST
#ifdef DO_TEX3_TEST
  {
    GLint redSize, greenSize, blueSize, alphaSize;
    GLint texWidth, texHeight, texDepth;

    myglTexImage3D(myGL_PROXY_TEXTURE_3D, 0, GL_RGBA, nx3d, ny3d, nz3d,
                   0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glGetTexLevelParameteriv(myGL_PROXY_TEXTURE_3D, 0, 
                             GL_TEXTURE_RED_SIZE, &redSize);
    glGetTexLevelParameteriv(myGL_PROXY_TEXTURE_3D, 0, 
                             GL_TEXTURE_GREEN_SIZE, &greenSize);
    glGetTexLevelParameteriv(myGL_PROXY_TEXTURE_3D, 0, 
                             GL_TEXTURE_BLUE_SIZE, &blueSize);
    glGetTexLevelParameteriv(myGL_PROXY_TEXTURE_3D, 0, 
                             GL_TEXTURE_ALPHA_SIZE, &alphaSize);
    glGetTexLevelParameteriv(myGL_PROXY_TEXTURE_3D, 0, 
                             GL_TEXTURE_WIDTH, &texWidth);
    glGetTexLevelParameteriv(myGL_PROXY_TEXTURE_3D, 0, 
                             GL_TEXTURE_HEIGHT, &texHeight);
    glGetTexLevelParameteriv(myGL_PROXY_TEXTURE_3D, 0, 
                             GL_TEXTURE_DEPTH_EXT, &texDepth);
    if(texDepth <= 0) {  /* to provide a space for a breakpoint inside scope */
      texDepth= 1;
    }
  }
#endif

  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
}
#pragma optimize( "", on )

void yglTexPoly(long nvert, float *verts, float *texverts)
{
  long i;

  /* Draw a 3Dtextured polygon. */
  glBindTexture(glCurrWin3d->myGL_TEXTURE_3D, texName3d);
  glBegin(GL_POLYGON);
  for(i= 0; i < nvert; i++) {
    glTexCoord3f(texverts[3*i], texverts[3*i+1], texverts[3*i+2]);
    glVertex3f(verts[3*i], verts[3*i+1], verts[3*i+2]);
  }
  glEnd();
}

void yglTexTris(long ntri, float *verts, float *texverts)
{
  long i, nd;

  /* Draw an array of 3Dtextured triangles. */
  glBindTexture(glCurrWin3d->myGL_TEXTURE_3D, texName3d);
  glBegin(GL_TRIANGLES);
  for(i= 0; i < ntri; i++) {
    nd= 9*i;  /* 3 coords per vertex, 3 vertices per tri */
    glTexCoord3fv(texverts+nd);
    glVertex3fv(verts+nd);
    glTexCoord3fv(texverts+nd+3);
    glVertex3fv(verts+nd+3);
    glTexCoord3fv(texverts+nd+6);
    glVertex3fv(verts+nd+6);
  }
  glEnd();
}

void slice_box(double s, double sbox[8], double *origin, 
               double *len, TriArrayGrp *tris)
{
  long dims[]= {4, 4, 4};
  long offsets[]= {0, 0, 0, 2, 2, 2};
  yPoint3D xyztmp[64], grdtmp[64];
  double vartmp[64], v2tmp[64], norigin[3];
  unsigned char flagtmp[64];

  /* This function finds the intersection of a plane and a 
     "rectangular" hexahedron.
     The caller must make tris large enough to hold 4
     triangles (the worst case).
     origin is the coordinate of the lowest numbered point
	 in the full array. The first point in the chunk is a guard
	 point and is thus outside the full array. 
  */
  norigin[0]= origin[0]-len[0];
  norigin[1]= origin[1]-len[1];
  norigin[2]= origin[2]-len[2];
  /* Contour on a 2-by-2-by-2 grid with no guard cells.
     All three sides have length 1.0.
  */
  ycInitCartPcen(dims, offsets, len, norigin, sbox, 0);
  /* iso-surface at the distance, s, of the current plane */
  ycContourTet_array(0, dims, s, 
            vartmp, v2tmp, xyztmp, grdtmp, flagtmp, tris);
}

void yglTex3dbox(double ds, double *origin, double *len)
{
  double s, smin, smax, facx, facy, facz, difs;
  double s0, s1, s2, s3, s4, s5, s6, s7;
  double sbox[8];
  yPoint3D *dverts;
  float *verts, *texverts;
  long i, j, ntri, maxtri, nvert, nslab;

  /* NOTE: len is the length of the sides of the whole volume,
     not the size of an individual cell.
  */
  maxtri= 4;
  facx= fracx3d/len[0];
  facy= fracy3d/len[1];
  facz= fracz3d/len[2];

  if(!tris_tex3d) {
    tris_tex3d= (TriArrayGrp *) malloc(sizeof(TriArrayGrp));
    tris_tex3d->next= 0;
    tris_tex3d->cellIDs= (long *) malloc(maxtri*sizeof(long));
    tris_tex3d->xyzverts= (yPoint3D *) malloc(3*maxtri*sizeof(yPoint3D));
    tris_tex3d->normals= (yPoint3D *) malloc(3*maxtri*sizeof(yPoint3D));
  }
  verts= (float *) malloc(sizeof(float)*3*3*maxtri);
  texverts= (float *) malloc(sizeof(float)*3*3*maxtri);
  dverts= tris_tex3d->xyzverts;

  s0= glCurrWin3d->view[0]*origin[0]+glCurrWin3d->view[1]*origin[1]+glCurrWin3d->view[2]*origin[2];
  sbox[0]= s0;
  smin= smax= s0;
  sbox[1]= s1= s0+glCurrWin3d->view[0]*len[0];
  if(s1 < smin) smin= s1;
  if(s1 > smax) smax= s1;
  sbox[2]= s2= s0+glCurrWin3d->view[1]*len[1];
  if(s2 < smin) smin= s2;
  if(s2 > smax) smax= s2;
  sbox[3]= s3= s0+glCurrWin3d->view[0]*len[0]+glCurrWin3d->view[1]*len[1];
  if(s3 < smin) smin= s3;
  if(s3 > smax) smax= s3;
  sbox[4]= s4= s0+glCurrWin3d->view[2]*len[2];
  if(s4 < smin) smin= s4;
  if(s4 > smax) smax= s4;
  sbox[5]= s5= s4+glCurrWin3d->view[0]*len[0];
  if(s5 < smin) smin= s5;
  if(s5 > smax) smax= s5;
  sbox[6]= s6= s4+glCurrWin3d->view[1]*len[1];
  if(s6 < smin) smin= s6;
  if(s6 > smax) smax= s6;
  sbox[7]= s7= s4+glCurrWin3d->view[0]*len[0]+glCurrWin3d->view[1]*len[1];
  if(s7 < smin) smin= s7;
  if(s7 > smax) smax= s7;
  nslab= (long) ((smax-smin)/ds);
  difs= (smax-smin)/nslab;
  s= smin+0.5*difs;

  yglPrepTex3d();

  for(i= 0; i < nslab; i++) {
    slice_box(s, sbox, origin, len, tris_tex3d);
    ntri= tris_tex3d->numTri;
    if(ntri > 0) {
      nvert= 3*ntri;
      for(j= 0; j < nvert; j++) {
        /* The texture coordinate is the fraction of the distance
           between the minimum and maximum world coord value.
           This is adjusted for any extra space added when making
           the texture a power of 2 in size.
        */
        verts[3*j]= (float) dverts[j].x;
        texverts[3*j]= (float) ((verts[3*j]-origin[0])*facx);
        verts[3*j+1]= (float) dverts[j].y;
        texverts[3*j+1]= (float) ((verts[3*j+1]-origin[1])*facy);
        verts[3*j+2]= (float) dverts[j].z;
        texverts[3*j+2]= (float) ((verts[3*j+2]-origin[2])*facz);
      }
      yglTexTris(ntri, verts, texverts);
    }
    s += difs;
  }
  yglEndTex3d();
}
