/*
 * $Id: glustub.c,v 1.1 2005-09-18 22:07:53 dhmunro Exp $
 * GLU is not always present alongside GL (e.g.- Windows)
 * here are the three GLU routines yorgl needs
 * -- these are loosely based on Mesa-4.0.1 sources
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "glcode.h"
#include <math.h>

static void gs_normalize(GLdouble x[]);

void
gluLookAt(GLdouble eyex, GLdouble eyey, GLdouble eyez,
	  GLdouble centerx, GLdouble centery, GLdouble centerz,
	  GLdouble upx, GLdouble upy, GLdouble upz)
{
  GLdouble m[16];
  GLdouble *x = &m[0];
  GLdouble *y = &m[1];
  GLdouble *z = &m[2];

  z[0] = eyex - centerx;
  z[4] = eyey - centery;
  z[8] = eyez - centerz;
  gs_normalize(z);

  x[0] = upy*z[8] - upz*z[4];
  x[4] = upz*z[0] - upx*z[8];
  x[8] = upx*z[4] - upy*z[0];
  gs_normalize(x);

  y[0] = z[4]*x[8] - z[8]*x[4];
  y[4] = z[8]*x[0] - z[0]*x[8];
  y[8] = z[0]*x[4] - z[4]*x[0];
  gs_normalize(y);

  m[3] = m[7] = m[11] =
    m[12] = m[13] = m[14] = 0.0;
  m[15] = 1.0;
  glMultMatrixd(m);

  glTranslated(-eyex, -eyey, -eyez);
}

static void
gs_normalize(GLdouble x[])
{
  double len = x[0]*x[0] + x[4]*x[4] + x[8]*x[8];
  if (len) {
    len = 1./sqrt(len);
    x[0] *= len;
    x[4] *= len;
    x[8] *= len;
  }
}

void
gluPerspective(GLdouble fovy, GLdouble aspect, GLdouble near_, GLdouble far_)
{
   GLdouble m[16], rdz;

   fovy = tan(fovy*3.14159265358979324/360.0);
   rdz = 1. / (far_ - near_);

   m[0] = 1. / (fovy*aspect);
   m[1] = m[2] = m[3] = m[4] = 0.0;
   m[5] = 1. / fovy;
   m[6] = m[7] = m[8] = m[9] = 0.0;
   m[10] = -(far_ + near_) * rdz;
   m[11] = -1.0;
   m[12] = m[13] = 0.0;
   m[14] = -2.0 * far_ * near_ * rdz;
   m[15] = 0.0;

   glMultMatrixd(m);
}

const GLubyte *
gluErrorString(GLenum code)
{
  if (code == GL_NO_ERROR)
    return (GLubyte *) "no error";
  else if (code == GL_INVALID_ENUM)
    return (GLubyte *) "GLenum argument out of range";
  else if (code == GL_INVALID_VALUE)
    return (GLubyte *) "numeric argument out of range";
  else if (code == GL_INVALID_OPERATION)
    return (GLubyte *) "operation illegal in current state";
  else if (code == GL_STACK_OVERFLOW)
    return (GLubyte *) "command would cause a stack overflow";
  else if (code == GL_STACK_UNDERFLOW)
    return (GLubyte *) "command would cause a stack underflow";
  else if (code == GL_OUT_OF_MEMORY)
    return (GLubyte *) "not enough memory to execute command";
  else
    return 0;
}
