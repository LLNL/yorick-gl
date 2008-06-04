/*
 * $Id: glTarray.c,v 1.3 2008-06-04 05:56:04 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
#include "glcode.h"
#include "glfunc.h"
#include "glStrips.h"
#include "glWrappers.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*
#define CHEK_ERROR(x)	yygl_chek_error(x)
*/
#define CHEK_ERROR(x)


void yglTarray(long smooth, long ntri, float *xyz, float *norm, float *colr, 
               long edge, long cpervrt, long emit)
{
  long i, base;
  float oldRGBA[4]= {-1.0, -1.0, -1.0, 1.0};
  float ambi[]= {1.0f, 1.0f, 1.0f, 1.0f};

  /* draw an array of triangles */
  if(ntri <= 0) return;
  if(alpha_pass) return;
  if(emit) {
    /*  this creates diffuse light not connected to any light source */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambi);
    /* turn off all light associated with the light source */
    glDisable(GL_LIGHT0);
  } else {
    if(smooth) {
      /* use smooth shading */
      yglSetShade(1);
    } else {
      /* use flat shading */
      yglSetShade(0);
    }
  }
  yglUpdateProperties();

  base= 0;
  glBegin(GL_TRIANGLES);
  if(emit) {
    /* NOTE: no normals are supplied when using emissive colors */
    if(cpervrt) {
      for(i= 0; i < ntri; i++) {
        glColor3fv(colr);
        colr += 3;
        glVertex3fv(xyz+base);
        glColor3fv(colr);
        colr += 3;
        glVertex3fv(xyz+base+3);
        glColor3fv(colr);
        colr += 3;
        glVertex3fv(xyz+base+6);
        base += 9;
      }
    } else {
      for(i= 0; i < ntri; i++) {
        if(colr[0] != oldRGBA[0] || colr[1] != oldRGBA[1] || 
           colr[2] != oldRGBA[2]) {
          oldRGBA[0]= colr[0];
          oldRGBA[1]= colr[1];
          oldRGBA[2]= colr[2];
          glColor3fv(oldRGBA);
        }
        colr += 3;
        /* NOTE: It makes no difference in local performance if the coords
           are converted to floats before being stored in the display list */
        glVertex3fv(xyz+base);
        glVertex3fv(xyz+base+3);
        glVertex3fv(xyz+base+6);
        base += 9;
      }
    }
  } else {
    if(cpervrt) {
      for(i= 0; i < ntri; i++) {
        if(smooth) {
          glColor3fv(colr);
          colr += 3;
          glNormal3fv(norm+base);
          glVertex3fv(xyz+base);
          glColor3fv(colr);
          colr += 3;
          glNormal3fv(norm+base+3);
          glVertex3fv(xyz+base+3);
          glColor3fv(colr);
          colr += 3;
          glNormal3fv(norm+base+6);
          glVertex3fv(xyz+base+6);
        } else {
          glColor3fv(colr);
          colr += 3;
          glNormal3fv(norm+3*i);
          glVertex3fv(xyz+base);
          glColor3fv(colr);
          colr += 3;
          glVertex3fv(xyz+base+3);
          glColor3fv(colr);
          colr += 3;
          glVertex3fv(xyz+base+6);
        }
        base += 9;
      }
    } else {
      for(i= 0; i < ntri; i++) {
        if(colr[0] != oldRGBA[0] || colr[1] != oldRGBA[1] || 
           colr[2] != oldRGBA[2]) {
          oldRGBA[0]= colr[0];
          oldRGBA[1]= colr[1];
          oldRGBA[2]= colr[2];
          glColor3fv(oldRGBA);
        }
        colr += 3;
        /* NOTE: It makes no difference in local performance if the coords
           are converted to floats before being stored in the display list */
        if(smooth) {
          glNormal3fv(norm+base);
          glVertex3fv(xyz+base);
          glNormal3fv(norm+base+3);
          glVertex3fv(xyz+base+3);
          glNormal3fv(norm+base+6);
          glVertex3fv(xyz+base+6);
        } else {
          glNormal3fv(norm+3*i);
          glVertex3fv(xyz+base);
          glVertex3fv(xyz+base+3);
          glVertex3fv(xyz+base+6);
        }
        base += 9;
      }
    }
  }
  glEnd();
  if(emit) {
    /* restore prior diffuse light not connected to any light source */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, glCurrWin3d->curr_ambient);
    glEnable(GL_LIGHT0);
  }
  CHEK_ERROR("yglTarray");
}

void yglTarrayAlpha(long smooth, long ntri, float *xyz, float *norm, 
                    float *colr, long edge, long cpervrt, long emit)
{
  long i, base;
  float oldRGBA[4]= {-1.0, -1.0, -1.0, 1.0};
  float ambi[]= {1.0f, 1.0f, 1.0f, 1.0f};

  /* draw an array of triangles */
  if(ntri <= 0) return;
{
  char msg[120];
  sprintf(msg, "in yglTarrayAlpha, alpha_pass is %d\n", alpha_pass);
  puts(msg);
}
  if(!alpha_pass) return;
{
  puts("drawing alpha tarray");
}
  if(emit) {
    /*  this creates diffuse light not connected to any light source */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambi);
    /* turn off all light associated with the light source */
    glDisable(GL_LIGHT0);
  } else {
    if(smooth) {
      /* use smooth shading */
      yglSetShade(1);
    } else {
      /* use flat shading */
      yglSetShade(0);
    }
  }
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);
  yglUpdateProperties();

  base= 0;
  glBegin(GL_TRIANGLES);
  if(emit) {
    /* NOTE: no normals are supplied when using emissive colors */
    if(cpervrt) {
      for(i= 0; i < ntri; i++) {
        glColor4fv(colr);
        colr += 4;
        glVertex3fv(xyz+base);
        glColor4fv(colr);
        colr += 4;
        glVertex3fv(xyz+base+3);
        glColor4fv(colr);
        colr += 4;
        glVertex3fv(xyz+base+6);
        base += 9;
      }
    } else {
      for(i= 0; i < ntri; i++) {
        if(colr[0] != oldRGBA[0] || colr[1] != oldRGBA[1] || 
           colr[2] != oldRGBA[2] || colr[3] != oldRGBA[3]) {
          oldRGBA[0]= colr[0];
          oldRGBA[1]= colr[1];
          oldRGBA[2]= colr[2];
          oldRGBA[3]= colr[3];
          glColor4fv(oldRGBA);
        }
        colr += 4;
        /* NOTE: It makes no difference in local performance if the coords
           are converted to floats before being stored in the display list */
        glVertex3fv(xyz+base);
        glVertex3fv(xyz+base+3);
        glVertex3fv(xyz+base+6);
        base += 9;
      }
    }
  } else {
    if(cpervrt) {
      for(i= 0; i < ntri; i++) {
        if(smooth) {
          glColor4fv(colr);
          colr += 4;
          glNormal3fv(norm+base);
          glVertex3fv(xyz+base);
          glColor4fv(colr);
          colr += 4;
          glNormal3fv(norm+base+3);
          glVertex3fv(xyz+base+3);
          glColor4fv(colr);
          colr += 4;
          glNormal3fv(norm+base+6);
          glVertex3fv(xyz+base+6);
        } else {
          glColor4fv(colr);
          colr += 4;
          glNormal3fv(norm+3*i);
          glVertex3fv(xyz+base);
          glColor4fv(colr);
          colr += 4;
          glVertex3fv(xyz+base+3);
          glColor4fv(colr);
          colr += 4;
          glVertex3fv(xyz+base+6);
        }
        base += 9;
      }
    } else {
      for(i= 0; i < ntri; i++) {
        if(colr[0] != oldRGBA[0] || colr[1] != oldRGBA[1] || 
           colr[2] != oldRGBA[2] || colr[3] != oldRGBA[3]) {
          oldRGBA[0]= colr[0];
          oldRGBA[1]= colr[1];
          oldRGBA[2]= colr[2];
          oldRGBA[3]= colr[3];
          glColor4fv(oldRGBA);
        }
        colr += 4;
        /* NOTE: It makes no difference in local performance if the coords
           are converted to floats before being stored in the display list */
        if(smooth) {
          glNormal3fv(norm+base);
          glVertex3fv(xyz+base);
          glNormal3fv(norm+base+3);
          glVertex3fv(xyz+base+3);
          glNormal3fv(norm+base+6);
          glVertex3fv(xyz+base+6);
        } else {
          glNormal3fv(norm+3*i);
          glVertex3fv(xyz+base);
          glVertex3fv(xyz+base+3);
          glVertex3fv(xyz+base+6);
        }
        base += 9;
      }
    }
  }
  glEnd();
  if(emit) {
    /* restore prior diffuse light not connected to any light source */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, glCurrWin3d->curr_ambient);
    glEnable(GL_LIGHT0);
  }
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  CHEK_ERROR("yglTarrayNoArrAlpha");
}

void yglTarrayCubeMapAlpha(long ntri, float *xyz, float *norm, float *colr, 
                           long cpervrt)
{
  long i, base;
  float old_red= -1.0;
  float old_green= -1.0;
  float old_blue= -1.0;
  float old_alpha= -1.0;

  /* Draw an array of triangles using a cube map texture 
     for lighting. Assumes the cube map texture has already
     been loaded and that cube map textures have been enabled. */
  if(ntri <= 0) return;
  if(!alpha_pass) return;

  if( !yglQueryTex3d(glCurrWin3d) ) return;
  if( !glCurrWin3d->hascubetex) return;
  base= 0;
  glBegin(GL_TRIANGLES);
  if(cpervrt) {
    /* one color per vertex */
    for(i= 0; i < ntri; i++) {
      glColor4fv(colr);
      colr += 4;
      glNormal3fv(norm+base);
      glVertex3fv(xyz+base);
      glColor4fv(colr);
      colr += 4;
      glNormal3fv(norm+base+3);
      glVertex3fv(xyz+base+3);
      glColor4fv(colr);
      colr += 4;
      glNormal3fv(norm+base+6);
      glVertex3fv(xyz+base+6);
      base += 9;
    }
  } else {
    /* one color per triangle */
    for(i= 0; i < ntri; i++) {
      if(colr[0] != old_red || colr[1] != old_green || 
         colr[2] != old_blue || colr[3] != old_alpha) {
        glColor4fv(colr);
        old_red= colr[0];
        old_green= colr[1];
        old_blue= colr[2];
        old_alpha= colr[3];
      }
      colr += 4;
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      glNormal3fv(norm+base);
      glVertex3fv(xyz+base);
      glNormal3fv(norm+base+3);
      glVertex3fv(xyz+base+3);
      glNormal3fv(norm+base+6);
      glVertex3fv(xyz+base+6);
      base += 9;
    }
  }
  glEnd();
  CHEK_ERROR("yglTarrayCubeMapAlpha");
}

void yglTarrayCubeMap(long ntri, float *xyz, float *norm, float *colr,
                      long cpervrt)
{
  long i, base;
  float old_colr[]= {-1.0, -1.0, -1.0, 1.0};

  /* Draw an array of triangles using a cube map texture 
     for lighting. Assumes the cube map texture has already
     been loaded and that cube map textures have been enabled. */
  if(ntri <= 0) return;
  if(alpha_pass) return;

  if( !yglQueryTexCube() ) return;
  yglLdCubeTex();
  yglPrepCubeTex();
  base= 0;
  glBegin(GL_TRIANGLES);
  if(cpervrt) {
    for(i= 0; i < ntri; i++) {
      glColor3fv(colr);
      colr += 3;
      glNormal3fv(norm+base);
      glVertex3fv(xyz+base);
      glColor3fv(colr);
      colr += 3;
      glNormal3fv(norm+base+3);
      glVertex3fv(xyz+base+3);
      glColor3fv(colr);
      colr += 3;
      glNormal3fv(norm+base+6);
      glVertex3fv(xyz+base+6);
      base += 9;
    }
  } else {
    for(i= 0; i < ntri; i++) {
      if(colr[0] != old_colr[0] || colr[1] != old_colr[1] || 
         colr[2] != old_colr[2]) {
        glColor4fv(old_colr);
        old_colr[0]= colr[0];
        old_colr[1]= colr[1];
        old_colr[2]= colr[2];
      }
      colr += 3;
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      glNormal3fv(norm+base);
      glVertex3fv(xyz+base);
      glNormal3fv(norm+base+3);
      glVertex3fv(xyz+base+3);
      glNormal3fv(norm+base+6);
      glVertex3fv(xyz+base+6);
      base += 9;
    }
  }
  glEnd();
  yglEndCubeTex();
  CHEK_ERROR("yglTarrayCubeMapAlpha");
}

void yglTarrayEmit(long do_alpha, long ntri, float *xyz, float *colr, 
                   long cpervrt)
{
  long i, base, colrsiz;
  float oldRGBA[4]= {-1.0, -1.0, -1.0, 1.0};
  float ambi[]= {1.0f, 1.0f, 1.0f, 1.0f};
  float def_ambi[]= {0.2f, 0.2f, 0.2f, 1.0f};

  /* draw an array of triangles */
  if(ntri <= 0) return;

  /*  this creates diffuse light not connected to any light */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambi);
  /* turn off all light associated with the light source */
  glDisable(GL_LIGHT0);
  if(do_alpha) {
    if(!alpha_pass) return;
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    colrsiz= 4;
  } else {
    if(alpha_pass) return;
    glDisable(GL_BLEND);
    colrsiz= 3;
  }

  base= 0;
  glBegin(GL_TRIANGLES);
  if(do_alpha) {
    if(cpervrt) {
      for(i= 0; i < ntri; i++) {
        glColor4fv(colr);
        colr += colrsiz;
        glVertex3fv(xyz+base);
        glColor4fv(colr);
        colr += colrsiz;
        glVertex3fv(xyz+base+3);
        glColor4fv(colr);
        colr += colrsiz;
        glVertex3fv(xyz+base+6);
        base += 9;
      }
    } else {
      for(i= 0; i < ntri; i++) {
        if(colr[0] != oldRGBA[0] || colr[1] != oldRGBA[1] || 
           colr[2] != oldRGBA[2]  || colr[3] != oldRGBA[3]) {
          oldRGBA[0]= colr[0];
          oldRGBA[1]= colr[1];
          oldRGBA[2]= colr[2];
          oldRGBA[3]= colr[3];
          glColor4fv(oldRGBA);
        }
        colr += colrsiz;
        /* NOTE: It makes no difference in local performance if the coords
           are converted to floats before being stored in the display list */
        glVertex3fv(xyz+base);
        glVertex3fv(xyz+base+3);
        glVertex3fv(xyz+base+6);
        base += 9;
      }
    }
  } else {
    if(cpervrt) {
      /* no alpha component */
      for(i= 0; i < ntri; i++) {
        glColor3fv(colr);
        colr += colrsiz;
        glVertex3fv(xyz+base);
        glColor3fv(colr);
        colr += colrsiz;
        glVertex3fv(xyz+base+3);
        glColor3fv(colr);
        colr += colrsiz;
        glVertex3fv(xyz+base+6);
        base += 9;
      }
    } else {
      /* no alpha component */
      for(i= 0; i < ntri; i++) {
        if(colr[0] != oldRGBA[0] || colr[1] != oldRGBA[1] || 
           colr[2] != oldRGBA[2]) {
          oldRGBA[0]= colr[0];
          oldRGBA[1]= colr[1];
          oldRGBA[2]= colr[2];
          glColor3fv(oldRGBA);
        }
        colr += colrsiz;
        /* NOTE: It makes no difference in local performance if the coords
           are converted to floats before being stored in the display list */
        glVertex3fv(xyz+base);
        glVertex3fv(xyz+base+3);
        glVertex3fv(xyz+base+6);
        base += 9;
      }
    }
  }
  glEnd();
  /* turn off diffuse light not connected to any light source */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, def_ambi);
  glEnable(GL_LIGHT0);
  if(do_alpha) {
    glDisable(GL_BLEND);
  }
  CHEK_ERROR("yglTarrayEmit");
}

#undef USE_ELEMENTS
#undef TVDBG
#undef TVDBG_BASE

#ifndef USE_ELEMENTS
void yglTvarray(long do_alpha, long cpervrt, long ntri, unsigned int *ptndx, float *xyz, 
                float *norm, float *colr)
{
  long i, ndx;
#ifdef TVDBG_BASE
  long ndmin=2000000000, ndmax=0;
  double nrmin=1.0e150, nrmax=-1.0e150, ntmp;
  double xmin=1.0e150, xmax=-1.0e150;
  puts("--- entering yglTvarray");
#endif

  if(do_alpha && !alpha_pass) return;
  if(!do_alpha && alpha_pass) return;
  yglUpdateProperties();

#ifdef TVDBG_BASE
  printf("Expanding arrays\nxyz pointer is %x, norm pointer is %x, colr pointer is %x\n", xyz, norm, colr);
#endif
  /* NOTE: there either must be a color per vertex or a single color
     for all triangles. There is NO support for a color per triangle.
     NOTE: for each enabled "array", a corresponding data vector 
     must be defined and have data filled in before a call 
     to glDrawArrays, glDrawElements, etc. because that call 
     will use data from each enabled array. */
  if(!cpervrt) {
    if(do_alpha) {
#ifdef TVDBG_BASE
printf("alpha, single color\n");
#endif
      glColor4fv(colr);
    } else {
#ifdef TVDBG_BASE
printf("NO alpha, single color\n");
#endif
      glColor3fv(colr);
    }
    glBegin(GL_TRIANGLES);
    for(i= 0; i < 3*ntri; i++) {
      glNormal3fv(norm+3*ptndx[i]);
      glVertex3fv(xyz+3*ptndx[i]);
    }
    glEnd();
  } else {
   if(do_alpha) {
#ifdef TVDBG_BASE
printf("alpha, color per vertex\n");
for(i= 0; i < 3*ntri; i++) {
  ndx= ptndx[i];
  if(xyz[3*ndx] != xyz[3*ndx]) printf("xyz is a NAN for i=%d\n", i);
  if(xyz[3*ndx+1] != xyz[3*ndx+1]) printf("xyz is a NAN for i=%d\n", i);
  if(xyz[3*ndx+2] != xyz[3*ndx+2]) printf("xyz is a NAN for i=%d\n", i);
  if(norm[3*ndx] != norm[3*ndx]) printf("norm is a NAN for i=%d\n", i);
  if(norm[3*ndx+1] != norm[3*ndx+1]) printf("norm is a NAN for i=%d\n", i);
  if(norm[3*ndx+2] != norm[3*ndx+2]) printf("norm is a NAN for i=%d\n", i);
  if(colr[4*ndx] != colr[4*ndx]) printf("colr is a NAN for i=%d\n", i);
  if(colr[4*ndx+1] != colr[4*ndx+1]) printf("colr is a NAN for i=%d\n", i);
  if(colr[4*ndx+2] != colr[4*ndx+2]) printf("colr is a NAN for i=%d\n", i);
  if(colr[4*ndx+3] != colr[4*ndx+3]) printf("colr is a NAN for i=%d\n", i);

  if(xyz[3*ndx]*xyz[3*ndx] > 1.0e30) printf("xyz is too big for i=%d\n", i);
  if(xyz[3*ndx+1]*xyz[3*ndx+1] > 1.0e30) printf("xyz is too big for i=%d\n", i);
  if(xyz[3*ndx+2]*xyz[3*ndx+2] > 1.0e30) printf("xyz is too big for i=%d\n", i);
  if(colr[4*ndx]*colr[4*ndx] > 1.0e30) printf("colr is too big for i=%d\n", i);
  if(colr[4*ndx+1]*colr[4*ndx+1] > 1.0e30) printf("colr is too big for i=%d\n", i);
  if(colr[4*ndx+2]*colr[4*ndx+2] > 1.0e30) printf("colr is too big for i=%d\n", i);
  if(colr[4*ndx+3]*colr[4*ndx+3] > 1.0e30) printf("colr is too big for i=%d\n", i);

  nrmin= norm[3*ndx]*norm[3*ndx] + norm[3*ndx+1]*norm[3*ndx+1] + norm[3*ndx+2]*norm[3*ndx+2];
  if(nrmin > 1.0e3) printf("norm is too big for i=%d\n", i);
  if(nrmin < 1.0e-3) printf("norm is too small for i=%d\n", i);
}
#endif
      glBegin(GL_TRIANGLES);
      for(i= 0; i < 3*ntri; i++) {
        ndx= ptndx[i];
#ifdef TVDBG
    if(ndx < ndmin) ndmin= ndx;
    if(ndx > ndmax) ndmax= ndx;
    if(xyz[3*ndx] < xmin) xmin= xyz[3*ndx];
    if(xyz[3*ndx] > xmax) xmax= xyz[3*ndx];
    if(xyz[3*ndx+1] < xmin) xmin= xyz[3*ndx+1];
    if(xyz[3*ndx+1] > xmax) xmax= xyz[3*ndx+1];
    if(xyz[3*ndx+2] < xmin) xmin= xyz[3*ndx+2];
    if(xyz[3*ndx+2] > xmax) xmax= xyz[3*ndx+2];
    ntmp= norm[3*ndx]*norm[3*ndx]+norm[3*ndx+1]*norm[3*ndx+1]+norm[3*ndx+2]*norm[3*ndx+2];
    if(ntmp < nrmin) nrmin= ntmp;
    if(ntmp > nrmax) nrmax= ntmp;
#endif
        glColor4fv(colr+4*ndx);
        glNormal3fv(norm+3*ndx);
        glVertex3fv(xyz+3*ndx);
      }
      glEnd();
#ifdef TVDBG
  printf("minimum index is %ld and max is %ld\n", ndmin, ndmax);
  printf("minimum x is %e and max is %e\n", xmin, xmax);
  printf("minimum norm is %e and max is %e\n", nrmin, nrmax);
#endif
    } else {
#ifdef TVDBG_BASE
printf("NO alpha, color per vertex\n");
for(i= 0; i < 3*ntri; i++) {
  ndx= ptndx[i];
  if(xyz[3*ndx] != xyz[3*ndx]) printf("xyz is a NAN for i=%d\n", i);
  if(xyz[3*ndx+1] != xyz[3*ndx+1]) printf("xyz is a NAN for i=%d\n", i);
  if(xyz[3*ndx+2] != xyz[3*ndx+2]) printf("xyz is a NAN for i=%d\n", i);
  if(norm[3*ndx] != norm[3*ndx]) printf("norm is a NAN for i=%d\n", i);
  if(norm[3*ndx+1] != norm[3*ndx+1]) printf("norm is a NAN for i=%d\n", i);
  if(norm[3*ndx+2] != norm[3*ndx+2]) printf("norm is a NAN for i=%d\n", i);
  if(colr[3*ndx] != colr[3*ndx]) printf("colr is a NAN for i=%d\n", i);
  if(colr[3*ndx+1] != colr[3*ndx+1]) printf("colr is a NAN for i=%d\n", i);
  if(colr[3*ndx+2] != colr[3*ndx+2]) printf("colr is a NAN for i=%d\n", i);

  if(xyz[3*ndx]*xyz[3*ndx] > 1.0e30) printf("xyz is too big for i=%d\n", i);
  if(xyz[3*ndx+1]*xyz[3*ndx+1] > 1.0e30) printf("xyz is too big for i=%d\n", i);
  if(xyz[3*ndx+2]*xyz[3*ndx+2] > 1.0e30) printf("xyz is too big for i=%d\n", i);
  if(norm[3*ndx]*norm[3*ndx] > 1.0e30) printf("norm is too big for i=%d\n", i);
  if(norm[3*ndx+1]*norm[3*ndx+1] > 1.0e30) printf("norm is too big for i=%d\n", i);
  if(norm[3*ndx+2]*norm[3*ndx+2] > 1.0e30) printf("norm is too big for i=%d\n", i);
  if(colr[3*ndx]*colr[3*ndx] > 1.0e30) printf("colr is too big for i=%d\n", i);
  if(colr[3*ndx+1]*colr[3*ndx+1] > 1.0e30) printf("colr is too big for i=%d\n", i);
  if(colr[3*ndx+2]*colr[3*ndx+2] > 1.0e30) printf("colr is too big for i=%d\n", i);

  if(norm[3*ndx]*norm[3*ndx] < 1.0e-3) printf("norm is too small for i=%d\n", i);
  if(norm[3*ndx+1]*norm[3*ndx+1] < 1.0e-3) printf("norm is too small for i=%d\n", i);
  if(norm[3*ndx+2]*norm[3*ndx+2] < 1.0e-3) printf("norm is too small for i=%d\n", i);
}
#endif
      glBegin(GL_TRIANGLES);
      for(i= 0; i < 3*ntri; i++) {
        ndx= ptndx[i];
#ifdef TVDBG
    if(ndx < ndmin) ndmin= ndx;
    if(ndx > ndmax) ndmax= ndx;
    if(xyz[3*ndx] < xmin) xmin= xyz[3*ndx];
    if(xyz[3*ndx] > xmax) xmax= xyz[3*ndx];
    if(xyz[3*ndx+1] < xmin) xmin= xyz[3*ndx+1];
    if(xyz[3*ndx+1] > xmax) xmax= xyz[3*ndx+1];
    if(xyz[3*ndx+2] < xmin) xmin= xyz[3*ndx+2];
    if(xyz[3*ndx+2] > xmax) xmax= xyz[3*ndx+2];
    ntmp= norm[3*ndx]*norm[3*ndx]+norm[3*ndx+1]*norm[3*ndx+1]+norm[3*ndx+2]*norm[3*ndx+2];
    if(ntmp < nrmin) nrmin= ntmp;
    if(ntmp > nrmax) nrmax= ntmp;
#endif
        glColor3fv(colr+3*ndx);
        glNormal3fv(norm+3*ndx);
        glVertex3fv(xyz+3*ndx);
      }
      glEnd();
#ifdef TVDBG
  printf("minimum index is %ld and max is %ld\n", ndmin, ndmax);
  printf("minimum x is %e and max is %e\n", xmin, xmax);
  printf("minimum norm is %e and max is %e\n", nrmin, nrmax);
#endif
    }
  }
  
  CHEK_ERROR("yglTvarray");
}
#else
void yglTvarray(long do_alpha, long cpervrt, long ntri, unsigned int *ptndx, float *xyz, 
                float *norm, float *colr)
{
  long i, ndx;
#ifdef TVDBG
  long ndmin=2000000000, ndmax=0;
  double nrmin=1.0e150, nrmax=-1.0e150, ntmp;
  double xmin=1.0e150, xmax=-1.0e150;
  puts("--- entering yglTvarray");
#endif

  if(do_alpha && !alpha_pass) return;
  if(!do_alpha && alpha_pass) return;
  yglUpdateProperties();

#ifdef TVDBG
  printf("Using EnableClientState\nxyz pointer is %x, norm pointer is %x, colr pointer is %x\n", xyz, norm, colr);
#endif
  /* NOTE: there either must be a color per vertex or a single color
     for all triangles. There is NO support for a color per triangle.
     NOTE: for each enabled "array", a corresponding data vector 
     must be defined and have data filled in before a call 
     to glDrawArrays, glDrawElements, etc. because that call 
     will use data from each enabled array. */
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_FLOAT, 0, norm);
  glVertexPointer(3, GL_FLOAT, 0, xyz);
  if(cpervrt) {
    glEnableClientState(GL_COLOR_ARRAY);
    if(do_alpha) {
      glColorPointer(4, GL_FLOAT, 0, colr);
    } else {
      glColorPointer(3, GL_FLOAT, 0, colr);
    }
  } else {
    if(do_alpha) {
      glColor4fv(colr);
    } else {
      glColor3fv(colr);
    }
  }
  glDrawElements(GL_TRIANGLES, 3*ntri, GL_UNSIGNED_INT, ptndx);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  CHEK_ERROR("yglTvarray");
}
#endif

void yglTivarray(long ntri, unsigned int *ptndx, void *ileave)
{
  /* draw an array of triangles */
  if(ntri <= 0) return;
  if(alpha_pass) return;

  /* use smooth shading */
  yglSetShade(1);
  yglUpdateProperties();

  /* NOTE: This function uses an interleaved vertex array. All of
     the data must have been stored into the array before this 
     function is called.
     NOTE: there isn't a mode with 3 color components, normals, 
     and vertices.
  */
#undef USE_INTERLEAVE
#ifdef USE_INTERLEAVE
  glInterleavedArrays(GL_C4F_N3F_V3F, 0, ileave);
  glDrawElements(GL_TRIANGLES, 3*ntri, GL_UNSIGNED_INT, ptndx);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_EDGE_FLAG_ARRAY);
  glDisableClientState(GL_INDEX_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
#else
  {
    float *now= (float *)ileave;
    long i;
    glBegin(GL_TRIANGLES);
    for(i= 0; i < ntri; i++) {
      glColor4fv(now);
      glNormal3fv(now+4);
      glVertex3fv(now+7);
      glColor4fv(now+10);
      glNormal3fv(now+14);
      glVertex3fv(now+17);
      glColor4fv(now+20);
      glNormal3fv(now+24);
      glVertex3fv(now+27);
      now += 30;
    }
    glEnd();
  }
#endif
  CHEK_ERROR("yglTivarray");
}

void yglQarray(long smooth, long nquad, float *xyz, float *norm, 
               float *colr, long edge, long cpervrt)
{
  long i, base;
  float old_red= -1.0;
  float old_green= -1.0;
  float old_blue= -1.0;

  /* draw an array of quadrilaterals */
  if(nquad <= 0) return;
  if(alpha_pass) return;
  if(smooth) {
    /* use smooth shading */
    yglSetShade(1);
  } else {
    /* use flat shading */
    yglSetShade(0);
  }
  yglUpdateProperties();

  base= 0;
  glBegin(GL_QUADS);
  if(cpervrt) {
    for(i= 0; i < nquad; i++) {
      if(smooth) {
        glColor3fv(colr+12*i);
        glNormal3fv(norm+base);
        glVertex3fv(xyz+base);
        glColor3fv(colr+12*i+3);
        glNormal3fv(norm+base+3);
        glVertex3fv(xyz+base+3);
        glColor3fv(colr+12*i+6);
        glNormal3fv(norm+base+6);
        glVertex3fv(xyz+base+6);
        glColor3fv(colr+12*i+9);
        glNormal3fv(norm+base+9);
        glVertex3fv(xyz+base+9);
      } else {
        glColor3fv(colr+12*i);
        glNormal3fv(norm+3*i);
        glVertex3fv(xyz+base);
        glColor3fv(colr+12*i+3);
        glVertex3fv(xyz+base+3);
        glColor3fv(colr+12*i+6);
        glVertex3fv(xyz+base+6);
        glColor3fv(colr+12*i+9);
        glVertex3fv(xyz+base+9);
      }
      base += 12;
    }
  } else {
    for(i= 0; i < nquad; i++) {
      if(colr[3*i] != old_red || colr[3*i+1] != old_green || 
         colr[3*i+2] != old_blue) {
        glColor3fv(colr+3*i);
        old_red= colr[3*i];
        old_green= colr[3*i+1];
        old_blue= colr[3*i+2];
      }
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      if(smooth) {
        glNormal3fv(norm+base);
        glVertex3fv(xyz+base);
        glNormal3fv(norm+base+3);
        glVertex3fv(xyz+base+3);
        glNormal3fv(norm+base+6);
        glVertex3fv(xyz+base+6);
        glNormal3fv(norm+base+9);
        glVertex3fv(xyz+base+9);
      } else {
        glNormal3fv(norm+3*i);
        glVertex3fv(xyz+base);
        glVertex3fv(xyz+base+3);
        glVertex3fv(xyz+base+6);
        glVertex3fv(xyz+base+9);
      }
      base += 12;
    }
  }
  glEnd();
  CHEK_ERROR("yglQarrayNoArr");
}

void yglQarrayAlpha(long smooth, long nquad, float *xyz, float *norm, 
                    float *colr, long edge, long cpervrt)
{
  long i, base;
  float old_red= -1.0;
  float old_green= -1.0;
  float old_blue= -1.0;
  float old_alpha= -1.0;
  long colrsiz= 4;

  /* draw an array of quadrilaterals */
  if(nquad <= 0) return;
  if(!alpha_pass) return;
  if(smooth) {
    /* use smooth shading */
    yglSetShade(1);
  } else {
    /* use flat shading */
    yglSetShade(0);
  }
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);
  yglUpdateProperties();

  base= 0;
  glBegin(GL_QUADS);
  if(cpervrt) {
    for(i= 0; i < nquad; i++) {
      if(smooth) {
        glColor3fv(colr+12*i);
        glNormal3fv(norm+base);
        glVertex3fv(xyz+base);
        glColor3fv(colr+12*i+3);
        glNormal3fv(norm+base+3);
        glVertex3fv(xyz+base+3);
        glColor3fv(colr+12*i+6);
        glNormal3fv(norm+base+6);
        glVertex3fv(xyz+base+6);
        glColor3fv(colr+12*i+9);
        glNormal3fv(norm+base+9);
        glVertex3fv(xyz+base+9);
      } else {
        glColor3fv(colr+12*i);
        glNormal3fv(norm+3*i);
        glVertex3fv(xyz+base);
        glColor3fv(colr+12*i+3);
        glVertex3fv(xyz+base+3);
        glColor3fv(colr+12*i+6);
        glVertex3fv(xyz+base+6);
        glColor3fv(colr+12*i+9);
        glVertex3fv(xyz+base+9);
      }
      base += 12;
    }
  } else {
    for(i= 0; i < nquad; i++) {
      if(colr[colrsiz*i] != old_red || colr[colrsiz*i+1] != old_green || 
         colr[colrsiz*i+2] != old_blue || colr[colrsiz*i+2] != old_alpha) {
        glColor3fv(colr+colrsiz*i);
        old_red= colr[colrsiz*i];
        old_green= colr[colrsiz*i+1];
        old_blue= colr[colrsiz*i+2];
        old_alpha= colr[colrsiz*i+2];
      }
      /* NOTE: It makes no difference in local performance if the coords
         are converted to floats before being stored in the display list */
      if(smooth) {
        glNormal3fv(norm+base);
        glVertex3fv(xyz+base);
        glNormal3fv(norm+base+3);
        glVertex3fv(xyz+base+3);
        glNormal3fv(norm+base+6);
        glVertex3fv(xyz+base+6);
        glNormal3fv(norm+base+9);
        glVertex3fv(xyz+base+9);
      } else {
        glNormal3fv(norm+3*i);
        glVertex3fv(xyz+base);
        glVertex3fv(xyz+base+3);
        glVertex3fv(xyz+base+6);
        glVertex3fv(xyz+base+9);
      }
      base += 12;
    }
  }
  glEnd();
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  CHEK_ERROR("yglQarrayNoArr");
}
