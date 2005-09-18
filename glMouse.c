/*
 * $Id: glMouse.c,v 1.1 2005-09-18 22:07:46 dhmunro Exp $
 * The functions in this file respond to mouse movements in 
 * an OpenGL window.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include "ydata.h"

#include "glMouse.h"
#include "glcode.h"
#include "glfunc.h"
#include "glWrappers.h"

/* #include <GL/glu.h> */
#include <math.h>

#define SMALL_RADIUS 1.0e-2

static MsPoint ms_down, ms_last;
int ygl_use_mouse_move= 1;
double ygl_ms_mov_val= 4;

static double ygl_view0[3], ygl_eye0[3], ygl_up0[3], ygl_center0[3], ygl_viewdist0;

/* ----------- TEST ----------- */
#undef DBG_MOUSE
#ifdef DBG_MOUSE
#include "play.h"
#include <stdio.h>
static char msg_buf[120];
#endif
/* ----------- TEST ----------- */

void DoLButtonDown(int flags, int x, int y) 
{
  ms_down.x= x;
  ms_down.y= y;
  ms_last= ms_down;
  ygl_view0[0]= glCurrWin3d->view[0];
  ygl_view0[1]= glCurrWin3d->view[1];
  ygl_view0[2]= glCurrWin3d->view[2];
  ygl_eye0[0]= glCurrWin3d->eye[0];
  ygl_eye0[1]= glCurrWin3d->eye[1];
  ygl_eye0[2]= glCurrWin3d->eye[2];
  ygl_up0[0]= glCurrWin3d->up[0];
  ygl_up0[1]= glCurrWin3d->up[1];
  ygl_up0[2]= glCurrWin3d->up[2];
  ygl_center0[0]= glCurrWin3d->center[0];
  ygl_center0[1]= glCurrWin3d->center[1];
  ygl_center0[2]= glCurrWin3d->center[2];
  ygl_viewdist0= glCurrWin3d->viewdist;
}

void DoLButtonUp(int flags, int x, int y, glWinProp *theWin3d) 
{
  /* carry out rotation, even for a small motion */
  new_mouse_pos(flags, x, y, 1, theWin3d);
}

static void ygl_set_idler(void);

/* flags=1 for left 2 for middle, 3 for right button */
void new_mouse_pos(int flags, int x, int y, int force, glWinProp *theWin3d)
{
  double distsq;
  int difx, dify;
  int action = glCurrWin3d->cursor_action;
  if (flags != 1) {
    /* for non-left button, do alternate action
     * glCurrWin3d->cursor_action == ACTION_ROTATE means
     *   ROTATE left, ZOOM right, PAN middle
     *   (also right is shift-left, middle is control-left)
     * left button always does glCurrWin3d->cursor_action
     *   this is swapped with ACTION_ROTATE for other buttons */
    if (action == ACTION_ROTATE) {
      if (flags == 3) action = ACTION_ZOOM;
      else action = ACTION_PAN;
    } else if (action == ACTION_ZOOM) {
      if (flags == 3) action = ACTION_ROTATE;
      else action = ACTION_PAN;
    } else {
      if (flags == 3) action = ACTION_ZOOM;
      else action = ACTION_ROTATE;
    }
  }

  /* has the mouse moved far enough to do something? */
  difx= x-ms_last.x;
  dify= y-ms_last.y;
  distsq= difx*difx+dify*dify;
  /* No need to do anything if the mouse hasn't moved since last time */
  if(distsq <= 0.0) return;
  /* Handle small residual mouse motion on button up (i.e. when force
     is turned on) */
  if(!force && (distsq < ygl_ms_mov_val*ygl_ms_mov_val) ) return;
#ifdef DBG_MOUSE
  sprintf(msg_buf, "distance moved squared is %e\n", distsq);
  p_stdout(msg_buf);
  sprintf(msg_buf, "x=%d, y=%d, ms_last.x=%d, ms_last.y=%d\n",
          x, y, ms_last.x, ms_last.y);
  p_stdout(msg_buf);
#endif

  /* The mouse has moved significantly since the last point
     or the transform must be updated even for a small move.
     Call the proper action depending on the current
     effect of mouse motions.
  */
  if (action == ACTION_ZOOM) {
    yglMouseZoom(x, y);
  } else if (action == ACTION_ROTATE) {
    yglMouseRot(x, y);
  } else {
    /* action == ACTION_PAN */
    yglMousePan(x, y);
  }

  ms_last.x= x;
  ms_last.y= y;
  if (ygl_use_mouse_move) {
    /* indicate that this window needs to be redrawn */
    theWin3d->dirty= 1;
#if 0
#ifdef IMMEDIATE_REDRAW
    yglDraw3d(theWin3d);
#else
    static long i_draw3_changes = -1;
    if (i_draw3_changes < 0)
      i_draw3_changes = Globalize("_draw3_changes", 0L);
    if (globTab[i_draw3_changes].ops == &dataBlockSym) {
      DataBlock *db = globTab[i_draw3_changes].value.db;
      globTab[i_draw3_changes].ops = &intScalar;
      Unref(db);
    } else if (globTab[i_draw3_changes].ops != &intScalar) {
      globTab[i_draw3_changes].ops = &intScalar;
    }
    globTab[i_draw3_changes].value.i = 1;
    ygl_set_idler();
#endif
#endif
  }
}

static void
ygl_set_idler(void)
{
  extern Function *y_idler_function;
  Function *idler;
  static long i_draw3_idler = -1;
  if (i_draw3_idler < 0) i_draw3_idler = Globalize("_draw3_idler", 0L);
  ASSERT( (globTab[i_draw3_idler].ops == &dataBlockSym &&
      globTab[i_draw3_idler].value.db->ops == &functionOps), 
      "_draw3_idler is not a function")
  idler = (Function *)globTab[i_draw3_idler].value.db;
  y_idler_function = Ref(idler);
}

void yglMouseRot(int x, int y)
{
  double rad, xx, yy;
  double r, csph, snph, snth, csth, r2, csph2, snph2, snth2, csth2;
  double p1[3], p2[3], b[3], aa[3], a[3], w[3];
  double amag, dot12, dotb2, dotu1, dotua, dotub, dotv1, dotva, dotvb;
  double xcen, ycen, dotuv, len, p10, p11, p12, p20, p21, p22, normv;

  /* Compute a new rotation angle based on the offset from
     the original point. 
     The assumption is that the mouse grabs the surface of a 
     virtual trackball and drags it around.
	 NOTE: the trackball is the size of the window, not the size
	 of the object being displayed. 
     A (theta,phi) on the sphere is computed from the 
     (x,y) of the mouse click. The first step is to get
     a cylindrical coord pair, moving the radius in so that
     it is no more than the shortest half-axis of the window.
     The value of theta must lie between zero and pi/2 because
     the mouse can only grab the front of the trackball.
  */
  /* Treat a return to the point where the button went down as a 
     special case and just restore the original viewing
	 transform. */
  if(x == ms_down.x && y == ms_down.y) {
    glCurrWin3d->view[0]= ygl_view0[0];
    glCurrWin3d->view[1]= ygl_view0[1];
    glCurrWin3d->view[2]= ygl_view0[2];
	glCurrWin3d->up[0]= ygl_up0[0];
	glCurrWin3d->up[1]= ygl_up0[1];
	glCurrWin3d->up[2]= ygl_up0[2];
	glCurrWin3d->eye[0]= ygl_eye0[0];
	glCurrWin3d->eye[1]= ygl_eye0[1];
	glCurrWin3d->eye[2]= ygl_eye0[2];
	/* NOTE: the center and viewing distance don't change for
	   rotations */
	return;
  }
  xcen= glCurrWin3d->width/2.0;
  ycen= glCurrWin3d->hite/2.0;
  rad= xcen;
  if(ycen < xcen) rad= ycen;
  xx= ms_down.x-xcen;
  yy= ms_down.y-ycen;
  r= sqrt(xx*xx+yy*yy);
  if(r < SMALL_RADIUS) {
    csph= 1.0;
    snph= 0.0;
  } else {
    csph= xx/r;
    snph= yy/r;
  }
  if(r > rad) r= rad;
  snth= r/rad;
  csth= sqrt(1.0-snth*snth);
  xx= x-xcen;
  yy= y-ycen;
  r2= sqrt(xx*xx+yy*yy);
  if(r2 < SMALL_RADIUS) {
    csph2= 1.0;
    snph2= 0.0;
  } else {
    csph2= xx/r2;
    snph2= yy/r2;
  }
  if(r2 > rad) r2= rad;
  snth2= r2/rad;
  csth2= sqrt(1.0-snth2*snth2);

  /* Let P1 be a unit vector pointing to (theta,phi) and 
     P2 be a unit vector pointing to (theta2,phi2).
     Let A be a vector pointing along the rotation axis.
     The axis is in the direction of P1-cross-P2.
     Let B be the cross product of A and P1.
     The rotation takes P! into P2 and does not alter A.
     In a coordinate system with P1, B, and A as axes,
     the rotation takes a vector S into
     S' = S.P1 P2 + S.B (P1.P2 B - B.P2 P1) + S.A A
     P2, B, and A are known in terms of world coords
     P1 = sin(theta1)*cos(phi1) W + sin(theta1)*sin(phi1) U + cos(theta1) V
     where V is the viewing direction, U is the up vector,
     and W is U-cross-V. There is a similar expression for
     P2 and A and B can be computed from P1 and P2.
  */
  p10= csph*snth;
  p11= snph*snth;
  p12= csth;
  p20= csph2*snth2;
  p21= snph2*snth2;
  p22= csth2;
  /* the vectors given above are in a coordinate system where
     U is the y-axis, V is the z-axis, and W is the x-axis.
     transform to global coordinates */
  w[0]= ygl_view0[1]*ygl_up0[2]-ygl_view0[2]*ygl_up0[1];
  w[1]= ygl_view0[2]*ygl_up0[0]-ygl_view0[0]*ygl_up0[2];
  w[2]= ygl_view0[0]*ygl_up0[1]-ygl_view0[1]*ygl_up0[0];
  p1[0]= p10*w[0]+p11*ygl_up0[0]+p12*ygl_view0[0];
  p1[1]= p10*w[1]+p11*ygl_up0[1]+p12*ygl_view0[1];
  p1[2]= p10*w[2]+p11*ygl_up0[2]+p12*ygl_view0[2];
  p2[0]= p20*w[0]+p21*ygl_up0[0]+p22*ygl_view0[0];
  p2[1]= p20*w[1]+p21*ygl_up0[1]+p22*ygl_view0[1];
  p2[2]= p20*w[2]+p21*ygl_up0[2]+p22*ygl_view0[2];
  aa[0]= p1[1]*p2[2]-p1[2]*p2[1];
  aa[1]= p1[2]*p2[0]-p1[0]*p2[2];
  aa[2]= p1[0]*p2[1]-p1[1]*p2[0];
  amag= sqrt(aa[0]*aa[0]+aa[1]*aa[1]+aa[2]*aa[2]+1.0e-20);
  a[0]= aa[0]/amag;
  a[1]= aa[1]/amag;
  a[2]= aa[2]/amag;
  b[0]= a[1]*p1[2]-a[2]*p1[1];
  b[1]= a[2]*p1[0]-a[0]*p1[2];
  b[2]= a[0]*p1[1]-a[1]*p1[0];
  dot12= p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
  dotb2= b[0]*p2[0]+b[1]*p2[1]+b[2]*p2[2];
  if(dot12 < -1.0) dot12= -1.0;
  if(dot12 > 1.0) dot12= 1.0;
  if(dotb2 < -1.0) dotb2= -1.0;
  if(dotb2 > 1.0) dotb2= 1.0;
	
  /* compute the rotated view vector */
  dotv1= (ygl_view0[0]*p1[0]+ygl_view0[1]*p1[1]+ygl_view0[2]*p1[2]);
  dotvb= (ygl_view0[0]*b[0] +ygl_view0[1]*b[1] +ygl_view0[2]*b[2]);
  dotva= (ygl_view0[0]*a[0] +ygl_view0[1]*a[1] +ygl_view0[2]*a[2]);
  if(dotv1 < -1.0) dotv1= -1.0;
  if(dotv1 > 1.0) dotv1= 1.0;
  if(dotvb < -1.0) dotvb= -1.0;
  if(dotvb > 1.0) dotvb= 1.0;
  if(dotva < -1.0) dotva= -1.0;
  if(dotva > 1.0) dotva= 1.0;
  glCurrWin3d->view[0]= dotv1*p2[0]+dotvb*(dot12*b[0]-dotb2*p1[0])+dotva*a[0];
  glCurrWin3d->view[1]= dotv1*p2[1]+dotvb*(dot12*b[1]-dotb2*p1[1])+dotva*a[1];
  glCurrWin3d->view[2]= dotv1*p2[2]+dotvb*(dot12*b[2]-dotb2*p1[2])+dotva*a[2];
  normv= sqrt(glCurrWin3d->view[0]*glCurrWin3d->view[0]+glCurrWin3d->view[1]*glCurrWin3d->view[1]
	     +glCurrWin3d->view[2]*glCurrWin3d->view[2]+1.0e-20);
  if(normv < 0.99 || normv > 1.01) {
    glCurrWin3d->view[0] /= normv;
    glCurrWin3d->view[1] /= normv;
    glCurrWin3d->view[2] /= normv;
  }
  /* compute the rotated up vector */
  dotu1= (ygl_up0[0]*p1[0]+ygl_up0[1]*p1[1]+ygl_up0[2]*p1[2]);
  dotub= (ygl_up0[0]*b[0] +ygl_up0[1]*b[1] +ygl_up0[2]*b[2]);
  dotua= (ygl_up0[0]*a[0] +ygl_up0[1]*a[1] +ygl_up0[2]*a[2]);
  if(dotu1 < -1.0) dotu1= -1.0;
  if(dotu1 > 1.0) dotu1= 1.0;
  if(dotub < -1.0) dotub= -1.0;
  if(dotub > 1.0) dotub= 1.0;
  if(dotua < -1.0) dotua= -1.0;
  if(dotua > 1.0) dotua= 1.0;
  glCurrWin3d->up[0]= dotu1*p2[0]+dotub*(dot12*b[0]-dotb2*p1[0])+dotua*a[0];
  glCurrWin3d->up[1]= dotu1*p2[1]+dotub*(dot12*b[1]-dotb2*p1[1])+dotua*a[1];
  glCurrWin3d->up[2]= dotu1*p2[2]+dotub*(dot12*b[2]-dotb2*p1[2])+dotua*a[2];
  /* make the up vector normal to the new viewing direction */
  dotuv= glCurrWin3d->up[0]*glCurrWin3d->view[0]+glCurrWin3d->up[1]*glCurrWin3d->view[1]
         +glCurrWin3d->up[2]*glCurrWin3d->view[2];
  glCurrWin3d->up[0] -= glCurrWin3d->view[0]*dotuv;
  glCurrWin3d->up[1] -= glCurrWin3d->view[1]*dotuv;
  glCurrWin3d->up[2] -= glCurrWin3d->view[2]*dotuv;
  len= sqrt(glCurrWin3d->up[0]*glCurrWin3d->up[0]+glCurrWin3d->up[1]*glCurrWin3d->up[1]+glCurrWin3d->up[2]*glCurrWin3d->up[2]+1.0e-20);
  glCurrWin3d->up[0] /= len;
  glCurrWin3d->up[1] /= len;
  glCurrWin3d->up[2] /= len;
  /* update the viewing position */
  glCurrWin3d->eye[0]= glCurrWin3d->center[0]+glCurrWin3d->viewdist*glCurrWin3d->view[0];
  glCurrWin3d->eye[1]= glCurrWin3d->center[1]+glCurrWin3d->viewdist*glCurrWin3d->view[1];
  glCurrWin3d->eye[2]= glCurrWin3d->center[2]+glCurrWin3d->viewdist*glCurrWin3d->view[2];
}

void yglMousePan(int x, int y)
{
  double xx, yy, right0, right1, right2, xcen, ycen, imsize;

  /* Translate the viewpoint perpendicular to the 
     viewing direction. 
     The extent of the window at the viewing distance can be determined
     from the field-of-view. The distance the mouse has moved can
     then be used to change the eye position and center of attention.
     xcen and ycen are the half-widths of the screen in pixels
     and imsize is the half-size of the window at the viewing distance
     in world coords.
  */
  xcen= glCurrWin3d->width/2.0;
  ycen= glCurrWin3d->hite/2.0;
  /* put a fudge factor into the image size to get the right
     motion on the screen */
  imsize= 0.42*glCurrWin3d->viewdist*tan(glCurrWin3d->fov*atan(1.0)/45.0);
  /* set the distance moved since the mouse down in pixels */
  xx= x-ms_down.x;
  yy= y-ms_down.y;

  /* ygl_view0 is the unit normal from the center of attention to the viewer.
     it does not change as the mouse is moved.
     ygl_up0 is the unit vector that points up in the window.
     form the unit vector that points right in the window.
  */
  right0= ygl_up0[1]*ygl_view0[2]-ygl_up0[2]*ygl_view0[1];
  right1= ygl_up0[2]*ygl_view0[0]-ygl_up0[0]*ygl_view0[2];
  right2= ygl_up0[0]*ygl_view0[1]-ygl_up0[1]*ygl_view0[0];

  /* move the center of attention by the fractions of the window
     represented by the mouse motions in the x and y windows
     (corresponding to right and ygl_up0)
     If xx is positive (moved right on the screen), the origin 
     needs to move left.
     If yy is positive, the mouse moved down and the origin needs to
     move up.
  */
  glCurrWin3d->center[0]= ygl_center0[0]+(-right0*xx/xcen+ygl_up0[0]*yy/ycen)*imsize;
  glCurrWin3d->center[1]= ygl_center0[1]+(-right1*xx/xcen+ygl_up0[1]*yy/ycen)*imsize;
  glCurrWin3d->center[2]= ygl_center0[2]+(-right2*xx/xcen+ygl_up0[2]*yy/ycen)*imsize;
  /* update the viewing position */
  glCurrWin3d->eye[0]= glCurrWin3d->center[0]+glCurrWin3d->viewdist*glCurrWin3d->view[0];
  glCurrWin3d->eye[1]= glCurrWin3d->center[1]+glCurrWin3d->viewdist*glCurrWin3d->view[1];
  glCurrWin3d->eye[2]= glCurrWin3d->center[2]+glCurrWin3d->viewdist*glCurrWin3d->view[2];
}

void yglMouseZoom(int x, int y)
{
  double yy, ycen, frac, zm;

  /* zooming only considers the up/down mouse motion.
     moving from the bottom of the screen to the top zooms in by a factor of 4.
     moving the other direction zooms out by a factor of 4.
  */
  ycen= glCurrWin3d->hite/2.0;
  yy= y-ms_down.y;
  /* scale to a fraction of the full screen */
  frac= yy/ycen;
  /* this fractions is the base 2 log of the zoom factor */
  zm= pow(2.0, frac);
  /* modify the lookat args */
  glCurrWin3d->eye[0]= glCurrWin3d->center[0]+zm*ygl_viewdist0*glCurrWin3d->view[0];
  glCurrWin3d->eye[1]= glCurrWin3d->center[1]+zm*ygl_viewdist0*glCurrWin3d->view[1];
  glCurrWin3d->eye[2]= glCurrWin3d->center[2]+zm*ygl_viewdist0*glCurrWin3d->view[2];
  glCurrWin3d->viewdist= zm*ygl_viewdist0;
}

void yglMouseMove(int flags, int x, int y, glWinProp *theWin3d) 
{
#ifdef DBG_MOUSE
  p_stdout("mouse moved\n");
#endif
  new_mouse_pos(flags, x, y, 0, theWin3d);
}

void yglMouseFunc3d(long val)
{
  /* Set the effect of mouse motions in the OpenGL window. */
  if (val == 1) glCurrWin3d->cursor_action = ACTION_ROTATE;
  else if (val == 2) glCurrWin3d->cursor_action = ACTION_ZOOM;
  else  if (val == 3) glCurrWin3d->cursor_action = ACTION_PAN;
  /* Ignore if value is bad */
}
