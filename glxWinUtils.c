/*  
 * $Id: glxWinUtils.c,v 1.1.1.1 2005-09-18 22:07:54 dhmunro Exp $
 * Windowing utilities for GLX.
 *
 * The intent here is to expose the policies behind visual selection
 * so you can tune them to fit your application.
 *
 * get_main_visual() selects a visual in the main planes, with the 
 * specified window attributes. get_colormap() uses the Xmu library 
 * to load a shared colormap based on the visual; if a shared 
 * colormap is not available for the visual that was selected 
 * then a private colormap is allocated.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <GL/gl.h>
#include <GL/glx.h>
#include "glxWinUtils.h"

static Atom wm_delete_window = None;

/*
 * Opens a GLX connection. If successful, this routine returns a
 * pointer to the display connection otherwise it returns NULL.
 */
Display *open_glx_connection( char *display)
{
  Display *dpy;
  int error, event;

  /* Open a connection to the X server */
  if ( !(dpy = XOpenDisplay( display)) ) {
    fprintf( stderr, "Cannot open display %s \n", XDisplayName( display));
    return NULL;
  }

  /* Make sure OpenGL is supported */
  if ( !glXQueryExtension( dpy, &error, &event) ) {
    fprintf( stderr, "No glx extension on %s \n", XDisplayName( display));
    return NULL;
  }

  return dpy;
}

/*  
 * Preference is given to the visual that supports the most colors. 
 */
static int score_visual_class(Display *dpy, XVisualInfo *visual)
{
  int red, green, blue;
  int score = 0;

  switch ( visual->class ) {
  case TrueColor:
  case DirectColor:
    glXGetConfig( dpy, visual, GLX_RED_SIZE, &red);
    glXGetConfig( dpy, visual, GLX_GREEN_SIZE, &green);
    glXGetConfig( dpy, visual, GLX_BLUE_SIZE, &blue);
    score= red+green+blue;
    break;
  }

  return( score);
}

/* 
 * Checks the visual attributes. If the specified attributes are not
 * supported then -1 is returned. If all the specified attributes are
 * supported then the visual is scored according to the depth of the
 * specified ancillary buffers, the depth of the color buffers, and
 * the type of visual. The score for the visual is returned.
 */
static int 
score_visual( Display *dpy, WinAttributes *attr, XVisualInfo *vis,
              WinAttributes *vis_attr)
{
  int val; 
  int score = 0;

  vis_attr->ancillary_buffers = 0;
   
  /* Color Index or RGBA? It must match the requested value */
  glXGetConfig( dpy, vis, GLX_RGBA, &vis_attr->color_mode);
  if ( attr->color_mode != WIN_DONTCARE &&
       attr->color_mode != vis_attr->color_mode ) return( -1);
   
  /* Stereo or mono visual? It must match the requested value */
  glXGetConfig( dpy, vis, GLX_STEREO, &vis_attr->stereo);
  if ( attr->stereo != WIN_DONTCARE &&
       attr->stereo != vis_attr->stereo ) return( -1);

  /* 
   * If double buffering is requested then visual must have a back buffer 
   * If single buffering is requested then it cannot have a back buffer   
   */
  glXGetConfig( dpy, vis, GLX_DOUBLEBUFFER, &val);
  vis_attr->ancillary_buffers |= (val) ? (WIN_DOUBLEBUFFER) : (WIN_SINGLEBUFFER);
  if ( (attr->ancillary_buffers & WIN_DOUBLEBUFFER) &&
       !(vis_attr->ancillary_buffers & WIN_DOUBLEBUFFER) ) return( -1);
  if ( (attr->ancillary_buffers & WIN_SINGLEBUFFER)  &&
       !(vis_attr->ancillary_buffers & WIN_SINGLEBUFFER) ) return( -1);
  /* 
   * Yorick uses alpha blending, but it only uses the alpha
   * of the new "fragment". This means there is no need for
   * an alpha buffer.
   */ 

  /* The depth buffer, if requested, is scored according to it's size */
  /* The visual is rejected if the size is zero                       */
  glXGetConfig( dpy, vis, GLX_DEPTH_SIZE, &val);
  if ( val > 0 ) 
    vis_attr->ancillary_buffers |= WIN_DEPTHBUFFER;
  if ( attr->ancillary_buffers & WIN_DEPTHBUFFER ) {
    if ( val <= 0 ) return( -1);
    score += val;
  }

  /* Score the visual class based on number of colors supported */
  score += score_visual_class(dpy, vis);

  return( score);
}


/* 
 * Finds the best main planes visual that supports all the requested attributes. 
 * (i.e., all requested ancillary buffers must be supported. Also 
 * if stereo, double buffering or single buffering is requested 
 * then the visual must support it. Note that, a double buffered visual will 
 * not be returned if single buffering is explicitly requested.) The caller 
 * should specify WIN_DONTCARE if they don't want care about the value of a
 * particular attribute.
 * 
 * If more than one visual satisfies the above criteria then preference
 * is given to the visual that supports the most colors and has the deepest
 * ancillary buffers. (But, as noted above, the number of available colors is 
 * considered more important than the depth of the requested ancillary
 * buffers.)
 */
static 
int get_main_visual( Display *dpy, WinAttributes *attr, XVisualInfo *vis_info, WinAttributes *vis_attr)
{
  int i, status, nvisuals, val, this_score, best_score;
  XVisualInfo *vis_list, *this_vis, *best_vis, sample_vis;
  WinAttributes this_vis_attr, best_vis_attr;

  /* Get list of visuals for this screen */
  sample_vis.screen = DefaultScreen( dpy);
  vis_list = XGetVisualInfo( dpy, VisualScreenMask, &sample_vis, &nvisuals);
  /* Loop through the visuals to find the one that has the best score */
  best_score = -1; best_vis = NULL;
  for ( i = 0; i < nvisuals; i++ ) {
    this_vis = &vis_list[i];
    this_score = 0;

    /* Visual must be supported by GLX */
    if ( glXGetConfig( dpy, this_vis, GLX_USE_GL, &val) ) continue;
    if ( !val ) continue;

    /* Visual must be in main planes which is level 0 */ 
    if ( glXGetConfig( dpy, this_vis, GLX_LEVEL, &val) ) continue;
    if ( val != 0 ) continue;
   
    /* 
     * Score the visual. If any of the requested buffers are not available 
     * the returned score will be negative.
     */
    this_score = score_visual( dpy, attr, this_vis, &this_vis_attr);
    if ( this_score < 0 ) continue;

    if ( this_score > best_score ) {
      best_score = this_score;
      best_vis = this_vis;
      best_vis_attr = this_vis_attr;
    }
  }

  status = GL_FALSE;
  if ( best_vis != NULL ) {
    *vis_info = *best_vis;
    *vis_attr = best_vis_attr;
    vis_attr->level = 0;
    vis_attr->transparent_type = WIN_OPAQUE;
    status = GL_TRUE;
  }

  if ( nvisuals > 0 ) XFree ( vis_list);
  return status;
}



/*
 * Given a visual, choose an appropriate colormap. We use the Xmu 
 * utilities to make sure a standard colormap is loaded. 
 *
 * There are several types of standard colormaps:
 *           
 *   RGB_DEFAULT_MAP: The range of colors may be less than for 
 *       RGB_BEST_MAP. For example on an 8-plane system, a typical 
 *       allocation is 6 reds, 6 greens and 6 blues, yielding 216
 *       uniformly distributed colors. This leaves 40 elements of a 
 *       256 element colormap available for special purpose colors.
 *       The RGB_DEFAULT_MAP is usually available for all main plane 
 *       visuals. 
 *   RGB_BEST_MAP: This may only be available for one visual--usually a
 *       true color or direct color visual. It may support more distinct
 *       colors than RGB_DEFAULT_MAP. 
 *   RGB_RED_MAP, RGB_GREEN_MAP, RGB_BLUE_MAP: These colormaps have all-red,
 *       all-green or all-blue components. This is useful if you want to
 *       make color separated images
 *   RGB_GRAY_MAP: A gray scale colormap.
 *
 * The code below always tries to load the default map. Why? The gray
 * red, green and blue maps, aren't generally desirable for 3D graphics.
 * And the best map--if it is supported on this visual--is probably the 
 * same as the default map. 
 * 
 * If a standard colormap isn't available for this visual, then we create
 * a colormap (Most likely, this will only happen in the case of overlay and 
 * underlay visuals)
 */
static 
Colormap get_colormap( Display *dpy, XVisualInfo *vis_info)
{
  Atom wm_property;
  XStandardColormap *std_cmaps;
  int num_cmaps, i;
  Colormap cmap;

  /* 
   * NOTE: You may want to check for RGB_BEST_MAP first and then, if it isn't
   * available, try to get the the default map
   */

  wm_property = XA_RGB_DEFAULT_MAP;

  if ( XmuLookupStandardColormap( dpy, vis_info->screen, vis_info->visualid,
                                  vis_info->depth, wm_property, False/*replace*/, True/*retain*/) ) {
    if ( XGetRGBColormaps( dpy, RootWindow( dpy, vis_info->screen),
                           &std_cmaps, &num_cmaps, wm_property) ) {
      for (i = 0; i < num_cmaps; i++) {
        if ( std_cmaps[i].visualid == vis_info->visualid ) {
          cmap = std_cmaps[i].colormap;
          XFree( std_cmaps);
          return cmap;
        }
      }
      if ( num_cmaps > 0 ) XFree( std_cmaps);
    }
  }

  return( XCreateColormap( dpy, RootWindow( dpy, vis_info->screen), 
                           vis_info->visual, AllocNone));
}


/* 
 * Finds the best visual in the overlay/underlay/main planes (depending
 * on the level attribute) and a colormap to match.
 */
int get_visual_and_colormap(Display *dpy, int mask, WinAttributes *attr, 
                            XVisualInfo *vis_info, Colormap *cmap)
{
  int status;
  WinAttributes all_attr;
  
  all_attr.color_mode = WIN_DONTCARE;
  all_attr.stereo = WIN_DONTCARE;
  all_attr.ancillary_buffers = 0;
  all_attr.aux_buffers = 0;
  all_attr.level = 0;
  all_attr.transparent_type = WIN_DONTCARE;

  if ( mask & WIN_COLOR_MODE_ATTR )
    all_attr.color_mode = attr->color_mode;
  if ( mask & WIN_STEREO_ATTR )
    all_attr.stereo = attr->stereo;
  if ( mask & WIN_ANCILLARY_BUFFERS_ATTR )
    all_attr.ancillary_buffers = attr->ancillary_buffers;

  status = get_main_visual( dpy, &all_attr, vis_info, attr);
  if(status == GL_FALSE) return status;

  *cmap = get_colormap( dpy, vis_info);
  return GL_TRUE;
}


/* 
 * Creates a window with the specified visual and colormap and parent. The 
 * specified position, size and title are used to set window property hints.
 */
Window create_window(Display *dpy, Window parent_win, char *title_string, 
                     XVisualInfo *vis_info, Colormap cmap, int x, int y,
                     unsigned int width, unsigned int height)
{
  Window win;
  XSizeHints win_hint; /* window manager hints */
  XSetWindowAttributes win_attr;

  win_hint.x = x; win_hint.y = y;
  win_hint.width = width; win_hint.height = height;
  win_hint.flags = PPosition | PSize;

  /* Create the window using the specified visual and colormap */
  /* Have to specify border pixel if visual is not the default */
  win_attr.colormap = cmap;
  win_attr.border_pixel = 0;
  win = XCreateWindow( dpy, parent_win,
                       win_hint.x, win_hint.y, win_hint.width, win_hint.height,
                       0 /*border width*/, vis_info->depth, InputOutput, vis_info->visual,
                       CWColormap | CWBorderPixel, &win_attr);

  /* Now set the window title and size hints */
  XSetStandardProperties( dpy, win, title_string, title_string, None, NULL,
                          0, &win_hint);

  /* Allow clean termination from window manager */
  if ( wm_delete_window == None )
    wm_delete_window = XInternAtom( dpy, "WM_DELETE_WINDOW", False);
  if ( wm_delete_window != None )
    XSetWMProtocols( dpy, win, &wm_delete_window, 1);

  return win;

}

int decode_client_message( XEvent *event)
{
  if (event->xclient.data.l[0] == wm_delete_window)
    return WIN_CLIENT_SHUTDOWN;
  return 0;
}

