/*
 * $Id: glxWinUtils.h,v 1.1.1.1 2005-09-18 22:07:54 dhmunro Exp $
 * DEC_COPYRIGHT
 */
/*
 * HISTORY
 * Log: win_utils.h,v
 * Revision 1.1.2.2  1994/12/21  15:49:04  John_Dennis
 * 	new contributed code from SGI
 * 	[1994/12/21  15:42:51  John_Dennis]
 *
 */

/* Data structure for SERVER_OVERLAY_VISUALS property */

typedef struct {
   VisualID visualid;
   int transparent_type; /* None, TransparentPixel, TransparentMask */
   int value;            /* Pixel value */
   int layer;            /* layer > 0 for overlay planes, < 0 for underlay planes */
} WinOverlayProperty;

/* GLX window attributes */

typedef struct {
    int color_mode;        /* WIN_DONTCARE, WIN_RGBA, or WIN_INDEX */
    int stereo;            /* WIN_DONTCARE, WIN_STEREO, or WIN_MONO */
    int ancillary_buffers; /* OR in desired ancillary buffers */
    int aux_buffers;       /* number of aux buffers desired */
    int level;             /* underlay < 0, overlay > 0, main planes 0 */
    int transparent_type;  /* WIN_DONTCARE, WIN_OPAQUE, or WIN_TRANSPARENT */
} WinAttributes;

#define WIN_COLOR_MODE_ATTR		1
#define WIN_STEREO_ATTR			2
#define WIN_ANCILLARY_BUFFERS_ATTR	4
#define WIN_AUX_BUFFERS_ATTR		8
#define WIN_LEVEL_ATTR			16
#define WIN_TRANSPARENT_TYPE_ATTR	32

#define WIN_DONTCARE		2

#define WIN_INDEX      	 	0  /* Matches value returned by glXGetConfig */
#define WIN_RGBA       	 	1  

#define WIN_MONO      	 	0  /* Matches value returned by glXGetConfig */
#define WIN_STEREO     	 	1  

#define WIN_OPAQUE         	0		
#define WIN_TRANSPARENT 	1

#define WIN_DOUBLEBUFFER       	1
#define WIN_SINGLEBUFFER       	2
#define WIN_ACCUMBUFFER         4
#define WIN_ALPHABUFFER       	8
#define WIN_DEPTHBUFFER        	16
#define WIN_STENCILBUFFER      	32

#define WIN_CLIENT_SHUTDOWN 	1

/* Function prototypes */

Display *open_glx_connection( char *);
int get_visual_and_colormap(Display *, int, WinAttributes *, XVisualInfo *, 
      Colormap *);
Window create_window(Display *, Window, char *, XVisualInfo *, Colormap, int, int, 
      unsigned int, unsigned int);
int decode_client_message( XEvent *);
