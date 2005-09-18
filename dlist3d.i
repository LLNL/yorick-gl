/*
 * $Id: dlist3d.i,v 1.1 2005-09-18 22:07:45 dhmunro Exp $
 * These functions insert objects into the yorgl display,
 * change the viewpoint, play the display list, etc.
 * Users ordinarily call functions in dlist3d.i so those 
 * functions are autoloaded by the yorgl plugin.
 * The compiled functions are connected via glfunc.i, so 
 * require that glfunc.i be "pulled in" any time something in dlist3d.i
 * is called.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

require,"tristruct.i";
require,"glfunc.i";

// by default do not use interleaved vertex arrays (i.e. separate vertex
// arrays for xyz, normal, and color)
use_interleave= 0;

func draw3d(called_as_idler)
/* DOCUMENT draw3d
     Draw the current 3D "cached" and "direct" display lists.
     (Ordinarily triggered automatically when the drawing changes.)
 */
{
  if (_draw3_changes) {
    /* draw "cached" and direct display lists */
    Draw3d;
    _draw3_changes= [];
  }
}

func draw3d_trigger
{
  extern _draw3_changes;
  /* arrange to call draw3 when everything else is finished */
  set_idler, _draw3_idler;
  _draw3_changes= 1;
}

func clear3d(void)
/* DOCUMENT clear3d
     Clear the 3D "cached" and "direct" display lists.
     Clear the screen.
 */
{
  ClearList3d;
}

func clear3d_direct(void)
/* DOCUMENT clear3d_direct
     Clear the 3D "direct" display list.
 */
{
  ClearDirectList3d;
}

func clear3d_cache(void)
/* DOCUMENT clear3d_cache
     Clear the 3D "cached" display list.
 */
{
  ClearCachedList3d;
}

func getback_rgb3d(void)
{
  rgb= array(0.0, 3);
  getbackrgb_raw3d,rgb;
  return rgb;
}

func getcage_rgb3d(void)
{
  rgb= array(0.0, 3);
  getcagergb_raw3d,rgb;
  return rgb;
}

func getgrid_rgb3d(void)
{
  rgb= array(0.0, 3);
  getgridrgb_raw3d,rgb;
  return rgb;
}

func map2color(ntri, cellids, var, vmin=, vmax=)
{
  /* map var(cellids) into RGB colors (where a color
     is a float triple between zero and one).
  */
  if(is_void(vmin)) vmin= min(var);
  if(is_void(vmax)) vmax= max(var);
  colr= array(float, 3, ntri);
  yglMap2ColorRaw3d, gl_ncolr, gl_rr, gl_gg, gl_bb, vmin, 
                   vmax, var, ntri, cellids, colr;
  return colr;
}

func palette3d(name)
/* DOCUMENT palette3d, palette_name
     Select a color palette to be used, for example,
	 to color slices through 3D meshes.
	 Uses the same palettes as the 2D graphics package.
   SEE ALSO: palette
 */
{
  extern gl_rr,gl_gg,gl_bb, gl_ncolr, gl_ctab;

  // NOTE: If name is a fully qualified path, should
  // try to open it and return if that fails.
  // Otherwise should try to find the file along the
  // "palette path"
  is_absolute= strpart(name, 1:1) == "/";
  if(is_absolute) {
    f= open(name, "r");
  } else {
    fnam= Y_SITE+"g/"+name;
    f= open(fnam, "r");
  }
  lin= rdline(f);
  while(1) {
    lin= rdline(f);
    res= strmatch(lin, "ncolors=");
    if(res) break;
  }
  gl_ncolr= 0;
  sread, lin, format="ncolors=%d", gl_ncolr;
  while(1) {
    lin= rdline(f);
    ASSERT, (lin), "EOF in palette file without finding rgb header line";
    lin2= strtok(lin);
    if(lin2(1) == "#") {
      lin2= strtok(lin2(2));
      if(lin2(1) == "r") {
        lin2= strtok(lin2(2));
    	if(lin2(1) == "g") {
    	  lin2= strtok(lin2(2));
    	  if(lin2(1) == "b") break;
    	}
      }
    }
  }
  gl_rr= gl_gg= gl_bb= array(char, gl_ncolr);
  read,f,gl_rr,gl_gg,gl_bb;
  gl_ctab= transpose([gl_rr,gl_gg,gl_bb]);
}

/* ------------------------------------------------------------------------ */

func light3d(ambient=,diffuse=,specular=,spower=,sdir=)
/* DOCUMENT light3d, ambient=a_level, diffuse=d_level,
            specular=s_level, spower=n, sdir=xyz
     Sets lighting properties for 3D shading effects.
     A surface will be shaded according to its orientation
     relative to the viewing direction.

     The ambient level A_LEVEL is a light level (arbitrary units)
     that is added to every surface independent of its orientation.

     The diffuse level D_LEVEL is a light level which is proportional
     to cos(theta), where theta is the angle between the surface
     normal and the viewing direction, so that surfaces directly
     facing the viewer are bright, while surfaces viewed edge on are
     unlit (and surfaces facing away, if drawn, are shaded as if they
     faced the viewer).

     The specular level S_LEVEL is a light level proportional to a high
     power spower=N of 1+cos(alpha), where alpha is the angle between
     the specular reflection angle and the viewing direction.  The light
     source for the calculation of alpha lies in the direction XYZ (a
     3 element vector) in the viewer's coordinate system at infinite
     distance.

   EXAMPLES:
     light3d, diffuse=.1, specular=1., sdir=[0,0,-1]
       (dramatic "tail lighting" effect)
     light3d, diffuse=.5, specular=1., sdir=[1,.5,1]
       (classic "over your right shoulder" lighting)

 */
{
  fn_name= "light3d";
  if (!is_void(ambient)) {
    ASSERT, (dimsof(ambient)(1) = 0), "ambient light level must be scalar in "+fn_name;
    ambient += 0.0;
  } else {
    ambient= 0.2;
  }
  if (!is_void(diffuse)) {
    ASSERT, (dimsof(diffuse)(1) == 0), "diffuse light level must be a scalar"+fn_name;
    diffuse += 0.0;
  } else {
    diffuse= 0.6;
  }
  if (!is_void(specular)) {
    ASSERT, (dimsof(specular)(1) == 0), "specular strength must be a scalar in "+fn_name;
    specular += 0.0;
  } else {
    specular= 1.0;
  }
  if (!is_void(spower)) {
    ASSERT, (dimsof(spower)(1) == 0), "specular power must be a scalar in "+fn_name;
    spower += 0.0;
  } else {
    spower= 2.0;
  }
  if (!is_void(sdir)) {
    dims= dimsof(sdir);
    ASSERT, (dims(1) == 1 && dims(2) == 3), "lighting direction must be 3 vector in "+fn_name;
    sdirval= [sdir(1), sdir(2), sdir(3), 0.0];
  } else {
    sdirval= [0.0, 0.0, 1.0, 0.0];
  }
  set_light3d, ambient, diffuse, specular, spower, sdirval;
  draw3d_trigger;
}


/* ------------------------------------------------------------------------ */

func stdview3d(dummy)
/* DOCUMENT stdview3d
     Sets the 3D viewing transform to a default value that
	 makes the whole scene visible.

   SEE ALSO: lookat3d, prtview3d
 */
{
  box= GetBounds3d();
  if(box(1,1) > box(2,1) || box(1,2) > box(2,2)||
     box(1,3) > box(2,3)) {
    // didn't find any objects, so set default values
    center= [0.0,0.0,0.0];
    delta= [1.0,1.0,1.0];
  } else {
    center= box(avg,);
    delta= abs(box(dif,));
  }
  /* The stdview3 function sets the viewpoint so that all
     objects will be visible, no matter how the scene is
     rotated. The field-of-view is for the width or height
     of the screen (not its diagonal), so sizing is done in
     a plane (y=0 is fine). 
     Take the longest edge of the bounding box and construct
     a square of that size and the same center in the y=0 plane.
     Draw the circumscribing circle (it has radius sqrt(2)/2
     times the longest edge). Draw a tangent line from the 
     viewpoint to the circle using the given field-of-view.
     The distance from the eye to the center of interest is
     the radius of the circle divided by the sine of the fov.
  */
  radius= sqrt(2.0)*0.5*max(delta);
  dist= radius/sin(pi/180.0*get_fov3d());
  eye= center+[dist, 0.0, 0.0];
  up= [0.0, 0.0, 1.0];
  lookat3d,eye,center,up;
}

func prtview3d
/* DOCUMENT prtview3d
     Print the parameters for the current viewing 
     transformation.

   SEE ALSO: lookat3d, stdview3
 */
{
  eye= center= up= array(0.0, 3);
  get_center3d,center;
  get_eye3d,eye;
  get_up3d,up;
  write,"[viewpoint], [center of interest], [up direction] are";
  write,format="[%e,%e,%e], \n", eye(1), eye(2), eye(3);
  write,format="[%e,%e,%e], \n", center(1), center(2), center(3);
  write,format="[%f,%f,%f] \n", up(1), up(2), up(3);
}

func lookat3d(eye,center,up)
/* DOCUMENT lookat3d(eye,center,up)
     Sets the 3D viewing transform so that the viewer is
     at "eye, the center of the view is "center", and 
     "up" points upward.

   SEE ALSO: stdview3, prtview3
 */
{
  fn_name= "lookat3d";
  dim_e= dimsof(eye);
  ASSERT, (dim_e(1) == 1 && dim_e(2) == 3), "eye must be a 3 element vector in "+fn_name;
  dim_c= dimsof(center);
  ASSERT, (dim_c(1) == 1 && dim_c(2) == 3), "center must be a 3 element vector in "+fn_name;
  if (is_void(up)) {
    vw= eye-center;
    crs= [-vw(2), vw(1), 0.0];
    lv= vw(rms);
    lc= crs(rms);
    if(lc < 1.0e-4*lv) {
      up= [0.0, 1.0, 0.0];
    } else {
      up= [vw(2)*crs(3)-vw(3)*crs(2), 
           vw(3)*crs(1)-vw(1)*crs(3), 
           vw(1)*crs(2)-vw(2)*crs(1)];
    }
  }
  /* lookat_raw3d will make compute an up vector that is orthonormal
     to the viewing direction. */
  lookat_raw3d,eye,center,up;
  draw3d_trigger;
}

func get_lims3d(dummy)
/* DOCUMENT get_lims3d
     Compute and return a bounding box containing everything 
     in the "cached" 3D display list.
 */
{
  lim3d= array(0.0, 2, 3);
  res= GetBounds3d(lim3d);
  /* if no objects, will return an empty box */
  return lim3d;
}

func get_normal3d(xyz, nxyz)
/* DOCUMENT get_normal3d(xyz, nxyz)
         or get_normal3d(xyz)

     return 3D normals for polygons with vertices XYZ.  If NXYZ is
     specified, XYZ should be 3-by-sum(nxyz), with NXYZ being the
     list of numbers of vertices for each polygon (as for the plfp
     function).  If NXYZ is not specified, XYZ should be a quadrilateral
     mesh, 3-by-ni-by-nj (as for the plf function).  In the first case,
     the return value is 3-by-numberof(NXYZ); in the second case, the
     return value is 3-by-(ni-1)-by-(nj-1).

     The normals are constructed from the cross product of the lines
     joining the midpoints of two edges which as nearly quarter the
     polygon as possible (the medians for a quadrilateral).  No check
     is made that these not be parallel; the returned "normal" is
     [0,0,0] in that case.  Also, if the polygon vertices are not
     coplanar, the "normal" has no precisely definable meaning.

   SEE ALSO: get_centroid3d
 */
{
  if (is_void(nxyz)) {
    /* if no polygon list is given, assume xyz is 2D mesh */
    /* form normal as cross product of medians */
    m1= xyz(,zcen,dif);
    m2= xyz(,dif,zcen);

  } else {
    /* with polygon list, more elaborate calculation required */
    frst= nxyz(psum)-nxyz+1;

    /* form normal by getting two approximate diameters
     * (reduces to above medians for quads) */
    n2= (nxyz+1)/2;
    zero= frst-1;
    c0= 0.5*(xyz(,zero+1)+xyz(,zero+2));
    i= zero+n2;  /* n2>=2, nxyz>=3 */
    c1= 0.5*(xyz(,i)+xyz(,i+1));
    i= 1+n2/2;
    c2= 0.5*(xyz(,zero+i)+xyz(,zero+i+1));
    i= min(i+n2, nxyz);
    c3= 0.5*(xyz(,zero+i)+xyz(,zero+i%nxyz+1));
    m1= c1 - c0;
    m2= c3 - c2;
  }

  /* poly normal is cross product of two medians (or diameters) */
  normal= m1;
  normal(1,..)= n1= m1(2,..)*m2(3,..) - m1(3,..)*m2(2,..);
  normal(2,..)= n2= m1(3,..)*m2(1,..) - m1(1,..)*m2(3,..);
  normal(3,..)= n3= m1(1,..)*m2(2,..) - m1(2,..)*m2(1,..);
  m1= abs(n1,n2,n3)(-,..);
  m1= m1 + (m1==0.0);
  normal/= m1;

  return normal;
}

func get_centroid3d(xyz, nxyz)
/* DOCUMENT get_centroid3d(xyz, nxyz)
         or get_centroid3d(xyz)

     return 3D centroids for polygons with vertices XYZ.  If NXYZ is
     specified, XYZ should be 3-by-sum(nxyz), with NXYZ being the
     list of numbers of vertices for each polygon (as for the plfp
     function).  If NXYZ is not specified, XYZ should be a quadrilateral
     mesh, 3-by-ni-by-nj (as for the plf function).  In the first case,
     the return value is 3-by-numberof(NXYZ); in the second case, the
     return value is 3-by-(ni-1)-by-(nj-1).

     The centroids are constructed as the mean value of all vertices
     of each polygon.

   SEE ALSO: get_normal3d, get_light3d
 */
{
  if(is_void(xyz)) return [0.0,0.0,0.0];
  if (is_void(nxyz)) {
    /* if no polygon list is given, assume xyz is 2D mesh */
    centroid= xyz(,zcen,zcen);

  } else {
    /* with polygon list, more elaborate calculation required */
    last= nxyz(psum);
    list= histogram(1+last)(1:-1);
    list(1)+= 1;
    list= list(psum);
    centroid= array(0.0, 3, numberof(nxyz));
    centroid(1,)= histogram(list, xyz(1,));
    centroid(2,)= histogram(list, xyz(2,));
    centroid(3,)= histogram(list, xyz(3,));
    centroid/= double(nxyz);
  }

  return centroid;
}

func plpoly3d(nv, xyz, color, norm, draw_edge=)
/* DOCUMENT plpoly3d(nv, xyz, color, norm, draw_edge=)

     Draws a set of filled polygons in 3D using OpenGL.

   SEE ALSO: plcell3d, lookat3d, light3d, win3d
 */
{
  /* This function accepts a wide range of inputs and coerces
     them into a more restricted set that are saved in yorick's
	 display list.
	 This function expects xyz, color, and norm to be doubles,
	 but should function correctly if they happen to be floats.
  */
  fn_name= "plpoly3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3), "xyz must be a 3-by-nverts array in "+fn_name;
  nvert= sum(nv);
  npoly= numberof(nv);
  dimc= dimsof(color);
  if(dimc(1) == 1) {
    ASSERT, (dimc(2) == 3), "color must be a 3-by-ncolor array in "+fn_name;
    /* broadcast color out to all polys */
    color= array(color, npoly);
  } else {
    ASSERT, (dimc(1) == 2 && dimc(2) == 3), "color must be a 3-by-ncolor array in "+fn_name;
    if(dimc(3) == 1) {
      /* broadcast color out to all polys */
      color= array(color(,1), npoly);
    } else {
      ASSERT, (dimc(3) == npoly), "number of colors must be 1 or npoly in "+fn_name;
    }
  }
  if(typeof(color) == "char") color= float(color/255.0);
  /* the user can request outlined edges, not filled polygons */
  do_edge= 0;
  if(!is_void(draw_edge) && draw_edge) do_edge= 1;
  /* Do not light the surface unless normals were supplied. 
     If there are as many normals as vertices, use smooth
     shading */
  do_light= 1;
  if(is_void(norm)) {
    do_light= 0;
    norm= 0;
  } else {
    dimsn= dimsof(norm);
    ASSERT, (dimsn(1) == 2 && dimsn(2) == 3), "norm must be a 3-by-n array in "+fn_name;
    /* normals must be one per poly */
    if(dimsn(3) == npoly) {
      smooth= 0; /* normal per polygon, so use flat shading */
    } else if(dimsn(3) == dimsx(3)) {
      /* normals must be one per vertex */
      smooth= 1; /* normal per vertex, so use smooth shading */
    } else {
      ASSERT, (0),"norm must be a 3-by-npoly or a 3-by-nvert array in "+fn_name;
    }
  }
  /* add the object to the display list */
  polys3d, npoly, nv, xyz, norm, color, do_edge, smooth, do_light;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plcell3d(corners, color)
/* DOCUMENT plcell3d(corners, color)

     Draws a cell array in 3D using OpenGL.

   SEE ALSO: plpoly3d, lookat3d, light3d, win3d
 */
{
  /* This function expects corner, color, and norm to be doubles on input,
     but should function correctly if they happen to be floats.
     No check is made for corners containing 3 co-linear points.
  */
  fn_name= "plcell3d";
  dimc= dimsof(color);
  ASSERT, (dimc(1) == 3), "color array must be 3-by-nx-by-ny in "+fn_name;
  if(dimc(2) == 3) {
    do_alpha= 0;
  } else if(dimc(2) == 4) {
    do_alpha= 1;
  } else {
    ASSERT, (0),"color array must have first dimension 3 or 4 in "+fn_name;
  }
  if(typeof(color) == "char") color= float(color/255.0);
  nx= dimc(3);
  ny= dimc(4);
  norm= get_normal3d(corners, [3]);
  /* add the object to the display list */
  cells3d, nx, ny, corners, norm, color, do_alpha;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plm3d(xyz, color)
/* DOCUMENT plm3d(xyz, color)

     Draws a logically 2D mesh in 3D using OpenGL.

   SEE ALSO: plpoly3d, plcell3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz and color to be doubles,
     but should function correctly if they happen to be floats.
  */
  /* are the arguments compatible lengths? */
  fn_name= "plm3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 3 && dimsx(2) == 3), "mesh vertices must be 3-by-nx-by-ny in "+fn_name;
  nx= dimsx(3);
  ny= dimsx(4);
  is_3_vector,color,"color must be a 3 element vector in "+fn_name;
  if(typeof(color) == "char") color= float(color/255.0);
  /* add the object to the display list */
  plm3d_raw, nx, ny, xyz, color;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plf3d(xyz, color)
/* DOCUMENT plf3d(xyz, color)

     Fills a logically 2D mesh in 3D using OpenGL.
     color is (3 or 4)-by-(nx-1)-by-(ny-1).
     The initial index of color is 3 (for RGB colors) or 4 (for 
     translucent RGBA colors).
     Each zone is filled with a single color.
     The surface does not "reflect" light, so there will not be any shiny highlights.

   SEE ALSO: plm3d, plsurf3d, plcolrsurf3d, plcell3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz and color to be doubles,
     but should function correctly if they happen to be floats.
  */
  /* are the arguments compatible lengths? */
  fn_name= "plf3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 3 && dimsx(2) == 3), "mesh vertices must be 3-by-nx-by-ny in "+fn_name;
  nx= dimsx(3);
  ny= dimsx(4);
  dimc= dimsof(color);
  ASSERT, (dimc(1) == 3), "color must be a 2D array of colors in "+fn_name;
  if(dimc(3) == nx-1 && dimc(4) == ny-1) {
    // proper zone centered colors
  } else {
    ASSERT, (0),"color is not a 3-by-(nx-1)-by-(ny-1) or 4-by-(nx-1)-by-(ny-1) array in "+fn_name;
  }
  if(typeof(color) == "char") color= float(color/255.0);
  // The compiled routine needs to have an alpha value for each color, so set alpha to one
  // if it was not supplied
  if(dimc(2) == 3) {
    colnu= array(float, 4, nx-1, ny-1);
    colnu(1:3,..)= color;
    colnu(4,..)= 1.0f;
    /* add the object to the display list */
    plf3d_raw, nx, ny, xyz, colnu;
  } else {
    /* add the object to the display list */
    plf3d_raw, nx, ny, xyz, color;
  }
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plsurf3d(xyz, norm, color, flip=)
/* DOCUMENT plsurf3d(xyz, norm, color, flip=)

     Draws a lit 2D mesh in 3D using OpenGL.

   SEE ALSO: plcolrsurf3d, plm3d, plpoly3d, plcell3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz, color, and norm to be doubles,
	 but should function correctly if they happen to be floats.
  */
  /* are the arguments compatible lengths? */
  fn_name= "plsurf3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 3 && dimsx(2) == 3), "mesh vertices must be 3-by-nx-by-ny in "+fn_name;
  nx= dimsx(3);
  ny= dimsx(4);
  ASSERT, allof(dimsx == dimsof(norm)),"coordinates and normals must have the same dimensions in "+fn_name;
  dimc= dimsof(color);
  if(dimc(2) == 3) {
    do_alpha= 0;
  } else if(dimc(2) == 4) {
    do_alpha= 1;
  } else {
    ASSERT, (0), "color must have 3 or 4 components in "+fn_name;
  }
  if(typeof(color) == "char") color= float(color/255.0);
  if(!is_void(flip)) {
    dimf= dimsof(flip);
    ASSERT, (dimf(1) == 1 && dimf(2) == 3), "flip must be a 3 element vector in "+fn_name;
    xyz2= xyz*flip;
    norm2= norm*flip;
    /* Must keep normal and circulation order of triangles consistent.
       Normals will be found by same reflection operation as points.
       Will reverse order of triangle vertices in tri-arrays,
       but multiply normals by an extra minus for other lists
       if the reflection has "odd parity". */
    nneg= numberof(where(flip < 0));
    if(nneg & 1) {
      xyz2= xyz2(,0:1:-1,..);
      norm2= norm2(,0:1:-1,..);
    }
    /* add the object to the display list */
    surf3d, do_alpha, nx, ny, xyz2, norm2, color;
  } else {
    /* add the object to the display list */
    surf3d, do_alpha, nx, ny, xyz, norm, color;
  }
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plcolrsurf3d(xyz, norm, color, flip=)
/* DOCUMENT plcolrsurf3d(xyz, norm, color, flip=)

     Draws a lit 2D mesh in 3D using OpenGL.
     Color varies across the mesh and is specified at mesh vertices.

   SEE ALSO: plsurf3d, plm3d, plpoly3d, plcell3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz, color, and norm to be doubles,
     but should function correctly if they happen to be floats.
  */
  /* are the arguments compatible lengths? */
  fn_name= "plcolrsurf3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 3 && dimsx(2) == 3), "mesh vertices must be 3-by-nx-by-ny in "+fn_name;
  nx= dimsx(3);
  ny= dimsx(4);
  ASSERT, allof(dimsx == dimsof(norm)),"coordinates and normals must have the same dimensions in "+fn_name;
  dimc= dimsof(color);
  if(dimc(2) == 3) {
    do_alpha= 0;
  } else if(dimc(2) == 4) {
    do_alpha= 1;
  } else {
    ASSERT, (0), "color must have 3 or 4 components in "+fn_name;
  }
  ASSERT, ( (dimc(1) == 3) && (dimsx(3) == dimc(3)) && (dimsx(4) == dimc(4)) ),
          "coordinates and colors must have the same dimensions in "+fn_name;
  if(typeof(color) == "char") color= float(color/255.0);
  if(!is_void(flip)) {
    dimf= dimsof(flip);
    ASSERT, (dimf(1) == 1 && dimf(2) == 3),"flip must be a 3 element vector in "+fn_name;
    xyz2= xyz*flip;
    norm2= norm*flip;
    /* Must keep normal and circulation order of triangles consistent.
       Normals will be found by same reflection operation as points.
       Will reverse order of triangle vertices in tri-arrays,
       but multiply normals by an extra minus for other lists
       if the reflection has "odd parity". */
    nneg= numberof(where(flip < 0));
    if(nneg & 1) {
      xyz2= xyz2(,0:1:-1,..);
      norm2= norm2(,0:1:-1,..);
      colr2= color(,0:1:-1,..);
    }
    /* add the object to the display list */
    colrsurf3d, do_alpha, nx, ny, xyz2, norm2, colr2;
  } else {
    /* add the object to the display list */
    colrsurf3d, do_alpha, nx, ny, xyz, norm, color;
  }
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func pllines3d(xyz, color)
/* DOCUMENT pllines3d(xyz, color)

     Draws a Polyline in 3D using OpenGL.

   SEE ALSO: plpoly3d, plcell3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz and color to be doubles,
	 but should function correctly if they happen to be floats.
  */
  /* are the arguments compatible lengths? */
  fn_name= "pllines3d";
  dimsx= dimsof(xyz);
  dimc= dimsof(color);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3), "vertices must be 3-by-nvert in "+fn_name;
  if(typeof(color) == "char") color= float(color/255.0);
  nvert= dimsx(3);
  ASSERT, (dimc(1) == 1 && dimc(2) == 3 ), "color array must be a 3 element vector in "+fn_name;
  /* add the object to the display list */
  lines3d, nvert, xyz, color;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plpoints3d(xyz, color)
/* DOCUMENT plpoints3d(xyz, color)

     Draws a set of Points in 3D using OpenGL.

   SEE ALSO: plpoly3d, plcell3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz and color to be doubles,
     but should function correctly if they happen to be floats.
  */
  /* are the arguments compatible lengths? */
  fn_name= "plpoints3d";
  dimsx= dimsof(xyz);
  dimc= dimsof(color);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3), "vertices must be 3-by-nvert in "+fn_name;
  if(typeof(color) == "char") color= float(color/255.0);
  nvert= dimsx(3);
  ASSERT, (dimc(1) == 2 && dimc(2) == 3 && dimc(3) == nvert ),
          "color array must be 3-by-nvert in "+fn_name;
  /* add the object to the display list */
  points3d, nvert, xyz, color;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plglyphs3d(xyz, scal, theta, phi, color)
/* DOCUMENT plglyphs3d(xyz, scal, theta, phi, color)

     Draws a set of Glyphs using OpenGL. The glyphs are located
     at xyz and have color "color". The base glyph size is multiplied
     scal and the glyph axis points in the direction (theta, phi).

   SEE ALSO: plpoints3d, lookat3d, light3d, win3d
 */
{
  /* This function expects the arguments to be doubles,
     but should function correctly if they happen to be floats.
  */
  /* are the arguments compatible lengths? */
  fn_name= "plglyphs3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3), "points must be 3-by-nglyph in "+fn_name;
  nglyph= dimsx(3);
  dims= dimsof(scal);
  ASSERT, (dims(1) == 1 && dims(2) == nglyph ), "scale vector must be nglyph long in "+fn_name;
  dimt= dimsof(theta);
  ASSERT, (dimt(1) == 1 && dimt(2) == nglyph ), "theta vector must be nglyph long in "+fn_name;
  dimp= dimsof(phi);
  ASSERT, (dimp(1) == 1 && dimp(2) == nglyph ), "phi vector must be nglyph long in "+fn_name;
  dimc= dimsof(color);
  ASSERT, (dimc(1) == 2 && dimc(2) == 3 && dimc(3) == nglyph ),
          "color array must be 3-by-nglyph in "+fn_name;
  if(typeof(color) == "char") color= float(color/255.0);
  /* add the object to the display list */
  glyphs3d, nglyph, xyz, scal, theta, phi, color;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func pltrilists3d(tris,flip=,offset=,cubemap=,emit=)
/* DOCUMENT pltrilists3d(tris,flip=,offset=,cubemap=,emit=)

     This function is used to display an arbitrary number
     of triangle strips or arrays in 3D. An example is the 
     output of an iso-surface routine.
     tris is either a TriStripGrp, a TriArrayGrp, or a TriVertexGrp.
     flip, if present, is a 3 element vector that multiplies
     the x, y, and z coords of every point and normal.
     Ordinarily, all elements are either +1 or -1.
     If cubemap is non-zero, perform specular lighting using
     cube map textures.
     A TriStripGrp is composed of several triangle strips.
     It has a pointer to the next TriStripGrp.
     A TriArrayGrp is an array of triangles with a pointer to the
     next TriArrayGrp. 
     A TriVertexGrp is like a TriArrayGrp except that the
     vertices and normals have been gathered into vectors 
     and each triangle is specified by three indices
     into these vectors.

   SEE ALSO: plcell3d, lookat3d, light3d, win3d
 */
{
  extern n_tri_3d;

  fn_name= "pltrilists3d";
  if(is_void(tris)) return;
  if(!is_void(flip)) {
    dimf= dimsof(flip);
    ASSERT, (dimf(1) == 1 && dimf(2) == 3), "flip must be a 3 element vector in "+fn_name;
    /* Must keep normal and circulation order of triangles consistent.
       Normals will be found by same reflection operation as points.
       Will reverse order of triangle vertices in tri-arrays,
       but multiply normals by an extra minus for other lists
       if the reflection has "odd parity". */
    nneg= numberof(where(flip < 0));
    if(nneg & 1) {
      /* negative number of flips */
      flipn= flip*[-1.0,-1.0,-1.0];
      odd_parity= 1;
    } else {
      flipn= flip;
      odd_parity= 0;
    }
  }
  if(!is_void(offset)) {
    dimo= dimsof(offset);
    ASSERT, (dimo(1) == 1 && dimo(2) == 3), "offset must b e a 3 element vector in "+fn_name;
  }
  if(is_void(cubemap)) {
    cubemap= 0;
  } else if(cubemap) {
    /* CURRENTLY ONLY SUPPORT CUBE MAP LIGHTING FOR TRIANGLE ARRAYS */
	if(structof(tris) != TriArrayGrp) cubemap= 0;
  }
  if(is_void(emit)) {
    emit= 0;
  }
  /* NOTE: Currently create an object in yorick's display list for
     every object in the list of triangle groups.
     Could consider collapsing all of them into one large 
     object.
  */
  /* run over the triangle list */
  nTriTot= nStripTot= nVertTot=0;
  triptr= &tris;
  if(structof(tris) == TriStripGrp) {
    while(1) {
      if(!triptr) break;
      if(is_void(*triptr)) break;
      nStrip= triptr->nStrip;
      nVert= triptr->nVert;
      if(nStrip) {
        ntri= nVert-2*nStrip;
        if(!is_void(flip)) {
          xyz2= (*(triptr->xyzverts))*flip;
          norm2= (*(triptr->normals))*flipn;
          if(!is_void(offset)) {
            xyz2 += offset;
          }
          /* normals for smooth shading must have been supplied */
          pltstrips3d, (*(triptr->triLen)), 
                          xyz2, (*(triptr->colors)), norm2;
        } else {
          /* normals for smooth shading must have been supplied */
          xyz2= (*(triptr->xyzverts));
          if(!is_void(offset)) {
            xyz2 += offset;
          }
          pltstrips3d, (*(triptr->triLen)), xyz2, 
		          (*(triptr->colors)), (*(triptr->normals));
        }
        nTriTot += ntri;
        nStripTot += nStrip;
        nVertTot += nVert;
      }
      triptr= triptr->next;
    } ;
    return [nTriTot, nStripTot, nVertTot];
  } else if(structof(tris) == TriArrayGrp) {
    while(1) {
      if(!triptr) break;
      if(is_void(*triptr)) break;
      ntri= triptr->numTri;
      if(ntri <= 0) break; 
      colr= *(triptr->colors);
      dimc= dimsof(colr);
      ASSERT, (dimc(1) > 0 && (dimc(2) == 3 || dimc(2) == 4) ),
              "a color must have leading length 3 or 4 in "+fn_name;
      if(dimc(1) == 1) {
        /* broadcast the color out to the proper length */
        colr= array(colr, ntri);
      } else if(dimc(1) == 2 && dimc(3) != ntri) {
        ASSERT, (0), "there must be ntri colors in "+fn_name;
      } else if(dimc(1) == 3 && (dimc(3) != 3 || dimc(4) != ntri) ) {
        ASSERT, (0), "there must be 3-by-ntri colors in "+fn_name;
      } else if(dimc(1) > 3) {
        ASSERT, (0), "there must be one color, a color per triangle, or a color per vertex in "+fn_name;
      }
      if(!is_void(flip)) {
        xyz2= (*(triptr->xyzverts))*flip;
        norm2= (*(triptr->normals))*flip;
        if(odd_parity) {
          xyz2= xyz2(,0:1:-1,..);
          norm2= norm2(,0:1:-1,..);
          if(dimc(1) == 3 && dimc(4) == ntri) {
            colr= colr(,0:1:-1,..);
          }
        }
        if(!is_void(offset)) {
          xyz2 += offset;
        }
        pltarray3d, xyz2, norm2, colr, ntri, cubemap, emit;
      } else {
        xyz2= (*(triptr->xyzverts));
        if(!is_void(offset)) {
          xyz2 += offset;
        }
        pltarray3d, xyz2, *(triptr->normals), colr, ntri, cubemap, emit;
      }
      nTriTot += ntri;
      triptr= triptr->next;
    }
    /* logically one triangle per "strip" */
    return [nTriTot, nTriTot, 3*nTriTot];
  } else if(structof(tris) == TriVertexGrp) {
    while(1) {
      if(!triptr) break;
      if(is_void(*triptr)) break;
      ntri= triptr->numTri;
      if(ntri <= 0) break; 
      nvert= triptr->numEdg
      colr= *(triptr->colors);
      dimc= dimsof(colr);
      ASSERT, (dimc(1) > 0 && (dimc(2) == 3 || dimc(2) == 4) ),
              "a color must have leading length 3 or 4 in "+fn_name;
      if(dimc(2) == 3) {
        // convert color from RGB to RGBA
        dimc(2)= 4;
        colrnew= array(1.0f, dimc);
        colrnew(1:3,..)= colr;
        colr= colrnew;    colrnew= [];
      }
      if(dimc(1) == 1) {
        /* broadcast the color out to the proper length */
        colr= array(colr, nvert);
      } else if(dimc(1) == 2) {
        ASSERT, (dimc(3) == nvert), "there must be nvert colors in "+fn_name;
      } else {
        write,"color has dimensions", dimc(1), dimc(2), dimc(3), dimc(4);
        ASSERT, (0), "there should be nvert colors in "+fn_name;
      }
      if(!is_void(flip)) {
        xyz2= (*(triptr->xyzverts))*flip;
        norm2= (*(triptr->normals))*flipn;
        if(!is_void(offset)) {
          xyz2 += offset;
        }
        if(use_interleave) pltivarray3d, *(triptr->ptndx), xyz2, norm2, colr, ntri, nvert;
        else pltvarray3d, *(triptr->ptndx), xyz2, norm2, colr, ntri, nvert;
      } else {
        xyz2= (*(triptr->xyzverts));
        if(!is_void(offset)) {
          xyz2= offset+(*(triptr->xyzverts));
          if(use_interleave) pltivarray3d, *(triptr->ptndx), xyz2, *(triptr->normals), colr, ntri, nvert;
          else pltvarray3d, *(triptr->ptndx), xyz2, *(triptr->normals), colr, ntri, nvert;
        } else {
          if(use_interleave) pltivarray3d, *(triptr->ptndx), *(triptr->xyzverts), *(triptr->normals), colr, ntri, nvert;
          else pltvarray3d, *(triptr->ptndx), *(triptr->xyzverts), *(triptr->normals), colr, ntri, nvert;
        }
      }
      nTriTot += ntri;
      triptr= triptr->next;
    }
    /* logically one triangle per "strip" */
    return [nTriTot, nTriTot, 3*nTriTot];
  } else if(structof(tris) == TriStripNdxGrp) {
    while(1) {
      if(!triptr) break;
      if(is_void(*triptr)) break;
      nVert= triptr->nVert;
      nStrip= triptr->nStrip;
      ntri= nVert-2*nStrip;
      numedg= triptr->numEdg;
      if(nStrip <= 0) break; 
      colr= *(triptr->colors);
      dimc= dimsof(colr);
      ASSERT, (dimc(1) == 1 && (dimc(2) == 3 || dimc(2) == 4) ),
              "color must be a 3 or 4 element vector in "+fn_name;
      if(!is_void(flip)) {
        xyz2= (*(triptr->xyzverts))*flip;
        norm2= (*(triptr->normals))*flipn;
        if(!is_void(offset)) {
          xyz2 += offset;
        }
        pltivstrips3d, *(triptr->triLen), *(triptr->ptndx), xyz2, 
		                  norm2, colr;
      } else {
        xyz2= (*(triptr->xyzverts));
        if(!is_void(offset)) {
          xyz2 += offset;
        }
        pltivstrips3d, *(triptr->triLen), *(triptr->ptndx), xyz2, 
		                  *(triptr->normals), colr;
      }
      nTriTot += ntri;
      nStripTot += nStrip;
      nVertTot += nVert;
      triptr= triptr->next;
    }
    return [nTriTot, nStripTot, nVertTot];
  } else {
    ASSERT, (0),"the tris argument is neither a TriArrayGrp, a TriVertexGrp, TriStripNdxGrp, nor a TriStripGrp in "+fn_name;
  }
}

func pltstrips3d(nv, xyz, color, norm, draw_edge=)
/* DOCUMENT pltstrips3d(nv, xyz, color, norm, draw_edge=)

     Draws a set of triangle-strips in 3D using OpenGL.
     pltrilists3d is normally called instead of this 
     lower level routine.
     xyz is 3-by-nvert.
     nv is the number of vertices per triangle strip (i.e. the length 
     of the run in the last index of xyz).
     nvert= sum(nv)
     nstrip= numberof(nv);
     ntri= nvert-2*nstrip;
     color is normally 3-by-ncolor. ncolor can be 1 (all tris the 
     same color), numberof(nv)  (one color per strip), or
     ntri (one color per triangle). color will always be increased in
     size to ntri before being stored.
     If color is 4-by-ncolor, the 4th element is the alpha value. 
     norm can be 3-by-ntri for flat shading or 3-by-nvert
     for smooth shading. If norm is not specified, lighting
     will not be used (the specified color will be used directly).


   SEE ALSO: pltrilists3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz, color, and norm to be doubles,
     but should function correctly if they happen to be floats.
  */
  /* This function accepts a wide range of inputs and coerces
     them into a more restricted set that are saved in yorick's
     display list. */
  fn_name= "pltstrips3d";
  nvert= sum(nv);
  nstrip= numberof(nv);
  ntri= nvert-2*nstrip;
  if(nstrip <= 0) return;
  /* The float conversions below makes sure that the right data
     type is passed to the compiled routine.
  */
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3 && dimsx(3) == nvert),
          "xyz must be a 3-by-nvert array in "+fn_name;
  dimc= dimsof(color);
  ASSERT, (dimc(2) == 3 || dimc(2) == 4), "first dimension of color must be 3 or 4 in "+fn_name;
  if(dimc(2) == 3) {
    do_alpha= 0;
  } else {
    do_alpha= 1;
  }
  if(typeof(color) == "char") color= float(color/255.0);
  if(dimc(1) == 1) {
    /* broadcast color out to all tris */
    color= array(color, ntri);
  } else {
    ASSERT, (dimc(1) == 2), "color can have at most 2 dimensions in "+fn_name;
    if(dimc(3) == 1) {
      /* broadcast color out to all tris */
	  color= array(color(,1), ntri);
    } else if(dimc(3) == nstrip) {
      clrnew= array(0.0, dimc(2), ntri);
      for(i= 1, base= 1; i <= nstrip; i++) {
        ntr= nv(i)-2;
        clrnew(,base:base+ntr-1)= color(,i);
        base += ntr;
      }
      /* NOTE: eq_nocopy doesn't work unless arg. 2 is a temporary
         (there is no copy involved here, just that there is a 
         temporary variable on the yorick stack).
      */
      eq_nocopy, color, *&clrnew;
      clrnew= [];
    } else {
      ASSERT, (dimc(3) == ntri), "number of colors must be 1, nstrip, or ntri in "+fn_name;
    }
  }
  /* the user can request outlined edges, not filled polygons */
  do_edge= 0;
  if(!is_void(draw_edge) && draw_edge) do_edge= 1;
  /* Do not light the surface unless normals were supplied. 
     If there are as many normals as vertices, use smooth
     shading */
  smooth= 0;
  do_light= 1;
  if(is_void(norm)) {
    do_light= 0;
    norm= 0;
  } else {
    dimsn= dimsof(norm);
    ASSERT, (dimsn(1) == 2 && dimsn(2) == 3), "norm must be a 3-by-n array in "+fn_name;
    if(dimsn(3) == dimsx(3)) {
      smooth= 1; /* normal per vertex, so use smooth shading */
    } else if(dimsn(3) == ntri) {
      smooth= 0; /* normal per triangle, so use flat shading */
    } else {
      ASSERT, (0),"norm must have length nvertices or ntris in "+fn_name;
    }
  }
  /* add the object to the display list */
  tstrips3d, nstrip, nv, xyz, norm, color, do_edge, smooth, do_light, do_alpha;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plqstrips3d(nv, xyz, color, norm, draw_edge=)
/* DOCUMENT plqstrips3d(nv, xyz, color, norm, draw_edge=)

     Draws a set of quad-strips in 3D using OpenGL.
     xyz is 3-by-2-by-nvert.
     nv is the number of vertices per quad strip (i.e. the length 
     of the run in the last index of xyz).
     nvert= sum(nv)
     nstrip= numberof(nv);
     nquad= nvert-nstrip;
     color is 3-by-ncolor. ncolor can be 1 (all quads the same color),
     numberof(nv)  (one color per strip), or
     nquad (one color per quad). color will always be increased in
     size to nquad before being stored.
     norm can be 3-by-nquad for flat shading or 3-by-2-by-nvert
     for smooth shading. If norm is not specified, lighting
     will not be used (the specified color will be used directly).

   SEE ALSO: pltrilists3d, pltstrips3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz, color, and norm to be doubles,
     but should function correctly if they happen to be floats.
  */
  /* This function accepts a wide range of inputs and coerces
     them into a more restricted set that are saved in yorick's
     display list. */
  fn_name= "plqstrips3d";
  nvert= sum(nv);
  nstrip= numberof(nv);
  nquad= nvert-nstrip;
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 3 && dimsx(2) == 3 && dimsx(3) == 2 && dimsx(4) == nvert),
          "xyz must be a 3-by-2-by-nvert array in "+fn_name;
  dimc= dimsof(color);
  ASSERT, (dimc(2) == 3 || dimc(2) == 4), "first dimension of color must be 3 or 4 in "+fn_name;
  if(dimc(2) == 3) {
    do_alpha= 0;
  } else {
    do_alpha= 1;
  }
  if(typeof(color) == "char") color= float(color/255.0);
  if(dimc(1) == 1) {
    /* broadcast color out to all quads */
    color= array(color, nquad);
  } else {
    ASSERT, (dimc(1) == 2), "color can have at most 2 dimensions in "+fn_name;
    if(dimc(3) == 1) {
      /* broadcast color out to all tris */
      color= array(color(,1), nquad);
    } else if(dimc(3) == nstrip) {
      clrnew= array(0.0, dimc(2), nquad);
      for(i= 1, base= 1; i <= nstrip; i++) {
        nq= nv(i)-1;
        clrnew(,base:base+nq-1)= color(,i);
        base += nq;
      }
      /* NOTE: eq_nocopy doesn't work unless arg. 2 is a temporary
         (there is no copy involved here, just that there is a 
    	 temporary variable on the yorick stack).
      */
      eq_nocopy, color, *&clrnew;
      clrnew= [];
    } else if(dimc(3) != nquad) {
      ASSERT, (0),"number of colors must be 1, nstrip, or ntri in "+fn_name;
    }
  }
  /* the user can request outlined edges, not filled polygons */
  do_edge= 0;
  if(!is_void(draw_edge) && draw_edge) do_edge= 1;
  /* Do not light the surface unless normals were supplied. 
     If there are as many normals as vertices, use smooth
     shading */
  do_light= 1;
  if(is_void(norm)) {
    do_light= 0;
    norm= 0;
    smooth= 0;
  } else {
    dimsn= dimsof(norm);
    if(dimsn(1) == 2) {
      /* normals must be one per quad */
      ASSERT, (dimsn(3) == nquad), "norm must be a 3-by-nquad array in "+fn_name;
      smooth= 0; /* normal per quad, so use flat shading */
    } else if(dimsn(1) == 3) {
      /* normals must be one per vertex */
      ASSERT, (dimsn(2) == 3 && dimsn(3) == 2 && dimsn(4) == dimsx(4)),
              "norm must a 3-by-2-by-nvert array in "+fn_name;
      smooth= 1; /* normal per vertex, so use smooth shading */
    } else {
      ASSERT, (0),"norm must be a 3-by-nquad or a 3-by-2-by-nvert array in "+fn_name;
    }
  }
  /* add the object to the display list */
  qstrips3d, nstrip, nv, xyz, norm, color, do_edge, smooth, do_light, do_alpha;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func pltivstrips3d(nv, ndx, xyz, norm, color, draw_edge=)
/* DOCUMENT pltivstrips3d(nv, ndx, xyz, norm, color, draw_edge=)

     Draws a set of triangle-strips in 3D using OpenGL.
     pltrilists3d is normally called instead of this 
     lower level routine.
     xyz and norm are 3-by-nvert.
     ndx holds the vertex numbers (indices into xyz or norm)
     of the vertices in the tristrips.
     nv is the number of vertices per triangle strip (i.e. the length 
     of the run in the last index of ndx).
     nvert= sum(nv)
     nstrip= numberof(nv);
     ntri= nvert-2*nstrip;
     color is normally 3-by-ncolor. ncolor can be 1 (all tris the 
     same color), numberof(nv)  (one color per strip), or
     ntri (one color per triangle). color will always be increased in
     size to ntri before being stored.
     If color is 4-by-ncolor, the 4th element is the alpha value. 

   SEE ALSO: pltrilists3d, pltristrips3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz, color, and norm to be doubles,
	 but should function correctly if they happen to be floats.
  */
  /* This function accepts a wide range of inputs and coerces
     them into a more restricted set that are saved in yorick's
     display list. */
  fn_name= "pltivstrips3d";
  nvert= sum(nv);
  nstrip= numberof(nv);
  ntri= nvert-2*nstrip;
  numedg= numberof(xyz(1,..));
  if(nstrip <= 0) return;
  /* The float conversions below makes sure that the right data
     type is passed to the compiled routine.
  */
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3 && dimsx(3) == numedg),
          "xyz must be a 3-by-numedg array in "+fn_name;
  if(is_void(norm)) {
      ASSERT, (0), "norm must be supplied in "+fn_name;
  } else {
    dimsn= dimsof(norm);
    ASSERT, (dimsn(1) == 2 && dimsn(2) == 3 && dimsx(3) == numedg),
            "norm must be a 3-by-numedg array in "+fn_name;
  }
  dimc= dimsof(color);
  ASSERT, (dimc(2) == 3 || dimc(2) == 4), "first dimension of color must be 3 or 4 in "+fn_name;
  if(dimc(2) == 3) {
    do_alpha= 0;
  } else {
    do_alpha= 1;
  }
  if(typeof(color) == "char") color= float(color/255.0);
  if(dimc(1) == 1) {
    /* broadcast color out to all tris */
    color= array(color, ntri);
  } else {
    ASSERT, (dimc(1) == 2), "color can have at most 2 dimensions in "+fn_name;
    if(dimc(3) == 1) {
      /* broadcast color out to all tris */
	  color= array(color(,1), ntri);
    } else if(dimc(3) == nstrip) {
      clrnew= array(0.0, dimc(2), ntri);
      for(i= 1, base= 1; i <= nstrip; i++) {
        ntr= nv(i)-2;
        clrnew(,base:base+ntr-1)= color(,i);
        base += ntr;
      }
      /* NOTE: eq_nocopy doesn't work unless arg. 2 is a temporary
         (there is no copy involved here, just that there is a 
         temporary variable on the yorick stack).
      */
      eq_nocopy, color, *&clrnew;
      clrnew= [];
    } else if(dimc(3) != ntri) {
      ASSERT, (0),"number of colors must be 1, nstrip, or ntri in "+fn_name;
    }
  }
  /* the user can request outlined edges, not filled polygons */
  do_edge= 0;
  if(!is_void(draw_edge) && draw_edge) do_edge= 1;
  /* add the object to the display list */
  tstripsndx3d, nstrip, numedg, ntri, nv, ndx, xyz, norm, color, do_edge, do_alpha;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func pltarray3d(xyz, norm, color, ntri, cubemap, emit, draw_edge=)
/* DOCUMENT pltarray3d(xyz, norm, color, ntri, cubemap, emit, draw_edge=)

     Draws a Triangle Array in 3D using OpenGL.
     pltrilists3d is normally called instead of this
     lower level routine.
     If cubemap is non-zero, specular lighting is handled by 
     cube map textures.

   SEE ALSO: pltrilists3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz, color, and norm to be doubles,
     but should function correctly if they happen to be floats.
  */
  /* make private copies of all input arrays so that they can 
     be used as OpenGL vertex arrays (whose storage must
     remain undisturbed until the picture is cleared)
  */
  fn_name= "pltarray3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 3 && dimsx(2) == 3 && dimsx(3) == 3 && dimsx(4) >= ntri),
          "triangle array vertices must be 3-by-3-by-ntri in "+fn_name;
  dimc= dimsof(color);
  ASSERT, (dimc(2) == 3 || dimc(2) == 4), "first dimension of color must be 3 or 4 in "+fn_name;
  if(dimc(2) == 3) {
    do_alpha= 0;
  } else {
    do_alpha= 1;
  }
  ASSERT, (dimc(1) <= 3), "color can have at most 3 dimensions in "+fn_name;
  if(typeof(color) == "char") color= float(color/255.0);
  if(dimc(1) == 1) {
    /* broadcast color out to all tris */
    color= array(color, ntri);
    colrpervrt= 0;
  } else {
    if(dimc(1) == 3) {
      ASSERT, (dimc(3) == 3 && dimc(4) == ntri), "number of colors must be 1 or ntri or 3*ntri in "+fn_name;
      cpervrt= 1;  /* triangles have a color per vertex */
    } else {
      if(dimc(3) == 1) {
        /* broadcast color out to all tris */
        color= array(color(,1), ntri);
      } else if(dimc(3) != ntri) {
        ASSERT, (0),"number of colors must be 1 or ntri in "+fn_name;
      }
      cpervrt= 0;
    }
  }
  if(is_void(norm)) {
    smooth= 1;
    norm= array(0.0, 3, 4, ntri);
    ln= array(3, ntri);
    nr= get_normal3d(xyz, ln);
    norm(,1,)= nr;
    norm(,2,)= nr;
    norm(,3,)= nr;
    nr= [];
  } else {
    dimsn= dimsof(norm);
    if(dimsn(1) == 2) {
      smooth= 0;
      ASSERT, (dimsn(2) == 3 && dimsn(3) >= ntri), "triangle array normals must be 3-by-ntri in "+fn_name;
      /* replicate to get one normal per vertex to make the
         array compatible with OpenGL 1.1 vertex arrays
      */
      norm= norm(,-:1:3,);
    } else if(dimsn(1) == 3) {
      smooth= 1;
      ASSERT, (dimsn(2) == 3 && dimsn(3) == 3 && dimsn(4) == ntri),
              "triangle array normals must be 3-by-3-by-ntri in "+fn_name;
    } else {
      ASSERT, (0),"triangle array normals must be 3-by-ntri or 3-by-3-by-ntri in "+fn_name;
    }
  }
  /* the user can request outlined edges, not filled polygons */
  do_edge= 0;
  if(!is_void(draw_edge) && draw_edge) do_edge= 1;
  do_light= 1;
  /* add the object to the display list */
  tarray3d, ntri, xyz, norm, color, do_edge, smooth, do_light, do_alpha, 
            cpervrt, cubemap, emit;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func plqarray3d(xyz, norm, color, nquad, draw_edge=)
/* DOCUMENT plqarray3d(xyz, norm, color, nquad, draw_edge=)

     Draws a Quadrangle Array in 3D using OpenGL.

   SEE ALSO: pltarray3d, plqstrips3d, lookat3d, light3d, win3d
 */
{
  /* This function expects xyz, color, and norm to be doubles,
	 but should function correctly if they happen to be floats.
  */
  /* make private copies of all input arrays so that they can 
     be used as OpenGL vertex arrays (whose storage must
     remain undisturbed until the picture is cleared)
  */
  fn_name= "plqarray3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 3 && dimsx(2) == 3 && dimsx(3) == 4 && dimsx(4) == nquad),
          "quad array vertices must be 3-by-4-by-nquad in "+fn_name;
  dimc= dimsof(color);
  ASSERT, (dimc(2) == 3 || dimc(2) == 4), "first dimension of color must be 3 or 4 in "+fn_name;
  if(dimc(2) == 3) {
    do_alpha= 0;
  } else {
    do_alpha= 1;
  }
  ASSERT, (dimc(1) > 3), "color can have at most 2 dimensions in "+fn_name;
  if(typeof(color) == "char") color= float(color/255.0);
  if(dimc(1) == 1) {
    /* broadcast color out to all tris */
    color= array(color, ntri);
    colrpervrt= 0;
  } else {
    if(dimc(1) == 3) {
      ASSERT, (dimc(3) == 4 && dimc(4) == nquad),
              "number of colors must be 1 or nquad or 4*nquad in "+fn_name;
      cpervrt= 1;  /* triangles have a color per vertex */
    } else {
      if(dimc(3) == 1) {
        /* broadcast color out to all quads */
        color= array(color(,1), nquad);
      } else if(dimc(3) != nqaud) {
        ASSERT, (0), "number of colors must be 1 or nquad in "+fn_name;
      }
      cpervrt= 0;
    }
  }
  /* normals could be either one per quad or one per vertex.
     one per vertex means smooth shading. */
  if(is_void(norm)) {
    smooth= 1;
    norm= array(0.0, 3, 4, nquad);
	ln= array(4, nquad);
    nr= get_normal3d(xyz, ln);
    norm(,1,)= nr;
    norm(,2,)= nr;
    norm(,3,)= nr;
    norm(,4,)= nr;
    nr= [];
  } else {
    dimsn= dimsof(norm);
    if(dimsn(1) == 2) {
      smooth= 0;
      ASSERT, (dimsn(2) == 4 && dimsn(3) == nquad),
              "quad array normals must be 4-by-nquad in "+fn_name;
      /* replicate to get one normal per vertex to make the
         array compatible with OpenGL 1.1 vertex arrays
      */
      norm= norm(,-:1:4,);
    } else if(dimsn(1) == 3) {
      smooth= 1;
      ASSERT, (dimsn(2) == 3 && dimsn(3) == 4 && dimsn(4) >= nquad),
              "triangle array normals must be 3-by-4-by-nquad in "+fn_name;
    } else {
      ASSERT, (0),"triangle array normals must be 3-by-nquad or 3-by-4-by-nquad in "+fn_name;
    }
  }
  /* the user can request outlined edges, not filled polygons */
  do_edge= 0;
  if(!is_void(draw_edge) && draw_edge) do_edge= 1;
  do_light= 1;
  /* add the object to the display list */
  qarray3d, nquad, xyz, norm, color, do_edge, smooth, do_light, do_alpha, cpervrt;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func pltivarray3d(ptndx, xyz, norm, color, ntri, nvert)
/* DOCUMENT pltivarray3d

     Draws a Triangle Array in 3D using OpenGL.
     Triangles are specified using indices into a list of vertices,
     normals, and colors.
     pltrilists3d is normally called instead of this
     lower level routine.

   SEE ALSO: pltarray3d, lookat3d, pltrilists3d
 */
{
  /* make a private interleaved array containing all input arrays 
     so that it can be used as OpenGL vertex array
  */
  fn_name= "pltivarray3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3 && dimsx(3) >= nvert),
          "triangle array vertices must be 3-by-n_edges_cut in "+fn_name;
  if(is_void(norm)) {
    ASSERT, (0), "normals must be supplied in "+fn_name;
  } else {
    dimsn= dimsof(norm);
    ASSERT, allof(dimsn == dimsx), "normals and xyz must be the same size in "+fn_name;
  }
  /* there must be an RGBA color for every vertex */
  dimc= dimsof(color);
  ASSERT, (dimc(1) == 2 && dimc(2) == 4 && dimc(3) >= nvert),
          "triangle array colors must be 4-by-n_edges_cut in "+fn_name;
  if(typeof(color) == "char") color= float(color/255.0);
  tivarray3d, ntri, nvert, ptndx, xyz, norm, colr;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func pltvarray3d(ptndx, xyz, norm, color, ntri, nvert)
/* DOCUMENT pltvarray3d

     Draws a Triangle Array in 3D using OpenGL.
     Triangles are specified using indices into a list of vertices,
     normals, and colors.
     pltrilists3d is normally called instead of this
     lower level routine.

   SEE ALSO: pltarray3d, pltivarray3d, lookat3d, pltrilists3d
 */
{
  fn_name= "pltvarray3d";
  dimsx= dimsof(xyz);
  ASSERT, (dimsx(1) == 2 && dimsx(2) == 3 && dimsx(3) >= nvert),
          "triangle array vertices must be 3-by-n_edges_cut in "+fn_name;
  if(is_void(norm)) {
    ASSERT, (0), "normals must be supplied in "+fn_name;
  } else {
    dimsn= dimsof(norm);
    ASSERT, allof(dimsn == dimsx), "normals and xyz must be the same size in "+fn_name;
  }
  /* there must be a RGB color for every vertex */
  dimc= dimsof(color);
  ASSERT, (dimc(2) == 3 || dimc(2) == 4),
          "triangle array colors must be RGB or RGBA in "+fn_name;
  if(dimc(2) == 4) do_alpha= 1;
  else do_alpha= 0;
  if(dimc(1) == 1) {
    cpervrt= 0;
  } else {
    cpervrt= 1;
    ASSERT, (dimc(1) == 2 && dimc(3) >= nvert),
          "must have n_edges_cut colors in "+fn_name;
  }
  if(typeof(color) == "char") color= float(color/255.0);
  tvarray3d, ntri, nvert, do_alpha, cpervrt, ptndx, xyz, norm, colr;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func getpix3d(nx,ny)
/* DOCUMENT getpix3d

     Read pixels from an OpenGL window.
*/
{
  /* if the input sizes are bigger than the window, reduce 
     them before reading the pixels */
  maxnx= get_width3d();
  if(is_void(nx)) {
    nx= maxnx
  } else {
    nx= min(nx, maxnx);
  }
  maxny= get_hite3d();
  if(is_void(ny)) {
    ny= maxny;
  } else {
    ny= min(ny, maxny);
  }
  pix= array(char, 3, nx, ny);
  grabpix3d,nx,ny,pix;
  return pix;
}

func putpix3d(arr)
/* DOCUMENT putpix3d

     Draws pixels into an OpenGL window.

   SEE ALSO: clear3
 */
{
  dima= dimsof(arr);
  nx= dima(3);
  ny= dima(4);
  /* add the object to the display list */
  plpix3d, nx, ny, arr;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func CollapseTri(tris)
/* DOCUMENT CollapseTri(tris)

     This function is used to collapse an arbitrary number
     of triangle arrays into a single triangle array.
     tris must be a TriArrayGrp (which can point to another
     TriArrayGrp, etc.).

   SEE ALSO: pltrilists3d, sortTri
 */
{
  local newtri, triptr, nTriTot;

  if(is_void(tris)) return [];
  /* run over the triangle list */
  ASSERT, (structof(tris) == TriArrayGrp), "CollapseTri requires a TriArrayGrp as its argument";
  /* count the number of triangles in all the lists
     and at the same determine whether colors are
     RGB or RGBA and whether there is a color per triangle
     or a color per array. */
  dimc= dimsof(*(tris.colors));
  ASSERT, (dimc(1) >= 1), "CollapseTri invalid color in a TriArrayGrp";
  ncomp= dimc(2);
  nTriTot= 0;
  triptr= &tris;
  while(triptr && !is_void(*triptr)) {
    nTriTot += triptr->numTri;
    triptr= triptr->next;
  }
  if(nTriTot <= 0) return [];
  /* colrtyp is the number of color components per triangle.
     It is negative if there is a single color for all
     triangles. The program is likely to CRASH if the number
     of color components varies in the list or if some
     parts have a color per triangle and some have
     a color per list. */
  colrtyp= ncomp; 
  if(dimc(1) <= 1 || dimc(3) <= 1) colrtyp= -colrtyp;
  color_per_vert= 0;
  if(dimc(1) == 3) {
    // apparently one color per triangle vertex
    if(dimc(3) != 3 || numberof((*tris.colors)(1,1,)) != numberof((*tris.colors)(1,1,)) ) {
      error,"Dimensions of colors are wrong for color per vertex in CollapseTris";
    } else {
      // set color per vertex flag
      color_per_vert= 1;
      if(colrtyp < 0) colrtyp -= 16;
      else colrtyp += 16;
    }
  }
  /* xyzverts and normals are really array(Point3D, 3, nTriTot) */
  newxyz= array(double, 3, 3, nTriTot );
  newnorm= array(double, 3, 3, nTriTot );
  if(color_per_vert) {
    newcolr= array(float, ncomp, 3, nTriTot );
  } else {
    newcolr= array(float, ncomp, nTriTot );
  }
  newids= array(long, nTriTot );
  triptr= &tris;
  if( triptr->var2 && !is_void(*(triptr->var2)) ) {
    nvar2= array(double, 3, nTriTot);
  } else {
    nvar2= nulvar;
  }
  nulvar= [];
  newtri= TriArrayGrp(numTri= nTriTot,
        xyzverts= &newxyz,
        normals= &newnorm,
        colors= &newcolr,
        cellIDs= &newids,
        var2= &nvar2,
        /* the next 3 are not used */
        nTris= &nulvar, triEdg= &nulvar, triStart= &nulvar,
        next= &nulvar);
  triptr= &tris;
  CollapseTriArrays3d, colrtyp, triptr, &newtri;
  return newtri;
}

func SortTri(tris)
/* DOCUMENT SortTri(tris)

     This function performs a depth sort of a triangle
     array. It uses the current viewpoint.
     The resulting triangle array should display correctly
     with translucency when viewed from the same point.
     Rotating the scene may lead to incorrect results
     if the triangles are translucent.
     tris must be a TriArrayGrp.

   SEE ALSO: pltrilists3d, collapseTri
 */
{
  local newxyz, newnorm, newcolr, newids;

  ntri= tris.numTri;
  dimc= dimsof(*(tris.colors));
  ncomp= dimc(2);
  newxyz= array(double, 3, 3, ntri);
  newnorm= array(double, 3, 3, ntri);
  newcolr= array(float, ncomp, ntri);
  newids= array(long, ntri);
  hasVar2= ( tris.var2 && !is_void(*(tris.var2)) );
  if(hasVar2) {
    nvar2= array(double, 3, nTriTot);
  } else {
    nvar2= nulvar;
  }
  nulvar= [];
  newtri= TriArrayGrp(numTri= ntri,
       xyzverts= &newxyz,
        normals= &newnorm,
        colors= &newcolr,
        cellIDs= &newids,
        /* the next 3 are not used except while forming strips */
        nTris= &nulvar,  
        triEdg= &nulvar,
        triStart= &nulvar,
        next= &nulvar);
  DoSortTri3d, ncomp, &tris, &newtri;
  return newtri;
}

func pltex2dvol(delta, texval)
/* DOCUMENT pltex2dvol

     Draws a volume visualization of a 3D cell array using OpenGL.
     Uses 2D textures.

   SEE ALSO: pltex3dvol, lookat3d, light3d, win3d, pltex3dvol
 */
{
  fn_name= "pltex2dvol";
  texval= char(texval);
  dimv= dimsof(texval);
  ASSERT, ((dimv(1) == 4) && (dimv(2) == 4)), "data array must be 4-nx-ny-nz in "+fn_name;
  nx= dimv(3);
  ny= dimv(4);
  nz= dimv(5);
  /* add the object to the display list */
  texcell2d, nx, ny, nz, delta, texval;
  /* increment counter so that later code knows to make a new OpenGL
     display list. */
  inc_seq3d;
  draw3d_trigger;
}

func pltex3dvol(nslab, boxsiz, texval, origin=)
/* DOCUMENT pltex3dvol

     Draws a volume visualization of a 3D cell array using OpenGL.
     Uses 3D textures.

   SEE ALSO: pltex2dvol, lookat3d, light3d, win3d
 */
{
  fn_name= "pltex3dvol";
  if(typeof(texval) != "char") texval= char(texval);
  if(is_void(boxsiz)) boxsiz= [1.0,1.0,1.0];
  dimv= dimsof(texval);
  if(dimv(1) == 4 && dimv(2) == 4) {
    nx= dimv(3);
    ny= dimv(4);
    nz= dimv(5);
  } else {
    ASSERT, (0),"data array must be 4-nx-ny-nz in "+fn_name;
  }
  if(is_void(origin)) origin= [0.0, 0.0, 0.0];
  ds= max(boxsiz)/nslab;
  /* NOTE: The 3D texture will be loaded now, which is good if the 
     object is rotated, but bad if it is not shown even once. */
  ldtex3d, nx, ny, nz, texval;
  /* add the object to the display list */
  tex3d, ds, origin, boxsiz;
  /* Increment counter so that later code knows to make a new OpenGL
     display list. 
     QUESTION - is this needed for "direct" objects? */
//  inc_seq3d;
  draw3d_trigger;
}
