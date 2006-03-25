/*
 * $Id: ContourTets3D.c,v 1.2 2006-03-25 03:12:29 dhmunro Exp $
 * Iso-surface routines based on a 6 tetrahedron decomposition of
 * each hexahedral cell.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

#include "Contour3D.h"
#include "pstdlib.h"
#include <stdlib.h>
#include <math.h>

#define TET_USE_POLYS
#define DBG_ISO

/* 
   These iso-surface routines use a 6 tetrahedron decomposition of each 
   hexahedral zone. The diagonal from the lowest numbered vertex 
   in the zone to the diagonally opposite vertex is used to 
   make the split. The remaining six vertices form six pairs that share an 
   edge of the hexahedron. The six tethedrons are formed by combining those 
   pairs with the two points on the diagonal. The two points on the main diagonal 
   are vertices of every tetrahedron.

   There are 12 edges on the hexahedron. These will be referred to as "real" edges.
   They are numbered zero through 11. There are 7 other edges that occur in 
   the tetrahedrons. Edge 12 is the long diagonal of the hexahedron and edges 
   13-18 are diagonals of the faces of the hexahedron.


   User data may be cell-centered or point-centered. It may have or 
   not have ghost cells. Handling all cases leads to an explosion of
   routines. The underlying iso-surface routines work with 
   point-centered values. Converting a full input array from 
   zone-centered to point-centered can chew up a lot of memory.

   CONCLUSION: The iso-surface routine will work on a chunk in a 
   (potentially) larger array. The iso-surface will be computed using
   point-centered data stored in a "chunk" sized array. The calling 
   routine will provide a function pointer to map data from the original 
   grid to the chunk, point-centering if needed. Gradients will be
   computed on the chunk. The chunk will always have a layer of
   ghost points. The "re-mapping" function will be able to handle ghost 
   cells if the request "over-reaches" the input array.

   Octree based iso-surfacing routines build a tree based on the data 
   range for a cell, then a 2-by-2-by-2 group of cells and so on 
   until they have the data range for a whole chunk. Computing the 
   octree requires point-centered variable values for all vertices 
   in the chunk. Generating an iso-surface requires point-centered 
   variable values and gradients at all corners of a zone cut by 
   the iso-surface, but doesn't need them for other zones. Gradients 
   (in particular) are expensive to compute, so it is faster to only 
   compute the ones that are needed.

   CONCLUSION: The iso-surface routine will take as input storage for 
   a point-centered variable (and storage for an auxiliary variable if
   requested) with one ghost point on all sides of the chunk. There is also 
   chunk-sized input storage for co-ordinates and gradients. There is 
   a chunk-sized input character array that contains bit flags for 
   whether co-ordinates, variables, auxiliary variables, and gradients 
   have been computed at each point. All these variables are computed 
   as needed and tracked with bit flags. If creating an octree, only 
   variable values will be computed. When generating an iso-surface, 
   all of them may be computed.

   The lower level iso-surface routine is now the same for any ijk-indexible
   grid. Differences between point-centered and zone-centered, 
   the presence or absence of ghost zones, and Cartesian versus
   curvi-linear grids are all handled by the "data gathering" 
   routines supplied as arguments.

   There are multiple low-level iso-surface routines, but that is to 
   provide different kinds of output triangles. One approach is 
   an array of triangles, another is a list of triangle strips, 
   another is an array of triangles described by indices into a 
   vertex array, and a fourth is triangle strips described by 
   index lists.

   There is a set of functions for each combination of grid type,
   variable centering, and for the case that a cell-centered 
   variable is the same size as the grid (i.e. the variable 
   size is one greater in each dimension than the size needed 
   to store the variable to allow the same 1D indices to be used for 
   the grid and the variable). An initialization function stores 
   the information needed to map between a point in a chunk 
   and the full array. The second function fetches a function 
   value at a point in the chunk (may involve centering), the
   third fetches the coordinates of a point in a chunk, and the 
   fourth computes the gradient at a point in a chunk.

    ---------------------------------------------------

   offsets[0:2] contains 0-based starting indices of this 
   chunk in the larger array. Yorick must convert to C-style
   indices before calling the compiled routines. These offsets 
   point to the first real cell of this chunk in the full array. 
   We CAN'T use the first ghost cell because the larger array
   may not have ghost cells. If the input array is point-centered, 
   the offset is to the first vertex of the first real cell.

   offsets[3:5] contains the number of CELLS in the larger array 
   (including ghost cells if present) in the "x", "y", and "z" directions. 
   This is used to compute STRIDES in the full array. Strides 
   for coordinate arrays (in the case of curvi-linear coordinates) 
   are obtained by adding one to these values. This choice was made 
   because there is always an input variable array, but there isn't 
   a coordinate array for Cartesian grids.

   sizes[0:2] contains the number of POINTS in each direction in 
   this chunk of the array including the extra point on each side 
   of this chunk. The extra point is a "ghost" from the viewpoint 
   of the chunk, but it may be a real point/cell in the larger array.
   The chunk always includes ghost points, even if they "overlay"
   points outside the full array. A chunk array is NEVER cell-centered.

   deltas[0:2] contains the physical size of each cell. 

   origin[0:2] is the co-ordinate of the first real point. This point
   is interior to the chunk volume and is located between the first
   and second cells of the input chunk.

   The var and var2 (optional) chunk arrays always have point-centered 
   values. The grd array has point-centered gradients. The chunk 
   always has room for ghost points on all sides. 

   The values in var and var2 are copied from the full arrays if they
   are point-centered. They are computed from the surrounding 
   cells in the full array if it iarrays zone-centered. The only time 
   extrapolation is needed is at the boundaries of the full array 
   when the full array doesn't have ghost cells.

   sizes(1:3) is the number of points in the chunk, including ghost
   points. The number of cells processed in one call is
   (sizes(1)-3)*(sizes(2)-1)*(sizes(3)-1).

   grd, done, and above are scratch arrays. These are point centered 
   arrays. 

*/

typedef struct Edge {
  int vert0, vert1;
} Edge ;

static int hex2tets[6][4] = { {0,6,5,1}, {0,6,1,2}, {0,6,2,3}, {0,6,3,7},
                              {0,6,7,4}, {0,6,4,5} };

static Edge edges[19] = {{0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7},
                           {7,4}, {0,4}, {1,5}, {2,6}, {3,7}, {0,6},
                           {0,5}, {1,6}, {0,2}, {3,6}, {0,7}, {4,6} };

static int tet_edges[6][6] = { {12,5,13,0,14,9}, {12,14,0,15,10,1}, {12,10,15,3,16,2}, 
                           {12,16,3,17,6,11}, {12,6,17,8,18,7}, {12,18,8,13,5,4} };

#ifdef OLD_ISO_STUFF
/* There may be up to two triangles (i.e. 6 edges cut) per tet.
   -1 means "not there". The edge numbering is relative to the tet. To get
   edges relative to the hex, use the edges[] array. */
static int tet_cases[16][6]= { {-1,-1,-1,-1,-1,-1}, {0,3,2,-1,-1,-1}, {0,1,4,-1,-1,-1}, 
                           {3,2,1,3,1,4}, {1,2,5,-1,-1,-1}, {3,5,1,3,1,0}, 
                           {2,5,4,2,4,0}, {3,5,4,-1,-1,-1}, {5,3,4,-1,-1,-1}, 
                           {5,2,0,5,0,4}, {5,3,0,5,0,1}, {5,2,1,-1,-1,-1}, 
                           {2,3,4,2,4,1}, {0,4,1,-1,-1,-1}, {3,0,2,-1,-1,-1}, 
                           {-1,-1,-1,-1,-1,-1} } ;
#endif
/* The iso-surface may also be described using polygons.
   In all cases, the first and second vertices (indices 0 and 1) are those on the first
   face of the tet. The poly_last vector indicates the index of the first
   vertex on the last face.
*/
static int tet_polys[16][4]= { {-1,-1,-1,-1}, {2,0,3,-1}, {0,1,4,-1}, 
                           {2,1,4,3}, {1,2,5,-1}, {1,0,3,5}, 
                           {0,2,5,4}, {5,4,3,-1}, {5,3,4,-1}, 
                           {2,0,4,5}, {0,1,5,3}, {2,1,5,-1}, 
                           {1,2,3,4}, {1,0,4,-1}, {0,2,3,-1}, 
                           {-1,-1,-1,-1} } ;
static int poly_first[16]= { -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, -1 } ;
static int poly_last[16]= { -1, 1, 2, 2, -1, 1, 3, 1, 1, 1, 3, -1, 2, 1, 2, -1 } ;
static int tri_count[16]= {  0, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1,  0} ;

/* The tets are ordered so that each one has a "first face" and a "last face".
   The last face of tet N is the same as ths first face of tet N+1.
   For each of the 16 cases, the edges on the first and last faces have been
   identified and stored in arrays first_edge[] and last_edge[].
   These edges are numbered relative to the tet, so they must be mapped through 
   tet_edges before they can be compared to the edges in the adjacent tet.
*/
static Edge first_edge[16]= { {-1,-1}, {2,0}, {0,1}, {2,1}, {1,2}, {1,0}, {0,2}, 
                            {-1,-1}, {-1,-1}, {2,0}, {0,1}, {2,1}, {1,2}, 
                            {1,0}, {0,2}, {-1,-1} } ;
static Edge last_edge[16]= { {-1,-1}, {0,3}, {4,0}, {4,3}, {-1,-1}, {0,3}, {4,0}, 
                            {4,3}, {3,4}, {0,4}, {3,0}, {-1,-1}, {3,4}, 
                            {0,4}, {3,0}, {-1,-1} } ;
static char face_code[16]= {0x03, 0x05, 0x09, 0x11, 0x22, 0x24, 0x28, 
                            0x30, 0x12, 0x06, 0x0c, 0x18 } ;

/* The worst case is two triangles per tet. Provide space for that, then collapse
   down to triangles with vertices only on the "real" edges */
typedef struct tri_data {
  int edge1, edge2, edge3;
  double x, y, z, nrmx, nrmy, nrmz;
} tri_data ;

/* The triangles will be assembled into polygons. After all tets have been 
   processed, vertices on internal edges will be removed and the polygons
   will be turned back into triangle strips.
   The worst case polygon will have 12 triangles per zone. This can happen when
   the vertices at the end of the main diagonal are both low and all other 
   points are high or vice versa.
   In the worst case there can't be more than 6 polys (most likely
   some polys have to patch together if there are 6 of them).
*/
typedef struct poly_data {
  int edges[36], num_edge, ndx_first, ndx_last;
} poly_data ;

typedef struct tstrip_data {
  long nvert, edges[12];
} tstrip_data ;

typedef struct tstrip_array {
  long nStrip;
  tstrip_data the_tris[6];
} tstrip_array ;

typedef struct zone_tris {
  long nStrip;
  long *lens;
  long *edges;
} zone_tris ;

zone_tris iso_cases[256];
int have_iso_cases= 0;

#ifdef OLD_ISO_STUFF
static tstrip_data the_tris[12];
static int tri_edges_found[12][3];
#endif

static poly_data the_polys[6];
static tstrip_data the_strips[6];
static int num_poly;

char vertflag[8];

#define XYZ_DONE 1
#define GRD_DONE 2
#define VAR_DONE 4
#define V2_DONE 8
#define ABOVE_DONE 16
#define VAR_ABOVE 32
typedef yPoint3D (*GET_XYZ)(long i, long j, long k);
typedef double (*GET_VAR)(long i, long j, long k);
typedef double (*GET_V2)(long i, long j, long k);

static GET_XYZ f_xyz;
static GET_VAR f_var;
static GET_V2 f_v2;

double *cntr_var, *cntr_v2;
double cntr_dx, cntr_dy, cntr_dz, cntr_x0, cntr_y0, cntr_z0;
yPoint3D *cntr_xyz, cntr_grad;
long cntr_iOrigin, cntr_jOrigin, cntr_kOrigin, cntr_iSize, cntr_jSize, cntr_kSize;
long cntr_xy_siz, cntr_x_siz, cntr_nzone, cntr_nexndx;

extern int ycContourTet_array(long make_strip, long sizes[3], double level, 
           double *var, double *var2, yPoint3D *xyz, yPoint3D *grd,  
           unsigned char *flag, TriArrayGrp *triangles);
extern int ycContourTet_array_ndx(long make_strip, long sizes[3], double level, 
           double *var, double *var2, yPoint3D *xyz, yPoint3D *grd,  
           unsigned char *flag, long *ndx, TriVertexGrp *triangles);
extern int ycContourTetHex(double level, long ifirst, long nzone, yPoint3D *xyz, 
           yPoint3D *grad, long hexndx[][8], 
           double *var, double *var2, TriArrayGrp *triangles);
extern int ycGradientChunk(long nx, long xy_siz, long i, long j, long k, long idx, yPoint3D *xyz, 
           double *var, yPoint3D *grd, unsigned char *flag);
extern int ycGetVarXyzCart(long nx, long xy_siz, long i, long j, long k, long idx, yPoint3D *xyz, 
           double *var, unsigned char *flag, GET_XYZ f_xyz, GET_VAR f_var);

extern int ycInitCartGrdPcen(long sizes[3], long offsets[6], 
           double deltas[3], double origin[3],  double *var, double *var2);
extern int ycInitCartGrdZcen(long sizes[3], long offsets[6], 
           double deltas[3], double origin[3],  double *var, double *var2);
extern int ycInitCartGrdPcenNdx(long sizes[3], long offsets[6], 
           double deltas[3], double origin[3],  double *var, double *var2);
extern int ycInitCartGrdZcenNdx(long sizes[3], long offsets[6], 
           double deltas[3], double origin[3],  double *var, double *var2);

extern int ycInitCrvGrdPcen(long sizes[3], long offsets[6], 
           yPoint3D *xyz, double *var, double *var2);
extern int ycInitCrvGrdZcen(long sizes[3], long offsets[6], 
           yPoint3D *xyz, double *var, double *var2);
extern int ycInitCrvGrdPcenNdx(long sizes[3], long offsets[6], 
           yPoint3D *xyz, double *var, double *var2);
extern int ycInitCrvGrdZcenNdx(long sizes[3], long offsets[6], 
           yPoint3D *xyz, double *var, double *var2);

extern int ycContourTet_OneZone(double level, long zon_ndx, int index,
           double s[8], double vars[8], yPoint3D pts[8],
           yPoint3D gradients[8], TriArrayGrp *triangles);

extern double ycContourGrdPcenVar(long i, long j, long k);
extern double ycContourGrdPcenV2(long i, long j, long k);

extern yPoint3D ycContourCartXyz(long i, long j, long k);
extern double ycContourCartZcenVar(long i, long j, long k);
extern double ycContourCartZcenV2(long i, long j, long k);
extern double ycContourCartGrdZcenVar(long i, long j, long k);
extern double ycContourCartGrdZcenV2(long i, long j, long k);

extern double ycContourPcenVar(long i, long j, long k);
extern double ycContourPcenV2(long i, long j, long k);
extern double ycContourPcenAllvar(long i, long j, long k, double *var);

extern yPoint3D ycContourCrvGrdXyz(long i, long j, long k);
extern double ycContourCrvGrdZcenVar(long i, long j, long k);
extern double ycContourCrvGrdZcenV2(long i, long j, long k);
extern double ycContourCrvGrdZcenAllvar(long i, long j, long k, double *var);


extern int ycTetIso_zone(double level, double *var, tstrip_array *t_strips);
extern int tetiso_zone(tstrip_data t_strips[6]);
extern void patch_poly(long pnum, long index, long num, long itet);
extern void patch_2polys(long ip1, long ip2);
extern void assemble_strip(int ndx, int ipn, tstrip_data *t_strips);

#ifdef _DEBUG
#define ASSERT(cond, msg) if(!(cond)) YError(msg);
#else
#define ASSERT(cond, msg) 
#endif

void patch_poly(long pnum, long index, long num, long itet)
{
  long num_edg, nd_old, nd_new, ii;

  nd_old= the_polys[pnum].ndx_last;
  num_edg= the_polys[pnum].num_edge;
  /* Shift vertices of existing polygon up to make room to insert new ones.
     If appending, the for does nothing. */
  for(ii= num_edg+num-1; ii > nd_old+num; ii--) {
    the_polys[pnum].edges[ii]= the_polys[pnum].edges[ii-num];
  }
  nd_new= nd_old+1;
  the_polys[pnum].num_edge += num;
  /* the first edge of the new tri or quad is at the start, so
     1 or 2 vertices need to be inserted into the existing poly
     between the vertices of its current final edge */
  the_polys[pnum].edges[nd_new]= tet_edges[itet][ tet_polys[index][2] ];
  if(num == 2) {
    the_polys[pnum].edges[nd_new+1]= tet_edges[itet][ tet_polys[index][3] ];
  }
  switch(poly_last[index]) {
    case 1:    /* the edge on the last face of the tet starts at tet edge 1 */
      the_polys[pnum].ndx_last= nd_new-1;
      break;
    case 2:    /* the edge on the last face of the tet starts at tet edge 2 */
      the_polys[pnum].ndx_last= nd_new;
      break;
    case 3:    /* the edge on the last face of the tet starts at tet edge 3 */
     the_polys[pnum].ndx_last= nd_new+1;
      break;
    default:   /* Impossible case? No, apparently not. */
      /*      ASSERT(0, "impossible case in patch_poly"); */
      break;
  }
}

void patch_2polys(long ip1, long ip2)
{
  long num_old, num_new, nd_new, ii, jj, nadd;
#ifdef CHEK_POLYS
  long bad= 0;
#endif
  
  /* This function should only be called when ip1 has a first and 
     no last and ip2 has a last but no first.
     ACTUALLY, it appears that this is not true. There are cases
     where the test fails, but the resulting iso-surfaces appear to
     be fine.
     NO EXPLANATION at this time. */
#ifdef CHEK_POLYS
  if(the_polys[ip1].ndx_last >= 0) {
    bad |= 1;
  }
  if(the_polys[ip1].ndx_first != 0) {
    bad |= 2;
  }
  if(the_polys[ip2].ndx_last < 0) {
    bad |= 4;
  }
  if(the_polys[ip2].ndx_first == 0) {
    bad |= 8;
  }
  if(bad) {
    return;  /* should never reach this line */
  }
#endif
  num_old= the_polys[ip1].num_edge;
  nd_new=  the_polys[ip2].ndx_last;
  num_new= the_polys[ip2].num_edge;
  nadd= num_new-2;
  /* Shift vertices of existing polygon up to make room to insert new ones. */
  for(ii= num_old-1; ii > 0; ii--) {
    the_polys[ip1].edges[ii+nadd]= the_polys[ip1].edges[ii];
  }
  /* copy vertices from the second polygon */
  for(ii= 1, jj= (nd_new+2)%num_new; ii < 1+nadd; ii++, jj= (jj+1)%num_new) {
    the_polys[ip1].edges[ii]= the_polys[ip2].edges[jj];
  }
  the_polys[ip1].num_edge += nadd;
  the_polys[ip1].ndx_last= the_polys[ip1].ndx_first;
  /* now remove ip2 (WARNING - this involves a lot of data movement!!) */
  for(ii= ip2; ii < num_poly; ii++) {
    the_polys[ii]= the_polys[ii+1];
  }
  num_poly--;
}

int ycInitCartPcen(long sizes[3], long offsets[6], 
            double deltas[3], double origin[3], double *var, double *var2)
{
  /* This is the version for point-centered data on a Cartesian grid.
     The input arrays do NOT include any guard cells.
     Validate arguments here. */
  if (!var || sizes[0] < 2 || sizes[1] < 2 || sizes[2] < 2) {
    return 0;
  }
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= deltas[0];
  cntr_dy= deltas[1];
  cntr_dz= deltas[2];
  cntr_x0= origin[0];
  cntr_y0= origin[1];
  cntr_z0= origin[2];
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points 
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCartXyz;
  f_var= ycContourPcenVar;
  f_v2= ycContourPcenV2;
  return 1;
}

int ycInitCartGrdPcen(long sizes[3], long offsets[6], 
            double deltas[3], double origin[3], double *var, double *var2)
{
  /* This is the version for point-centered data on a Cartesian grid.
     The input arrays include guard cells on all sides. */
  if (!var || sizes[0] < 4 || sizes[1] < 4 || sizes[2] < 4) {
    return 0;
  }
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= deltas[0];
  cntr_dy= deltas[1];
  cntr_dz= deltas[2];
  cntr_x0= origin[0];
  cntr_y0= origin[1];
  cntr_z0= origin[2];
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of cells including ghost cells
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCartXyz;
  f_var= ycContourGrdPcenVar;
  f_v2= ycContourGrdPcenV2;
  return 1;
}

int ycInitCartZcen(long sizes[3], long offsets[6], 
            double deltas[3], double origin[3], double *var, double *var2)
{
  /* This is the version for zone-centered data on a Cartesian grid.
     The input arrays do not have any guard cells. */
  if (!var || sizes[0] < 3 || sizes[1] < 3 || sizes[2] < 3) {
    return 0;
  }
  cntr_xyz= 0;    /* this is not used for Cartesian grids */
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= deltas[0];
  cntr_dy= deltas[1];
  cntr_dz= deltas[2];
  cntr_x0= origin[0];
  cntr_y0= origin[1];
  cntr_z0= origin[2];
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points including ghost points
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCartXyz;
  f_var= ycContourCartZcenVar;
  f_v2= ycContourCartZcenV2;
  return 1;
}

int ycInitCartGrdZcen(long sizes[3], long offsets[6], 
            double deltas[3], double origin[3], double *var, double *var2)
{
  /* This is the version for zone-centered data on a Cartesian grid.
     The input arrays include guard cells on all sides. */
  if (!var || sizes[0] < 3 || sizes[1] < 3 || sizes[2] < 3) {
    return 0;
  }
  cntr_xyz= 0;    /* this is not used for Cartesian grids */
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= deltas[0];
  cntr_dy= deltas[1];
  cntr_dz= deltas[2];
  cntr_x0= origin[0];
  cntr_y0= origin[1];
  cntr_z0= origin[2];
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points including ghost points
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCartXyz;
  f_var= ycContourCartGrdZcenVar;
  f_v2= ycContourCartGrdZcenV2;
  return 1;
}

int ycInitCartGrdPcenNdx(long sizes[3], long offsets[6], 
            double deltas[3], double origin[3], double *var, double *var2)
{
  /* This is the version for point-centered data on a Cartesian grid.
     The input arrays include guard cells on all sides.
     The result is triangles given by indices into a point list. */
  if (!var || sizes[0] < 4 || sizes[1] < 4 || sizes[2] < 4) {
    return 0;
  }
  cntr_xyz= 0;    /* this is not used for Cartesian grids */
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= deltas[0];
  cntr_dy= deltas[1];
  cntr_dz= deltas[2];
  cntr_x0= origin[0];
  cntr_y0= origin[1];
  cntr_z0= origin[2];
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points including ghost points
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCartXyz;
  f_var= ycContourGrdPcenVar;
  f_v2= ycContourGrdPcenV2;
  return 1;
}

int ycInitCartGrdZcenNdx(long sizes[3], long offsets[6], 
            double deltas[3], double origin[3], double *var, double *var2)
{
  /* This is the version for zone-centered data on a Cartesian grid.
     The input arrays include guard cells on all sides.
     The result is triangles given by indices into a point list. */
  if (!var || sizes[0] < 3 || sizes[1] < 3 || sizes[2] < 3) {
    return 0;
  }
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= deltas[0];
  cntr_dy= deltas[1];
  cntr_dz= deltas[2];
  cntr_x0= origin[0];
  cntr_y0= origin[1];
  cntr_z0= origin[2];
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points including ghost points
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCartXyz;
  f_var= ycContourCartGrdZcenVar;
  f_v2= ycContourCartGrdZcenV2;
  return 1;
}



int ycInitCrvGrdPcen(long sizes[3], long offsets[6], 
            yPoint3D *xyz, double *var, double *var2)
{
  /* This is the version for point-centered data on a Curvilinear grid.
     The input arrays include guard cells on all sides. */
  if (!var || sizes[0] < 4 || sizes[1] < 4 || sizes[2] < 4) {
    return 0;
  }
  cntr_xyz= xyz; 
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= 0.0;
  cntr_dy= 0.0;
  cntr_dz= 0.0;
  cntr_x0= 0.0;
  cntr_y0= 0.0;
  cntr_z0= 0.0;
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points including ghost points
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCrvGrdXyz;
  f_var= ycContourGrdPcenVar;
  f_v2= ycContourGrdPcenV2;
  return 1;
}

int ycInitCrvGrdZcen(long sizes[3], long offsets[6], 
            yPoint3D *xyz, double *var, double *var2)
{
  /* This is the version for zone-centered data on a Curvilinear grid.
     The input arrays include guard cells on all sides. */
  if (!var || sizes[0] < 3 || sizes[1] < 3 || sizes[2] < 3) {
    return 0;
  }
  cntr_xyz= xyz; 
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= 0.0;
  cntr_dy= 0.0;
  cntr_dz= 0.0;
  cntr_x0= 0.0;
  cntr_y0= 0.0;
  cntr_z0= 0.0;
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points including ghost points
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCrvGrdXyz;
  f_var= ycContourCrvGrdZcenVar;
  f_v2= ycContourCrvGrdZcenV2;
  return 1;
}

int ycInitCrvGrdPcenNdx(long sizes[3], long offsets[6], 
            yPoint3D *xyz, double *var, double *var2)
{
  /* This is the version for point-centered data on a Curvilinear grid.
     The input arrays include guard cells on all sides.
     The result is triangles given by indices into a point list. */
  if (!var || sizes[0] < 4 || sizes[1] < 4 || sizes[2] < 4) {
    return 0;
  }
  cntr_xyz= xyz; 
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= 0.0;
  cntr_dy= 0.0;
  cntr_dz= 0.0;
  cntr_x0= 0.0;
  cntr_y0= 0.0;
  cntr_z0= 0.0;
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points including ghost points
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCrvGrdXyz;
  f_var= ycContourGrdPcenVar;
  f_v2= ycContourGrdPcenV2;
  return 1;
}

int ycInitCrvGrdZcenNdx(long sizes[3], long offsets[6], 
            yPoint3D *xyz, double *var, double *var2)
{
  /* This is the version for zone-centered data on a Curvilinear grid.
     The input arrays include guard cells on all sides.
     The result is triangles given by indices into a point list. */
  if (!var || sizes[0] < 3 || sizes[1] < 3 || sizes[2] < 3) {
    return 0;
  }
  cntr_xyz= xyz; 
  cntr_var= var;
  cntr_v2= var2;
  cntr_dx= 0.0;
  cntr_dy= 0.0;
  cntr_dz= 0.0;
  cntr_x0= 0.0;
  cntr_y0= 0.0;
  cntr_z0= 0.0;
  /* on input, the offsets are 1-based, so convert to 0-based now */
  cntr_iOrigin= offsets[0]-1;
  cntr_jOrigin= offsets[1]-1;
  cntr_kOrigin= offsets[2]-1;
  /* iSize and jSize are numbers of points including ghost points
     in the full array. Leave them this way. */
  cntr_iSize= offsets[3];
  cntr_jSize= offsets[4];
  cntr_kSize= offsets[5];
  f_xyz= ycContourCrvGrdXyz;
  f_var= ycContourCrvGrdZcenVar;
  f_v2= ycContourCrvGrdZcenV2;
  return 1;
}


yPoint3D ycContourCartXyz(long i, long j, long k)
{
  yPoint3D xyz_out;

  /* i,j,k are indices relative to the current chunk. Use the 
     location of the corner of this chunk to get xyz. 
     Same function for point centered and zone centered 
     Cartesian grids. Works whether there are guard points or not*/
  xyz_out.x= cntr_x0+i*cntr_dx;
  xyz_out.y= cntr_y0+j*cntr_dy;
  xyz_out.z= cntr_z0+k*cntr_dz;
  return xyz_out;
}

double ycContourPcenAllvar(long i, long j, long k, double *var)
{
  long idxg, i0, j0, k0;

  /* NOTE: This function is for the case of point-centered variables
     without guard cell data. The iso-surface is set up to act as if
     guard cells are present. This means the iso-surface may generate
     indices outside the full array (e.g. a global index of -1).
     Force the indices back in range. */
  /* i,j,k are indices relative to the current chunk. Convert into indices
     into the global array. Same function for Cartesian and Curvilinear grids. 
     For the case of no guard points. 
     For now, do not try to extrapolate outside of the data range.
     This will produce curvature at the boundaries, but is very easy
     to program. If the effect is visually unpleasing, the user should
     generate guard points. */
  /* get indices relative to the full grid */
  i0= i+cntr_iOrigin;
  j0= j+cntr_jOrigin;
  k0= k+cntr_kOrigin;
  if(i0 < 0 || j0 < 0 || k0 < 0 || i0 > cntr_iSize-1 || j0 > cntr_jSize-1
        || k0 > cntr_kSize-1 ) {
    /* the point is outside the grid. */
    if(i0 < 0) i0= 0;
    if(j0 < 0) j0= 0;
    if(k0 < 0) k0= 0;
    if(i0 > cntr_iSize-1) i0= cntr_iSize-1;
    if(j0 > cntr_iSize-1) j0= cntr_jSize-1;
    if(k0 > cntr_iSize-1) k0= cntr_jSize-1;
    idxg= i0+j0*cntr_iSize+k0*cntr_iSize*cntr_jSize;
  } else {
    /* normal case - the point is in the grid */
    idxg= i0+j0*cntr_iSize+k0*cntr_iSize*cntr_jSize;
  }
  return var[idxg];
}

double ycContourPcenVar(long i, long j, long k)
{
  /* call with the contour variable */
  return ycContourPcenAllvar(i, j, k, cntr_var);
}

double ycContourPcenV2(long i, long j, long k)
{
  /* call with the auxiliary variable */
  return ycContourPcenAllvar(i, j, k, cntr_v2);
}

double ycContourGrdPcenVar(long i, long j, long k)
{
  long idxg;

  /* i,j,k are indices relative to the current chunk. Convert into indices
     into the global array. Same function for Cartesian and Curvilinear. */
  idxg= i+cntr_iOrigin+(j+cntr_jOrigin)*cntr_iSize
        +(k+cntr_kOrigin)*cntr_iSize*cntr_jSize;
  return cntr_var[idxg];
}

double ycContourGrdPcenV2(long i, long j, long k)
{
  long idxg;

  /* i,j,k are indices relative to the current chunk. Convert into indices
     into the global array. Same function for Cartesian and Curvilinear. */
  idxg= i+cntr_iOrigin+(j+cntr_jOrigin)*cntr_iSize
        +(k+cntr_kOrigin)*cntr_iSize*cntr_jSize;
  return cntr_v2[idxg];
}

double ycContourCartZcenAllvar(long i, long j, long k, double *var)
{
  long ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi;
  long ijSize= (cntr_iSize-1)*(cntr_jSize-1);
  long var_iSize= (cntr_iSize-1);
  double vcen;

  /* i,j,k are indices relative to the current chunk. Convert into indices
     into the global array. The variable is zone centered so need to interpolate
     except on the boundaries. The contributing cells are ig and ig-1 etc.
     NOTE: the var array is one smaller in each direction than given by 
     cntr_isize etc. */
  ig= i+cntr_iOrigin;
  jg= j+cntr_jOrigin;
  kg= k+cntr_kOrigin;
  if(ig <= 0) {
    ilo= ihi= 0;
  } else if(ig > cntr_iSize-2) {  /* cntr_iSize-2 is the last ZONE */
    ilo= ihi= cntr_iSize-2;
  } else {
    ilo= ig-1;
    ihi= ig;
  }
  if(jg <= 0) {
    jlo= jhi= 0;
  } else if(jg > cntr_jSize-2) {
    jlo= jhi= cntr_jSize-2;
  } else {
    jlo= jg-1;
    jhi= jg;
  }
  if(kg <= 0) {
    klo= khi= 0;
  } else if(kg > cntr_kSize-2) {
    klo= khi= cntr_kSize-2;
  } else {
    klo= kg-1;
    khi= kg;
  }

  vcen=   var[ilo+jlo*var_iSize+klo*ijSize]
        + var[ihi+jlo*var_iSize+klo*ijSize]
        + var[ilo+jhi*var_iSize+klo*ijSize]
        + var[ihi+jhi*var_iSize+klo*ijSize]
        + var[ilo+jlo*var_iSize+khi*ijSize]
        + var[ihi+jlo*var_iSize+khi*ijSize]
        + var[ilo+jhi*var_iSize+khi*ijSize]
        + var[ihi+jhi*var_iSize+khi*ijSize];
  vcen *= 0.125;
  return vcen;
}

double ycContourCartZcenVar(long i, long j, long k)
{
  /* call with the contour variable */
  return ycContourCartZcenAllvar(i, j, k, cntr_var);
}

double ycContourCartZcenV2(long i, long j, long k)
{
  /* call with the auxiliary variable */
  return ycContourCartZcenAllvar(i, j, k, cntr_v2);
}

double ycContourCartGrdZcenVar(long i, long j, long k)
{
  /* call with the contour variable */
  return ycContourCartZcenAllvar(i, j, k, cntr_var);
}

double ycContourCartGrdZcenV2(long i, long j, long k)
{
  /* call with the auxiliary variable */
  return ycContourCartZcenAllvar(i, j, k, cntr_v2);
}


yPoint3D ycContourCrvGrdXyz(long i, long j, long k)
{
  long idxg;

  /* i,j,k are indices relative to the current chunk. Use the 
     location of the corner of this chunk to get xyz. 
     Same function for point centered and zone centered 
     Cartesian grids. */
  idxg= i+cntr_iOrigin+(j+cntr_jOrigin)*cntr_iSize
        +(k+cntr_kOrigin)*cntr_iSize*cntr_jSize;
  return cntr_xyz[idxg];
}

double ycContourCrvGrdZcenAllvar(long i, long j, long k, double *var)
{
  long ig, jg, kg, ilo, ihi, jlo, jhi, klo, khi;
  long ijSize= (cntr_iSize-1)*(cntr_jSize-1);
  long var_iSize= (cntr_iSize-1);
  double vcen;

  /* i,j,k are indices relative to the current chunk. Convert into indices
     into the global array. The variable is zone centered so need to interpolate
     except on the boundaries. The contributing cells are ig and ig-1 etc.
     NOTE: the var array is one smaller in each direction than given by 
     cntr_isize etc. */
  ig= i+cntr_iOrigin;
  jg= j+cntr_jOrigin;
  kg= k+cntr_kOrigin;
  if(ig <= 0) {
    ilo= ihi= 0;
  } else if(ig > cntr_iSize-2) {  /* cntr_iSize-2 is the last ZONE */
    ilo= ihi= cntr_iSize-2;
  } else {
    ilo= ig-1;
    ihi= ig;
  }
  if(jg <= 0) {
    jlo= jhi= 0;
  } else if(jg > cntr_jSize-2) {
    jlo= jhi= cntr_jSize-2;
  } else {
    jlo= jg-1;
    jhi= jg;
  }
  if(kg <= 0) {
    klo= khi= 0;
  } else if(kg > cntr_kSize-2) {
    klo= khi= cntr_kSize-2;
  } else {
    klo= kg-1;
    khi= kg;
  }

  vcen=   var[ilo+jlo*var_iSize+klo*ijSize]
        + var[ihi+jlo*var_iSize+klo*ijSize]
        + var[ilo+jhi*var_iSize+klo*ijSize]
        + var[ihi+jhi*var_iSize+klo*ijSize]
        + var[ilo+jlo*var_iSize+khi*ijSize]
        + var[ihi+jlo*var_iSize+khi*ijSize]
        + var[ilo+jhi*var_iSize+khi*ijSize]
        + var[ihi+jhi*var_iSize+khi*ijSize];
  vcen *= 0.125;
  return vcen;
}

double ycContourCrvGrdZcenVar(long i, long j, long k)
{
  /* call with the contour variable */
  return ycContourCrvGrdZcenAllvar(i, j, k, cntr_var);
}

double ycContourCrvGrdZcenV2(long i, long j, long k)
{
  /* call with the auxiliary variable */
  return ycContourCrvGrdZcenAllvar(i, j, k, cntr_v2);
}

int ycGradientChunk(long nx, long xy_siz, long i, long j, long k, long idx,
                    yPoint3D *xyz, double *var, yPoint3D *grd,
                    unsigned char *flag)
{
  double del2, delx, dely, delz, delv, delv2;
  long ii, jj, kk, id;

  /* WARNING: this doesn't work next to the domain boundary for zone
     centered data. The variable on the guard points is not set
     properly. The consistent solution is to use a different 
     gradient calculation using only zone centered values. */
  /* compute gradients at all 8 corners */
  for(kk= 0; kk <= 1; kk++) {
    for(jj= 0; jj <= 1; jj++) {
      for(ii= 0; ii <= 1; ii++) {
        id= idx+ii+jj*nx+kk*xy_siz;
        if(!(flag[id] & GRD_DONE)) {
          /* i-direction */
          delv= var[id+1]  -var[id-1];
          delx= xyz[id+1].x-xyz[id-1].x;
          dely= xyz[id+1].y-xyz[id-1].y;
          delz= xyz[id+1].z-xyz[id-1].z;
          del2= delx*delx+dely*dely+delz*delz+1.0e-80;
          delv2= delv/del2;
          grd[id].x = delv2*delx;
          grd[id].y = delv2*dely;
          grd[id].z = delv2*delz;

          /* j-direction */
          delv= var[id+nx]  -var[id-nx];
          delx= xyz[id+nx].x-xyz[id-nx].x;
          dely= xyz[id+nx].y-xyz[id-nx].y;
          delz= xyz[id+nx].z-xyz[id-nx].z;
          del2= delx*delx+dely*dely+delz*delz+1.0e-80;
          delv2= delv/del2;
          grd[id].x += delv2*delx;
          grd[id].y += delv2*dely;
          grd[id].z += delv2*delz;

          /* k-direction */
          delv= var[id+xy_siz]  -var[id-xy_siz];
          delx= xyz[id+xy_siz].x-xyz[id-xy_siz].x;
          dely= xyz[id+xy_siz].y-xyz[id-xy_siz].y;
          delz= xyz[id+xy_siz].z-xyz[id-xy_siz].z;
          del2= delx*delx+dely*dely+delz*delz+1.0e-80;
          delv2= delv/del2;
          grd[id].x += delv2*delx;
          grd[id].y += delv2*dely;
          grd[id].z += delv2*delz;
      	  flag[id] |= GRD_DONE;
        }
      }
    }
  }

  return 0;
}

int ycGetVarXyzCart(long nx, long xy_siz, long i, long j, long k, long idx,
                    yPoint3D *xyz, double *var, unsigned char *flag,
                    GET_XYZ f_xyz, GET_VAR f_var)
{
  long id;

  /* idx is an index relative to the current chunk. 
     First make sure point-centered coordinates and variable values are in
     the chunk-sized arrays for all vertices of this cell and all the 
     vertices directly connected to the cell's vertices. */
  /* NOTE: already have variable values at the corners of the cell from 
     deciding which points were above and below the contour level. */
  id= idx     -xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j,  k-1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i,  j,  k-1); flag[id] |= VAR_DONE;}
  id= idx+1   -xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j,  k-1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+1,j,  k-1); flag[id] |= VAR_DONE;}

  id= idx  +nx-xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j+1,k-1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i,  j+1,k-1); flag[id] |= VAR_DONE;}
  id= idx+1+nx-xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j+1,k-1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+1,j+1,k-1); flag[id] |= VAR_DONE;}


  id= idx  -nx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j-1,k  ); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i,  j-1,k  ); flag[id] |= VAR_DONE;}
  id= idx+1-nx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j-1,k  ); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+1,j-1,k  ); flag[id] |= VAR_DONE;}

  id= idx-1;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i-1,j,  k  ); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i-1,j,  k  ); flag[id] |= VAR_DONE;}
  id= idx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j,  k  ); flag[id] |= XYZ_DONE;}
  id= idx+1;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j,  k  ); flag[id] |= XYZ_DONE;}
  id= idx+2;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+2,j,  k  ); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+2,j,  k  ); flag[id] |= VAR_DONE;}

  id= idx-1+nx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i-1,j+1,k  ); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i-1,j+1,k  ); flag[id] |= VAR_DONE;}
  id= idx  +nx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j+1,k  ); flag[id] |= XYZ_DONE;}
  id= idx+1+nx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j+1,k  ); flag[id] |= XYZ_DONE;}
  id= idx+2+nx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+2,j+1,k  ); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+2,j+1,k  ); flag[id] |= VAR_DONE;}

  id= idx  +2*nx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j+2,k  ); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i,  j+2,k  ); flag[id] |= VAR_DONE;}
  id= idx+1+2*nx;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j+2,k  ); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+1,j+2,k  ); flag[id] |= VAR_DONE;}


  id= idx  -nx+xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j-1,k+1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i,  j-1,k+1); flag[id] |= VAR_DONE;}
  id= idx+1-nx+xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j-1,k+1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+1,j-1,k+1); flag[id] |= VAR_DONE;}

  id= idx-1   +xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i-1,j  ,k+1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i-1,j  ,k+1); flag[id] |= VAR_DONE;}
  id= idx     +xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j  ,k+1); flag[id] |= XYZ_DONE;}
  id= idx+1   +xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j  ,k+1); flag[id] |= XYZ_DONE;}
  id= idx+2   +xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+2,j  ,k+1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+2,j  ,k+1); flag[id] |= VAR_DONE;}

  id= idx-1+nx+xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i-1,j+1,k+1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i-1,j+1,k+1); flag[id] |= VAR_DONE;}
  id= idx  +nx+xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i  ,j+1,k+1); flag[id] |= XYZ_DONE;}
  id= idx+1+nx+xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j+1,k+1); flag[id] |= XYZ_DONE;}
  id= idx+2+nx+xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+2,j+1,k+1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+2,j+1,k+1); flag[id] |= VAR_DONE;}

  id= idx  +2*nx+xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j+2,k+1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i,  j+2,k+1); flag[id] |= VAR_DONE;}
  id= idx+1+2*nx+xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j+2,k+1); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+1,j+2,k+1); flag[id] |= VAR_DONE;}


  id= idx     +2*xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i,  j,  k+2); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i,  j,  k+2); flag[id] |= VAR_DONE;}
  id= idx+1   +2*xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j,  k+2); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+1,j,  k+2); flag[id] |= VAR_DONE;}

  id= idx  +nx+2*xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i  ,j+1,k+2); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i,  j+1,k+2); flag[id] |= VAR_DONE;}
  id= idx+1+nx+2*xy_siz;
  if(!(flag[id] & XYZ_DONE)) {xyz[id] = (*f_xyz)(i+1,j+1,k+2); flag[id] |= XYZ_DONE;}
  if(!(flag[id] & VAR_DONE)) {var[id] = (*f_var)(i+1,j+1,k+2); flag[id] |= VAR_DONE;}

  return 0;
}


int ycContourTet_OneZone(double level, long zonidx, int index, double s[8],
                         double vars[8], yPoint3D pts[8], 
                         yPoint3D gradients[8], TriArrayGrp *triangles)
{
  int ie, vert0, vert1;
  Edge the_edg;
  long *cellIDs, numTri, numVert, ipn, nedg, indbase;
  double t, v1, v2, *v2vals;
  yPoint3D *xyzverts, *normals, n, x1, x2, n1, n2;

  /* The caller allocates storage for the triangles making
     up the iso-surface. 

     The input is a hexahedral zone. 
     s is the variable values.
     vars is the auxiliary variable (optional).
     pts is the coordinates of the corner of the hex.
     gradients is the gradient of s at the corners.
     index is the 8 bit code for the high and low corners.
  */

  cellIDs = triangles->cellIDs;
  xyzverts = triangles->xyzverts;
  normals = triangles->normals;
  v2vals= triangles->var2;
  numTri= triangles->numTri;

  /*
     Generate triangles and point gradients using a 6 tet decomposition.
     The contour level is known to cut this zone.
  */  
  num_poly= iso_cases[index].nStrip;
  /* For all tri-strips, find the intersection points.
     Break them back up into triangles for now. */
  for(ipn= 0, indbase= 0; ipn < num_poly; ipn++) {
    int phase= 1, ij, ind;
    nedg= (iso_cases[index].lens)[ipn];
    for(ie= 0; ie < nedg-2; ie++) {
      for(ij= 0; ij < 3; ij++) {
        if(phase == 0) {
          ind= ie+ij;
        } else {
          ind= ie+2-ij;
        }
        /* get the two vertices at the ends of the edge on which the
           point lies */
        the_edg= edges[ (iso_cases[index].edges)[indbase+ind] ];
        vert0= the_edg.vert0;
        vert1= the_edg.vert1;
        t = (double)(level - s[vert0]) / (s[vert1] - s[vert0]);
        x1 = pts[vert0];
        x2 = pts[vert1];
        n1 = gradients[vert0];
        n2 = gradients[vert1];
        numVert= 3*numTri+ij;
        xyzverts[numVert].x = x1.x + t * (x2.x - x1.x);
        xyzverts[numVert].y = x1.y + t * (x2.y - x1.y);
        xyzverts[numVert].z = x1.z + t * (x2.z - x1.z);
        if(vars) {
          v1 = vars[vert0];
          v2 = vars[vert1];
          v2vals[numVert] = v1 + t * (v2-v1);
        }
        n.x = n1.x + t * (n2.x - n1.x);
        n.y = n1.y + t * (n2.y - n1.y);
        n.z = n1.z + t * (n2.z - n1.z);
        ycNormalize(&n);
        normals[numVert]= n;
      }
      cellIDs[numTri]= zonidx;  /* the C-style index of this cell */
      numTri++;
      if(phase == 0) phase= 1;
      else phase= 0;
    }
    indbase += nedg;
  }
  triangles->numTri= numTri;
  return 1;
}


int ycContourTetHex(double level, long ifirst, long nzone, yPoint3D *xyz, 
            yPoint3D *grad, long hexndx[][8], 
            double *var, double *var2, TriArrayGrp *triangles)
{
  int mask, index, ii, nlo, nhi;
  long i, ndx0, ndx1, ndx2, ndx3, ndx4, ndx5, ndx6, ndx7;
  double s[8], vars[8];
  yPoint3D pts[8], gradients[8];

  /* This is the version for point-centered data on an arbitrarily
     connected grid of hexes.
     hexndx is an 8 by nzone array of indices into the xyz, 
     grad, var, and var2 arrays. 
     Do not use any chunking. The caller is urged to arrange the 
     ordering of nexndx so that zones that appear together in hexndx 
     are nearby spatially.
     To avoid making triangles very large, pass chunks of hexndx. 
     The input does NOT include guard cells because the gradient
     is supplied. */

  /* The caller allocates storage for the triangles making
     up the iso-surface. 
  */
  if (!var || nzone < 1) {
    return 0;
  }

  /* The input is an array of hexahedral zones. Each zone is specified
     by 8 indices into arrays of coordinates, gradients, variables, 
     and (optionally) auxiliary variables.
  */
  ifirst--;  /* convert from yorick-style index to C-style */

  triangles->numTri= 0;

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /*
     Traverse all cells, generating triangles and point gradients
     using a 6 tet decomposition of each hexahedron.
  */  
  for (i= ifirst; i < ifirst+nzone; i++) {
    /* Get data values at the corners. Note that this circulates around 
       the zone, which is not in storage order. Note also that hexndx 
       has yorick style indices. */
    ndx0 = hexndx[i][0]-1;
    ndx1 = hexndx[i][1]-1;
    ndx2 = hexndx[i][2]-1;
    ndx3 = hexndx[i][3]-1;
    ndx4 = hexndx[i][4]-1;
    ndx5 = hexndx[i][5]-1;
    ndx6 = hexndx[i][6]-1;
    ndx7 = hexndx[i][7]-1;
    s[0] = var[ndx0];
    s[1] = var[ndx1];
    s[2] = var[ndx2];
    s[3] = var[ndx3];
    s[4] = var[ndx4];
    s[5] = var[ndx5];
    s[6] = var[ndx6];
    s[7] = var[ndx7];

    /* Build the case table */
    for ( ii=0, nlo= 0, nhi= 0, mask= 1, index= 0; ii < 8; ii++, mask += mask) {
      if ( s[ii] >= level ) {
        index |= mask;
        nhi++;
      } else {
        nlo++;
      }
    }

    if ( nlo == 0 || nhi == 0 ) continue; /* no surface */

    /* pick up the coordinates of the corners of the cell */
    pts[0] = xyz[ndx0];
    pts[1] = xyz[ndx1];
    pts[2] = xyz[ndx2];
    pts[3] = xyz[ndx3];
    pts[4] = xyz[ndx4];
    pts[5] = xyz[ndx5];
    pts[6] = xyz[ndx6];
    pts[7] = xyz[ndx7];

    /* pick up the gradients at the corners of the cell */
    gradients[0] = grad[ndx0];
    gradients[1] = grad[ndx1];
    gradients[2] = grad[ndx2];
    gradients[3] = grad[ndx3];
    gradients[4] = grad[ndx4];
    gradients[5] = grad[ndx5];
    gradients[6] = grad[ndx6];
    gradients[7] = grad[ndx7];

    if(var2) {
      vars[0] = var[ndx0];
      vars[1] = var[ndx1];
      vars[2] = var[ndx2];
      vars[3] = var[ndx3];
      vars[4] = var[ndx4];
      vars[5] = var[ndx5];
      vars[6] = var[ndx6];
      vars[7] = var[ndx7];
      ycContourTet_OneZone(level, i, index, s, vars, pts, gradients, triangles);
    } else {
      ycContourTet_OneZone(level, i, index, s, 0,    pts, gradients, triangles);
    }

  } /* for i */

  if(triangles->numTri > 0) return 1;
  else return 0;
}


int ycContourTet_array(long make_strip, long sizes[3], 
                       double level, double *var, double *var2, yPoint3D *xyz,
                       yPoint3D *grd, unsigned char *flag,
                       TriArrayGrp *triangles)
{
  int i, j, k, xy_siz, mask, index;
  int jOffset, kOffset, idxv, idxcg, idxnow, ii, nlo, nhi;
  long numTri, ihi, jhi, khi, *cellIDs, *nTris, *triStart;
  long nx= sizes[0];
  long ny= sizes[1];
  long nz= sizes[2];
  double s[8], vars[8], *v2vals;
  yPoint3D *xyzverts, *normals, pts[8], gradients[8];

  /* The caller allocates storage for the triangles making
     up the iso-surface. 

     The input array has dimensions sizes[3].
     There are ghost cells on all sides so that gradients
     can be computed without testing for boundaries.
     This means there are sizes[i]-3 real cells in each direction.
     NOTE: origin is the coordinate of the 0,0,0 point.
  */

  cellIDs = triangles->cellIDs;
  xyzverts = triangles->xyzverts;
  normals = triangles->normals;
  v2vals= triangles->var2;
  triStart= triangles->triStart;
  nTris = triangles->nTris;
  xy_siz = nx * ny;
  /* set stride for y and z in a chunk array. */
  cntr_x_siz=  nx;
  cntr_xy_siz= nx*ny;

  /* ihi, jhi, khi are the indices of the upper vertex
     to be processed in this call */
  ihi= nx-2;
  jhi= ny-2;
  khi= nz-2;
  numTri= 0;
  /* This loop runs over all vertices in the chunk, including
     guard points. */
  for(k= 0; k < nz; k++) {
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        idxv = i + j*nx + k*xy_siz;
        flag[idxv]= 0;
      }
    }
  }

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /*
     Traverse all cells, generating triangles and point gradients
     using a 6 tet decomposition of each hexahedron.
     The upper loop limit is a cell index, not a vertex index,
     so need a less than.
  */  
  for ( k=1; k < khi; k++) {
    kOffset = k*xy_siz;
    for ( j=1; j < jhi; j++) {
      jOffset = j*nx;
      for ( i=1; i < ihi; i++) {
        /* get data values at the corners. idxv is a zero-based
           POINT index */
        idxv = i + jOffset + kOffset;
        if(make_strip) {
          /* record the first triangle for this cell */
          triStart[idxv]= numTri;
          nTris[idxv]= 0;
        }

        /* get data values at the corners */
        idxnow= idxv            ; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i,  j,  k  ); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv+1          ; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i+1,j,  k  ); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv  +nx       ; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i,  j+1,k  ); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv+1+nx       ; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i+1,j+1,k  ); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv     +xy_siz; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i,  j,  k+1); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv+1   +xy_siz; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i+1,j,  k+1); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv  +nx+xy_siz; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i,  j+1,k+1); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv+1+nx+xy_siz; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i+1,j+1,k+1); flag[idxnow] |= VAR_DONE;}
        /* note that this circulates around the zone, not counts up in i,j,k */
        s[0] = var[idxv];
        s[1] = var[idxv+1];
        s[2] = var[idxv+1+nx];
        s[3] = var[idxv  +nx];
        s[4] = var[idxv     +xy_siz];
        s[5] = var[idxv+1   +xy_siz];
        s[6] = var[idxv+1+nx+xy_siz];
        s[7] = var[idxv  +nx+xy_siz];

        /* Build the case table */
        for ( ii=0, nlo= 0, nhi= 0, mask= 1, index= 0; ii < 8; ii++, mask += mask) {
          if ( s[ii] >= level ) {
            index |= mask;
            nhi++;
          } else {
            nlo++;
          }
        }

        if ( nlo == 0 || nhi == 0 ) continue; /* no surface */
        if(cntr_v2) {
          idxnow= idxv            ; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i,  j,  k  ); flag[idxnow] |= V2_DONE;}
          idxnow= idxv+1          ; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i+1,j,  k  ); flag[idxnow] |= V2_DONE;}
          idxnow= idxv  +nx       ; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i,  j+1,k  ); flag[idxnow] |= V2_DONE;}
          idxnow= idxv+1+nx       ; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i+1,j+1,k  ); flag[idxnow] |= V2_DONE;}
          idxnow= idxv     +xy_siz; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i,  j,  k+1); flag[idxnow] |= V2_DONE;}
          idxnow= idxv+1   +xy_siz; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i+1,j,  k+1); flag[idxnow] |= V2_DONE;}
          idxnow= idxv  +nx+xy_siz; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i,  j+1,k+1); flag[idxnow] |= V2_DONE;}
          idxnow= idxv+1+nx+xy_siz; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i+1,j+1,k+1); flag[idxnow] |= V2_DONE;}
          vars[0] = var2[idxv];
          vars[1] = var2[idxv+1];
          vars[2] = var2[idxv+1+nx];
          vars[3] = var2[idxv  +nx];
          vars[4] = var2[idxv     +xy_siz];
          vars[5] = var2[idxv+1   +xy_siz];
          vars[6] = var2[idxv+1+nx+xy_siz];
          vars[7] = var2[idxv  +nx+xy_siz];
        }

        /* Get coordinates and variable values necessary for computing the gradient. 
           After this call will have coordinates for all corners of this cell. */
        ycGetVarXyzCart(nx, xy_siz, i, j, k, idxv, xyz, var, flag, f_xyz, f_var);
        /* pick up the coordinates of the corners of the cell */
        pts[0] = xyz[idxv];
        pts[1] = xyz[idxv+1];
        pts[2] = xyz[idxv+1+nx];
        pts[3] = xyz[idxv  +nx];
        pts[4] = xyz[idxv     +xy_siz];
        pts[5] = xyz[idxv+1   +xy_siz];
        pts[6] = xyz[idxv+1+nx+xy_siz];
        pts[7] = xyz[idxv  +nx+xy_siz];
        /* Compute gradients (if not already known) using values on the chunk, not 
           from the global array. */
        ycGradientChunk(nx, xy_siz, i, j, k, idxv, xyz, var, grd, flag);
        gradients[0] = grd[idxv];
        gradients[1] = grd[idxv+1];
        gradients[2] = grd[idxv+1+nx];
        gradients[3] = grd[idxv  +nx];
        gradients[4] = grd[idxv     +xy_siz];
        gradients[5] = grd[idxv+1   +xy_siz];
        gradients[6] = grd[idxv+1+nx+xy_siz];
        gradients[7] = grd[idxv  +nx+xy_siz];

        num_poly= iso_cases[index].nStrip;
        /* For all tri-strips, find the intersection points.
           Break them back up into triangles for now. */
        idxcg= i+cntr_iOrigin+(j+cntr_jOrigin)*(cntr_iSize-1)+(k+cntr_kOrigin)*(cntr_iSize-1)*(cntr_jSize-1);
        extract_tris_tet(index, idxcg, cntr_v2, &numTri, level, s, pts,
                gradients, vars, cellIDs, xyzverts, normals, v2vals);
      } /* for i */
    } /* for j */
  } /* for k */

  triangles->numTri= numTri;
  if(numTri > 0) return 1;
  else return 0;
}

int ycContourTet_array_ndx(long make_strip, long sizes[3], 
                           double level, double *var, double *var2,
                           yPoint3D *xyz, yPoint3D *grd, unsigned char *flag,
                           long *ndx, TriVertexGrp *triangles)
{
  int i, j, k, xy_siz, mask, index, num_idx;
  int jOffset, kOffset, idxv, idxc, idxcg, idxnow, ii, nlo, nhi;
  long numTri, numEdg, ihi, jhi, khi, *cellIDs, *nTris, *ptndx;
  long *triStart, edg_offset[12];
  long nx= sizes[0];
  long ny= sizes[1];
  long nz= sizes[2];
  double s[8], vars[8], *v2vals;
  yPoint3D *xyzverts, *normals, pts[8], gradients[8];
  
  /* The caller allocates storage for the triangles making
     up the iso-surface. 
     The triangles in the output are given by indices into vectors of 
     vertices, gradients, and associated variables. 
  */

  /* The input array has dimensions sizes[3].
     There are ghost cells on all sides so that gradients
     can be computed without testing for boundaries.
     This means there are sizes[i]-3 real cells in each direction.
  */

  cellIDs = triangles->cellIDs;
  xyzverts = triangles->xyzverts;
  normals = triangles->normals;
  ptndx= triangles->ptndx;
  v2vals= triangles->var2;
  triStart= triangles->triStart;
  nTris = triangles->nTris;
  xy_siz = nx * ny;
  /* set stride for y and z in a chunk array. */
  cntr_x_siz=  nx;
  cntr_xy_siz= nx*ny;
  edg_offset[0]= 0;
  edg_offset[1]= 4;
  edg_offset[2]= 3*nx;
  edg_offset[3]= 1;
  edg_offset[4]= 3*xy_siz;
  edg_offset[5]= 3*xy_siz+4;
  edg_offset[6]= 3*xy_siz+3*nx;
  edg_offset[7]= 3*xy_siz+1;
  edg_offset[8]= 2;
  edg_offset[9]= 5;
  edg_offset[10]= 3*nx+5;
  edg_offset[11]= 3*nx+2;

  /* ihi, jhi, khi are the indices of the upper vertex
     to be processed in this call */
  ihi= nx-2;
  jhi= ny-2;
  khi= nz-2;
  /* This loop runs over all vertices in the chunk, including
     guard points. */
  for(k= 0; k < nz; k++) {
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        idxv = i + j*nx + k*xy_siz;
        flag[idxv]= 0;
      }
    }
  }
  num_idx= 3*nx*ny*nz;
  for(i= 0; i < num_idx; i++) {
    ndx[i]= -1;
  }
  numTri= 0;
  numEdg= 0;  /* number of unique edges cut by the iso-surface */

  /* prepare the case table for the iso-surfaces, if not already done */
  if(!have_iso_cases) ycPrepIsoTet();

  /*
     Traverse all cells, generating triangles and point gradients
     using a 6 tet decomposition of each hexahedron.
     The upper loop limit is a cell index, not a vertex index,
     so need a less than.
  */  
  for ( k=1; k < khi; k++) {
    kOffset = k*xy_siz;
    for ( j=1; j < jhi; j++) {
      jOffset = j*nx;
      for ( i=1; i < ihi; i++) {
        /* get data values at the corners. idxv is a zero-based
           POINT index */
        idxv = i + jOffset + kOffset;
        idxc = i + j*(nx-1) + k*(nx-1)*(ny-1);
        if(make_strip) {
          /* record the first triangle for this cell */
          triStart[idxv]= numTri;
          nTris[idxv]= 0;
        }

        /* get data values at the corners */
        idxnow= idxv            ; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i,  j,  k  ); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv+1          ; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i+1,j,  k  ); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv  +nx       ; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i,  j+1,k  ); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv+1+nx       ; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i+1,j+1,k  ); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv     +xy_siz; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i,  j,  k+1); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv+1   +xy_siz; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i+1,j,  k+1); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv  +nx+xy_siz; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i,  j+1,k+1); flag[idxnow] |= VAR_DONE;}
        idxnow= idxv+1+nx+xy_siz; if(!(flag[idxnow] & VAR_DONE)) {var[idxnow] = (*f_var)(i+1,j+1,k+1); flag[idxnow] |= VAR_DONE;}
        /* note that this circulates around the zone, not counts up in i,j,k */
        s[0] = var[idxv];
        s[1] = var[idxv+1];
        s[2] = var[idxv+1+nx];
        s[3] = var[idxv  +nx];
        s[4] = var[idxv     +xy_siz];
        s[5] = var[idxv+1   +xy_siz];
        s[6] = var[idxv+1+nx+xy_siz];
        s[7] = var[idxv  +nx+xy_siz];

        /* Build the case table */
        for ( ii=0, nlo= 0, nhi= 0, mask= 1, index= 0; ii < 8; ii++, mask += mask) {
          if ( s[ii] >= level ) {
            index |= mask;
            nhi++;
          } else {
            nlo++;
          }
        }

        if ( nlo == 0 || nhi == 0 ) continue; /* no surface */
        if(cntr_v2) {
          idxnow= idxv            ; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i,  j,  k  ); flag[idxnow] |= V2_DONE;}
          idxnow= idxv+1          ; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i+1,j,  k  ); flag[idxnow] |= V2_DONE;}
          idxnow= idxv  +nx       ; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i,  j+1,k  ); flag[idxnow] |= V2_DONE;}
          idxnow= idxv+1+nx       ; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i+1,j+1,k  ); flag[idxnow] |= V2_DONE;}
          idxnow= idxv     +xy_siz; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i,  j,  k+1); flag[idxnow] |= V2_DONE;}
          idxnow= idxv+1   +xy_siz; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i+1,j,  k+1); flag[idxnow] |= V2_DONE;}
          idxnow= idxv  +nx+xy_siz; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i,  j+1,k+1); flag[idxnow] |= V2_DONE;}
          idxnow= idxv+1+nx+xy_siz; if(!(flag[idxnow] & V2_DONE)) {var2[idxnow] = (*f_v2)(i+1,j+1,k+1); flag[idxnow] |= V2_DONE;}
          vars[0] = var2[idxv];
          vars[1] = var2[idxv+1];
          vars[2] = var2[idxv+1+nx];
          vars[3] = var2[idxv  +nx];
          vars[4] = var2[idxv     +xy_siz];
          vars[5] = var2[idxv+1   +xy_siz];
          vars[6] = var2[idxv+1+nx+xy_siz];
          vars[7] = var2[idxv  +nx+xy_siz];
        }

        /* Get coordinates and variable values necessary for computing the gradient. 
           After this call will have coordinates for all corners of this cell. */
        ycGetVarXyzCart(nx, xy_siz, i, j, k, idxv, xyz, var, flag, f_xyz, f_var);
        /* pick up the coordinates of the corners of the cell */
        pts[0] = xyz[idxv];
        pts[1] = xyz[idxv+1];
        pts[2] = xyz[idxv+1+nx];
        pts[3] = xyz[idxv  +nx];
        pts[4] = xyz[idxv     +xy_siz];
        pts[5] = xyz[idxv+1   +xy_siz];
        pts[6] = xyz[idxv+1+nx+xy_siz];
        pts[7] = xyz[idxv  +nx+xy_siz];
        /* Compute gradients (if not already known) using values on the chunk, not 
           from the global array. */
        ycGradientChunk(nx, xy_siz, i, j, k, idxv, xyz, var, grd, flag);
        gradients[0] = grd[idxv];
        gradients[1] = grd[idxv+1];
        gradients[2] = grd[idxv+1+nx];
        gradients[3] = grd[idxv  +nx];
        gradients[4] = grd[idxv     +xy_siz];
        gradients[5] = grd[idxv+1   +xy_siz];
        gradients[6] = grd[idxv+1+nx+xy_siz];
        gradients[7] = grd[idxv  +nx+xy_siz];

        num_poly= iso_cases[index].nStrip;
        /* For all tri-strips, find the intersection points.
           Break them back up into triangles for now. */
        idxcg= i+cntr_iOrigin+(j+cntr_jOrigin)*(cntr_iSize-1)+(k+cntr_kOrigin)*(cntr_iSize-1)*(cntr_jSize-1);
        extract_tris_tet_ndx(index, idxv, idxcg, cntr_v2, &numTri, &numEdg, level, s, pts,
                gradients, vars, edg_offset, cellIDs, ptndx, ndx, xyzverts, normals, v2vals);
      } /* for i */
    } /* for j */
  } /* for k */

  triangles->numTri= numTri;
  triangles->numEdg= numEdg;
  if(numTri > 0) return 1;
  else return 0;
}

int tetiso_zone(tstrip_data t_strips[6])
{
  int poly_add, first_vert;
  int index, mask, ntri, num, ie;
  Edge first;
  long itet, vrt, ip, ipn, nedg, vrt0, vrt1, kk;

  ntri= 0; /* no triangles found so far in this zone */
  num_poly= 0; /* no polygons constructed so far for this cell */
  /* run through all tethedra for this cell. */
  for(itet= 0; itet < 6; itet++) {
    for(vrt= 0, mask= 1, index= 0; vrt < 4; vrt++, mask *= 2) {
      if(vertflag[hex2tets[itet][vrt]]) index |= mask;
    }
    /* The tet has been identified as one of the 16 possible cases */
    num= tri_count[index];
    if(num <= 0) {
      /* no polygons in this tet, so proceed to the next tet */
      continue;
    }
    ntri += num;
    /* There is only one polygon per tet, although there might
       be two triangles. Use hex relative edge numbers. */
    first.vert0= tet_edges[itet][ first_edge[index].vert0 ];
    first.vert1= tet_edges[itet][ first_edge[index].vert1 ];
    poly_add= 1;
    if(num_poly > 0) {
      /* determine whether this triangle can be patched onto 
	     a previous polygon */
      for(ip= 0; ip < num_poly; ip++) {
        vrt0= the_polys[ip].edges[ the_polys[ip].ndx_last ];
        vrt1= the_polys[ip].edges[ (the_polys[ip].ndx_last+1) 
              % the_polys[ip].num_edge ];
        if(vrt0 == first.vert1 && vrt1 == first.vert0) {
          /* patch this triangle onto the existing polygon */
          poly_add= 0;  /* don't make a new polygon */
          patch_poly(ip, index, num, itet);
          /* the triangle has been matched, so break out of the for loop */
          break;
        }
      }
    }
    if(poly_add) {
      /* The new polygon doesn't patch onto an old one.
         Keep the index into this polygon of the initial vertex
         of the last edge of the polygon.
         The first index is really a flag for whether the poly
         starts on the first face of the first tet. */
      the_polys[num_poly].ndx_first= poly_first[index];
      the_polys[num_poly].ndx_last= poly_last[index];
      the_polys[num_poly].edges[0]= tet_edges[itet][ tet_polys[index][0] ];
      the_polys[num_poly].edges[1]= tet_edges[itet][ tet_polys[index][1] ];
      the_polys[num_poly].edges[2]= tet_edges[itet][ tet_polys[index][2] ];
      if(num == 1) {
        the_polys[num_poly].num_edge= 3;
      } else {
        the_polys[num_poly].num_edge= 4;
        the_polys[num_poly].edges[3]= tet_edges[itet][ tet_polys[index][3] ];
      }
      num_poly++;
    }
  }
  /* Match polygons with an open last to those with an open first. */
  for(ip= 0; ip < num_poly-1; ip++) {
    poly_data *oldpoly= &(the_polys[ip]);
    long ip2, oldlast= oldpoly->ndx_last;
    if(oldlast >= 0 && oldpoly->ndx_first != 0) {
      /* the last edge of this polygon may patch onto
         the first edge of another polygon. */
      for(ip2= ip+1; ip2 < num_poly; ip2++) {
        if(the_polys[ip2].edges[0] != 0) continue;
        if(the_polys[ip2].edges[1] == oldpoly->edges[oldlast] &&
          the_polys[ip2].edges[0] == oldpoly->edges[ (oldlast+1)%oldpoly->num_edge] ) {
          /* patch the two polygons together, then adjust first and last */
          patch_2polys(ip2, ip);
          break;
        }
      }
    } else if(oldpoly->ndx_first == 0 && oldpoly->ndx_last < 0) {
      /* the first edge of this polygon may patch onto
         the last edge of another polygon. */
      for(ip2= ip+1; ip2 < num_poly; ip2++) {
        long ip2last= the_polys[ip2].ndx_last;
        if(ip2last < 0) continue;
        if(the_polys[ip2].edges[ip2last] == oldpoly->edges[1] &&
           the_polys[ip2].edges[(ip2last+1)%the_polys[ip2].num_edge] 
           == oldpoly->edges[0] ) {
          /* patch the two polygons together, then adjust first and last */
          patch_2polys(ip, ip2);
          break;
        }
      }
    } else if(oldpoly->ndx_first == 0 && oldpoly->ndx_last >= 0) {
      /* This polygon has both a first and last edge. See if it
         patches onto any other polygon. */
      for(ip2= ip+1; ip2 < num_poly; ip2++) {
        long ip2last= the_polys[ip2].ndx_last;
        if(ip2last < 0) continue;
        if(the_polys[ip2].edges[ip2last] == oldpoly->edges[1] &&
           the_polys[ip2].edges[(ip2last+1)%the_polys[ip2].num_edge] 
           == oldpoly->edges[0] ) {
          /* patch the two polygons together, then adjust first and last */
          patch_2polys(ip, ip2);
          break;
        }
      }
    }
  }
  /* find all open polygons and try to close them */
  for(ip= 0; ip < num_poly; ip++) {
    if(the_polys[ip].ndx_last == the_polys[ip].ndx_first) {
      /* this should be the result of patching two polys together in the 
         previous loop */
      continue;
    }
    if(the_polys[ip].ndx_last >= 0) {
      /* the last edge of this polygon should patch onto
         the first edge */
      if(the_polys[ip].edges[0] == the_polys[ip].edges[2]) {
        /* remove the first two vertices */
        for(kk= 2; kk < the_polys[ip].num_edge; kk++) {
          the_polys[ip].edges[kk-2]= the_polys[ip].edges[kk];
        }
        the_polys[ip].num_edge -= 2;
      } else if(the_polys[ip].edges[0] == the_polys[ip].edges[the_polys[ip].num_edge-2]) {
        /* remove the final two points */
        the_polys[ip].num_edge -= 2;
      } else {
        /* impossible error?? */
      }
    }
  }
  /* For all polygons, remove all points that are not on real edges.  */
  for(ipn= 0; ipn < num_poly; ipn++) {
    int iout;
    nedg= the_polys[ipn].num_edge;
    for(ie= 0, iout= 0; ie < nedg; ie++) {
      if(the_polys[ipn].edges[ie] >= 12) continue;
      if(iout < ie) {
        /* move the real edge down in the polygon */
        the_polys[ipn].edges[iout]= the_polys[ipn].edges[ie];
      }
      iout++;
    }
    the_polys[ipn].num_edge= iout;
  }
  /* For all polygons, remove duplicate first and last points.  */
  for(ipn= 0; ipn < num_poly; ipn++) {
    nedg= the_polys[ipn].num_edge;
    if(the_polys[ipn].edges[0] == the_polys[ipn].edges[nedg-1]) the_polys[ipn].num_edge--;
  }
  /* Convert all polygons into triangle strips. Special steps are
      required for cases where all 4 edges on a face are cut by the 
     iso-surface. */
  for(ipn= 0; ipn < num_poly; ipn++) {
    int fn, msk, fnum[6], nfours, last4;
    nfours= 0;
    nedg= the_polys[ipn].num_edge;
    for(fn= 0, msk= 1; fn < 6; fn++) {
      num= 0;
      for(ie= 0; ie < nedg; ie++) {
        if(face_code[ the_polys[ipn].edges[ie] ] & msk) num++;
      }
      if(num >= 3) {
        nfours++;
        last4= msk;
#undef BADNUM_CHEK
#ifdef BADNUM_CHEK
        /*        ASSERT( (num != 3), "bad count in tetiso_zone"); */
        if(num == 3) {
          int badone= 1;
        }
#endif
      }
      fnum[fn]= num;
      msk += msk;
    }
    switch(nedg) {
      case 0:
      default:
        break;
      case 3:
      case 4:
      case 5:
        /* tri-strip is simple because there aren't more than two edges 
           cut on any one face. */
        assemble_strip(0, ipn, t_strips);
        break;
      case 6:
        if(nfours > 0) {
          /* Find the first edge not on the face with four edges cut.
             The tri-strip should be started at the next vertex, which 
             will be on the face that is cut 4 times. */
          for(ie= 0; ie < nedg; ie++) {
            if( !(face_code[ the_polys[ipn].edges[ie] ] & last4) ) {
              /* this edge is not on the "special" face. Start the
                 tri-strip at the next vertex in the polygon */
              first_vert= (ie+1) % nedg;
              break;
            }
          }
          assemble_strip(first_vert, ipn, t_strips);
        } else {
          /* tri-strip can be started anywhere */
          assemble_strip(0, ipn, t_strips);
        }
        break;
      case 7:
      case 8:
        if(nfours > 0) {
          /* Find the first edge not on the face with four edges cut.
             If the next edge is on the face that is cut 4 times, start 
             the tri-strip at that next vertex. If not, skip one more vertex
             and start the tri-strip there (the start will now be on the 
             face with 4 edges cut). */
          for(ie= 0; ie < nedg; ie++) {
            if( !(face_code[ the_polys[ipn].edges[ie] ] & last4) ) {
              first_vert= (ie+1) % nedg;
              if( !(face_code[ the_polys[ipn].edges[ie+1] ] & last4) ) {
                first_vert= (ie+2) % nedg;
              }
              break;
            }
          }
          assemble_strip(first_vert, ipn, t_strips);
        } else {
          /* tri-strip can be started anywhere */
          assemble_strip(0, ipn, t_strips);
        }
        break;
      case 9:
        if(nfours > 0) {
          /* There are three faces that have all 4 edges cut by the iso-surface.
             The tri-strip should be started at the first edge that is shared 
             between two faces. Find that edge. */
          for(ie= 0; ie < nedg; ie++) {
            int nfaces= 0;
            for(fn= 0; fn < 6; fn++) {
              if(fnum[fn] >= 3) nfaces++;
            }
            if(nfaces > 1) {
              assemble_strip(ie, ipn, t_strips);
              break;
            }
          }
#undef BAD_NEDG_CHEK
#ifdef BAD_NEDG_CHEK
          /* should never be able to drop out of the loop without finding
             the desired edge */
          if(ie >= nedg) {
            int nbad= 1;
          }
#endif
        } else {
          /* impossible error */
          assemble_strip(0, ipn, t_strips);
        }
        break;
      case 12:
        if(nfours > 0) {
          /* There are six faces that have all 4 edges cut by the iso-surface.
             The tri-strip should be started at the first edge on face zero. */
          for(ie= 0; ie < nedg; ie++) {
            if(face_code[ the_polys[ipn].edges[ie] ] & 1) {
              /* have to make sure there isn't a problem with wrap-around */
              if(ie == 0 && !(face_code[ the_polys[ipn].edges[ie+1] ] & 1) ) {
                assemble_strip(nedg-1, ipn, t_strips);
                break;
              } else {
                assemble_strip(ie, ipn, t_strips);
                break;
              }
            }
          }
#ifdef BAD_NEDG_CHEK
          /* should never be able to drop out of the loop without finding
             the desired edge */
          if(ie >= nedg) {
            int nbad= 1;
          }
#endif
        } else {
          /* impossible error */
          assemble_strip(0, ipn, t_strips);
        }
        break;
      case 1:
      case 2:
      case 10:
      case 11:
        /* impossible!! */
#ifdef BAD_NEDG_CHEK
        {
          int badone= 1;
        }
#endif
        break;
    }
  }
#ifdef BAD_FACE_CHEK
  /* For all tri-strips, find cases with three consecutive edges on the same face.  */
  for(ipn= 0; ipn < num_poly; ipn++) {
    int i1, i2;
    nedg= t_strips[ipn].nvert;
    for(ie= 0; ie < nedg; ie++) {
      i1= (ie+1) % nedg;
      i2= (ie+2) % nedg;
      if(face_code[ t_strips[ipn].edges[ie] ] & 
         face_code[ t_strips[ipn].edges[i1] ] & 
         face_code[ t_strips[ipn].edges[i2] ] ) {
        int badone= 1;
      }
    }
  }
#endif
  return num_poly;
}

int ycTetIso_one_zone(double level, double *var, tstrip_array *t_strips)
{
  int ii, nstrip;

  /* Build the case table */
  for ( ii=0; ii < 8; ii++) {
    if ( var[ii] >= level ) {
      vertflag[ii]= 1;
    } else {
      vertflag[ii]= 0;
    }
  }
  nstrip= tetiso_zone( t_strips->the_tris );
  t_strips->nStrip= nstrip;
  return nstrip;
}

int ycPrepIsoTet(void)
{
  int ii, jj, kk, nstrip, mask;
  long *lens, *edges, lentot, now;

  /* release storage if necessary */
  if(have_iso_cases) {
    for(jj= 0; jj < 256; jj++) {
      if(iso_cases[jj].lens) p_free(iso_cases[jj].lens);
      if(iso_cases[jj].edges) p_free(iso_cases[jj].edges);
    }
    have_iso_cases= 0;
  }
  /* Run over all possible hi-lo combinations for the vertices
     of a hexahedron */
  for(jj= 0; jj < 256; jj++) {
    /* Build the case table */
    for ( ii=0, mask= 1; ii < 8; ii++, mask += mask) {
      if(jj & mask) {
        vertflag[ii]= 1;
      } else {
        vertflag[ii]= 0;
      }
    }
    nstrip= tetiso_zone( the_strips );
    iso_cases[jj].nStrip= nstrip;
    if(nstrip) {
      lens= (long *) p_malloc(nstrip*sizeof(long));
      iso_cases[jj].lens= lens;
      lentot= 0;
      for(ii= 0; ii < nstrip; ii++) {
        lens[ii]= the_strips[ii].nvert;
        lentot += lens[ii];
      }
      edges= (long *) p_malloc(lentot*sizeof(long));
      iso_cases[jj].edges= edges;
      now= 0;
      for(ii= 0; ii < nstrip; ii++) {
        for(kk= 0; kk < lens[ii]; kk++) {
          edges[now++]= the_strips[ii].edges[kk];
        }
      }
    } else {
      iso_cases[jj].lens= 0;
      iso_cases[jj].edges= 0;
    }
  }
  have_iso_cases= 1;
  return 0;
}

void assemble_strip(int ndx, int ipn, tstrip_data *t_strips)
{
  int phase= 0, lo, hi, nu, nedg, ie;

  /* start a triangle strip at index ndx */
  nedg= the_polys[ipn].num_edge;
  lo= ndx;
  hi= (lo-1);
  if(hi < 0) hi += nedg;
  nu= (lo+1);
  if(nu >= nedg) nu -= nedg;
  phase= 0;
  t_strips[ipn].edges[0]= the_polys[ipn].edges[lo];
  t_strips[ipn].edges[1]= the_polys[ipn].edges[hi];
  for(ie= 2; ie < nedg; ie++) {
    t_strips[ipn].edges[ie]= the_polys[ipn].edges[nu];
    if(phase == 0) {
      /* the new point in the strip is currently at the "low" side, 
         so it will be at the "high" side next time */
      lo= nu;
      nu= hi-1;
      if(nu < 0) nu += nedg;
      phase= 1;
    } else {
      /* the new point moves from the "high" side to the "low" side */
      hi= nu;
      nu= lo+1;
      if(nu >= nedg) nu -= nedg;
      phase= 0;
    }
  }
  t_strips[ipn].nvert= nedg;
}

void extract_tris_tet(int case_index, long idxcg, double *cntr_v2, long *numTri,
                      double lev, double s[8], yPoint3D pts[8],
                      yPoint3D gradients[8], double vars[8], long *cellIDs,
                      yPoint3D *xyzverts, yPoint3D *normals, double *v2vals)
{
  int indbase, ie, vert0, vert1, ij, ind;
  long ipn, nedg, numVert;
  double t, v1, v2;
  yPoint3D n, x1, x2, n1, n2;
  Edge the_edg;
  
  num_poly= iso_cases[case_index].nStrip;
  /* For all tri-strips, find the intersection points.
     Break them back up into triangles for now. */
  for(ipn= 0, indbase= 0; ipn < num_poly; ipn++) {
    int phase= 1;
    nedg= (iso_cases[case_index].lens)[ipn];
    for(ie= 0; ie < nedg-2; ie++) {
      for(ij= 0; ij < 3; ij++) {
        if(phase == 0) {
          ind= ie+ij;
        } else {
          ind= ie+2-ij;
        }
        /* get the two vertices at the ends of the edge on which the
           point lies */
        the_edg= edges[ (iso_cases[case_index].edges)[indbase+ind] ];
        vert0= the_edg.vert0;
        vert1= the_edg.vert1;
        t = (double)(lev - s[vert0]) / (s[vert1] - s[vert0]);
        x1 = pts[vert0];
        x2 = pts[vert1];
        n1 = gradients[vert0];
        n2 = gradients[vert1];
        numVert= 3*(*numTri)+ij;
        xyzverts[numVert].x = x1.x + t * (x2.x - x1.x);
        xyzverts[numVert].y = x1.y + t * (x2.y - x1.y);
        xyzverts[numVert].z = x1.z + t * (x2.z - x1.z);
        if(cntr_v2) {
          v1 = vars[vert0];
          v2 = vars[vert1];
          v2vals[numVert] = v1 + t * (v2-v1);
        }
        n.x = n1.x + t * (n2.x - n1.x);
        n.y = n1.y + t * (n2.y - n1.y);
        n.z = n1.z + t * (n2.z - n1.z);
        ycNormalize(&n);
        normals[numVert]= n;
      }
      cellIDs[*numTri]= idxcg;  /* the C-style global index of this cell */
      (*numTri)++;
      if(phase == 0) phase= 1;
      else phase= 0;
    }
    indbase += nedg;
  }
}

void extract_tris_tet_ndx(int case_index, long idxv, long idxcg, double *cntr_v2, long *numTri,
             long *numEdg, double lev, double s[8], yPoint3D pts[8],
             yPoint3D gradients[8], double vars[8], long edg_offset[12], long *cellIDs,
             long *ptndx, long *ndx, yPoint3D *xyzverts, yPoint3D *normals, double *v2vals)
{
  int indbase, ie, vert0, vert1, ij, ind;
  long ipn, nedg, numVert, itmp, edgnum, currEdg, currTri;
  double t, v1, v2;
  yPoint3D n, x1, x2, n1, n2;
  Edge the_edg;
  
  num_poly= iso_cases[case_index].nStrip;
  /* For all tri-strips, find the intersection points.
     Break them back up into triangles for now. */
  currEdg= *numEdg;
  currTri= *numTri;
  for(ipn= 0, indbase= 0; ipn < num_poly; ipn++) {
    int phase= 1;
    nedg= (iso_cases[case_index].lens)[ipn];
    for(ie= 0; ie < nedg-2; ie++) {
      for(ij= 0; ij < 3; ij++) {
        if(phase == 0) {
          ind= ie+ij;
        } else {
          ind= ie+2-ij;
        }
        itmp= (iso_cases[case_index].edges)[indbase+ind];
        edgnum= 3*idxv+edg_offset[itmp];
        if(ndx[edgnum] < 0) {
          /* get the two vertices at the ends of the edge on which the
             point lies */
          the_edg= edges[itmp];
          ndx[edgnum]= currEdg;
          vert0= the_edg.vert0;
          vert1= the_edg.vert1;
          t = (double)(lev - s[vert0]) / (s[vert1] - s[vert0]);
          x1 = pts[vert0];
          x2 = pts[vert1];
          xyzverts[currEdg].x = x1.x + t * (x2.x - x1.x);
          xyzverts[currEdg].y = x1.y + t * (x2.y - x1.y);
          xyzverts[currEdg].z = x1.z + t * (x2.z - x1.z);
          n1 = gradients[vert0];
          n2 = gradients[vert1];
          n.x = n1.x + t * (n2.x - n1.x);
          n.y = n1.y + t * (n2.y - n1.y);
          n.z = n1.z + t * (n2.z - n1.z);
          ycNormalize(&n);
          normals[currEdg]= n;
          if(cntr_v2) {
            v1 = vars[vert0];
            v2 = vars[vert1];
            v2vals[currEdg] = v1 + t * (v2-v1);
          }
          currEdg++;
        }
        /* Insert vertex into the triangle array.
           Note: vertex can be on any of the edges already cut. */
        numVert= 3*currTri+ij;
        ptndx[numVert]= ndx[edgnum];
      }
      cellIDs[currTri]= idxcg;  /* the C-style globalindex of this cell */
      currTri++;
      if(phase == 0) phase= 1;
      else phase= 0;
    }
    indbase += nedg;
  }
  *numEdg= currEdg;
  *numTri= currTri;
}

void extract_slicetris_tet(int case_index, long idxcg, double *cntr_v2, long *numTri,
                  double s[8], yPoint3D pts[8], double vars[8],
                  long *cellIDs, yPoint3D *xyzverts, double *v2vals)
{
  long numVert, ipn, indbase, nedg, ie, ij, ind, vert0, vert1;
  double t, var_1, var_2;
  yPoint3D x1, x2;
  Edge the_edg;
  
  num_poly= iso_cases[case_index].nStrip;
  /* For all tri-strips, find the intersection points.
     Break them back up into triangles for now. */
  for(ipn= 0, indbase= 0; ipn < num_poly; ipn++) {
    int phase= 1;
    nedg= (iso_cases[case_index].lens)[ipn];
    for(ie= 0; ie < nedg-2; ie++) {
      for(ij= 0; ij < 3; ij++) {
        if(phase == 0) {
          ind= ie+ij;
        } else {
          ind= ie+2-ij;
        }
        /* get the two vertices at the ends of the edge on which the
           point lies */
        the_edg= edges[ (iso_cases[case_index].edges)[indbase+ind] ];
        vert0= the_edg.vert0;
        vert1= the_edg.vert1;
        /* NOTE: if the program is functionally correctly, this
           statement can't result in a divide by zero. */
        t = (double)(0.0 - s[vert0]) / (s[vert1] - s[vert0]);
        x1 = pts[vert0];
        x2 = pts[vert1];
        numVert= 3*(*numTri)+ij;
        xyzverts[numVert].x = x1.x + t * (x2.x - x1.x);
        xyzverts[numVert].y = x1.y + t * (x2.y - x1.y);
        xyzverts[numVert].z = x1.z + t * (x2.z - x1.z);
        /* note: no normal because slicing plane */
        if(cntr_v2) {
          var_1= vars[vert0];
          var_2= vars[vert1];
          v2vals[numVert]= var_1 + t * (var_2-var_1);
        }
      }
      /* the C-style index of this cell in the global array */
      cellIDs[*numTri]= idxcg;  
      (*numTri)++;
    } /* for each triangle */
  }
}
