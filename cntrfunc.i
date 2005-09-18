/*
 * $Id: cntrfunc.i,v 1.1 2005-09-18 22:07:59 dhmunro Exp $
 * This file connects yorick to the compiled functions in 
 * the OpenGL based 3D graphics package.
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */

if (!is_void(plug_in)) plug_in, "yorgl";

/* set default chunk size for iso-surface and slice plane
   generation */
CHKSIZ= [20,16,10];

extern MakeSliceTreeCrv;
/* xxDOCUMENT  MakeSliceTreeCrv(xyz, tree)
    Build an octree for use in extracting slicing planes from an ijk grid.
    tree holds the size of the overall grid, the size of the chunk to process
    on this call, and the point where the chunk starts within the larger array.
    array. xyz is the coordinates of the grid points.
    The octree is returned in tree.
*/
/* PROTOTYPE
   int ycMakeSliceTreeCrv(double array xyz, pointer tree)
*/

extern SliceTreeCrv;
/* xxDOCUMENT  SliceTreeCrv(point, normal, xyz, var, triangles, tree)
    Compute a slicing plane on an ijk grid using a pre-computed octree.
    tree holds the size of the overall grid, the size of the chunk to process
    on this call, and the point where the chunk starts within the larger array.
    array. xyz is the coordinates of the grid points.
    var is an optional point centered variable to be interpolated to
    the vertices on the slicing plane.
    point lies on the slicing plane and normal is perpendicular
    to the plane.
    triangles is a struct with space in which to return the result.
    Makes a triangle array.
*/
/* PROTOTYPE
   int ycSliceTreeCrv(double array point, double array normal, 
                     double array xyz, pointer var, pointer triangles, 
                     pointer tree)
*/

extern SliceTree;
/* xxDOCUMENT  SliceTree(maxdepth, sizes, chunk, start, point, normal, 
                       deltas, origin, var2, triangles)
    Compute a slicing plane on a regular grid using bisection to
    efficiently find the cells cut by the plane. The number of
    bisections is maxdepth.
    The full grid has dimensions sizes=[nx,ny,nz]. This call handles
    chunk vertices and the chunk starts at location start in
    the full grid. 
    deltas is the size of a zone and origin is the coordinate of the
    first point.
    point lies on the slicing plane and normal is perpendicular
    to the plane.
    var2 (if not null) is an auxiliary variable whose value is returned
    at each vertex of the slice plane.
    triangles is a struct with space in which to return the result.
    Makes a triangle array.
*/
/* PROTOTYPE
   int ycSliceTree(long maxdepth, long array sizes, long array chunk, long array start, 
                   double array point,  double array normal, 
                   double array deltas, double array origin, 
                   pointer var2, pointer triangles)
*/

extern MakeContourTree;
/* xxDOCUMENT  MakeContourTree(var, tree)
    Build an octree for use in extracting iso-surfaces from a regular grid.
    tree holds the size of the overall grid, the size of the chunk to process
    on this call, and the point where the chunk starts within the larger array.
    array. var has cell centered values.
    The octree is returned in tree.
*/
/* PROTOTYPE
   int ycMakeContourTree(double array var, pointer tree)
*/

extern ContourTree;
/* xxDOCUMENT  ContourTree(deltas, origin, level, var, triangles, tree)
    Compute an iso-surface at the specified level on a regular grid
    using a pre-computed octree.
    tree holds the size of the overall grid, the size of the chunk to process
    on this call, and the point where the chunk starts within the larger array.
    array. var has cell centered values.
    level is the contour level and var is the array to contour.
    triangles is a struct with space in which to return the result.
    Makes a triangle array.
*/
/* PROTOTYPE
   int ycContourTree(double array deltas, double array origin, double level, 
                     double array var, pointer triangles, pointer tree)
*/

extern ContourTree2;
/* xxDOCUMENT  ContourTree2(deltas, origin, level, var, var2, triangles, tree)
    Compute an iso-surface at the specified level on a regular grid
    using a pre-computed octree.
    var2 is interpolated to the vertices of the iso-surface.
    tree holds the size of the overall grid, the size of the chunk to process
    on this call, and the point where the chunk starts within the larger array.
    array. var has cell centered values.
    level is the contour level and var is the array to contour.
    triangles is a struct with space in which to return the result.
    Makes a triangle array.
*/
/* PROTOTYPE
   int ycContourTree2(double array deltas, double array origin, double level, 
                     double array var, pointer var2, pointer triangles, pointer tree)
*/

extern ContourTreeCrv;
/* xxDOCUMENT  ContourTreeCrv(level, xyz, var, triangles, tree)
    Compute an iso-surface at the specified level on an ijk grid
    using a pre-computed octree.
    tree holds the size of the overall grid, the size of the chunk to process
    on this call, and the point where the chunk starts within the larger array.
    array. var has cell centered values.
    level is the contour level and var is the array to contour.
    triangles is a struct with space in which to return the result.
    Makes a triangle array.
*/
/* PROTOTYPE
   int ycContourTreeCrv(double level, double array xyz, 
                     double array var, pointer triangles, pointer tree)
*/

extern ContourTreeCrv2;
/* xxDOCUMENT  ContourTreeCrv2(level, xyz, var, var2, triangles, tree)
    Compute an iso-surface at the specified level on an ijk grid
    using a pre-computed octree.
    var2 is interpolated to the vertices of the iso-surface.
    tree holds the size of the overall grid, the size of the chunk to process
    on this call, and the point where the chunk starts within the larger array.
    array. var has cell centered values.
    level is the contour level and var is the array to contour.
    triangles is a struct with space in which to return the result.
    Makes a triangle array.
*/
/* PROTOTYPE
   int ycContourTreeCrv2(double level, double array xyz, double array var, 
             pointer var2, pointer triangles, pointer tree)
*/


extern ContourInitCartPcen;
/* xxDOCUMENT  ContourInitCartPcen(sizes, offsets, deltas, origin, var, var2)
    Prepare to compute an iso-surface at the specified level on a regular grid.
    The grid does NOT include any guard cells.
    var has point centered values.
    var2 is an auxiliary point centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first real cell in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth are the numbers of cells in the i, j, and k 
    directions in the global grid.
*/
/* PROTOTYPE
   int ycInitCartPcen(long array sizes, long array offsets, 
        double array deltas, double array origin, double array var, pointer var2)
 */

extern ContourInitCartGrdPcen;
/* xxDOCUMENT  ContourInitCartGrdPcen(sizes, offsets, deltas, origin, var, var2)
    Prepare to compute an iso-surface at the specified level on a regular grid.
    var has point centered values.
    var2 is an auxiliary point centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first real cell in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth are the numbers of cells in the i, j, and k 
    directions in the global grid.
*/
/* PROTOTYPE
   int ycInitCartGrdPcen(long array sizes, long array offsets, 
        double array deltas, double array origin, double array var, pointer var2)
 */

extern ContourInitCartZcen;
/* xxDOCUMENT  ContourInitCartZcen(offsets, deltas, origin, var, var2)
    Prepare to compute an iso-surface at the specified level on a regular grid.
    var has zone centered values.
    var2 is an auxiliary zone centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first real cell in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth arguemnts are the numbers of cells in 
    the i, j, and k directions in the global grid.
*/
/* PROTOTYPE
   int ycInitCartZcen(long array sizes, long array offsets, 
        double array deltas, double array origin, double array var, pointer var2)
 */

extern ContourInitCartGrdZcen;
/* xxDOCUMENT  ContourInitCartGrdZcen(offsets, deltas, origin, var, var2)
    Prepare to compute an iso-surface at the specified level on a regular grid.
    var has zone centered values.
    var2 is an auxiliary zone centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first real cell in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth arguments are the numbers of cells in 
    the i, j, and k directions in the global grid.
*/
/* PROTOTYPE
   int ycInitCartGrdZcen(long array sizes, long array offsets, 
        double array deltas, double array origin, double array var, pointer var2)
 */

extern ContourInitCartGrdPcenNdx;
/* xxDOCUMENT  ContourInitCartGrdPcenNdx(sizes, offsets, deltas, origin, var, var2)
    Prepare to compute an iso-surface at the specified level on a regular grid.
    var has point centered values.
    var2 is an auxiliary point centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first real cell in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth arguments are the numbers of cells in 
    the i, j, and k directions in the global grid.
    The output triangles are specified by indices into vectors of coordinates, 
    gradients, and auxiliary variables.
*/
/* PROTOTYPE
   int ycInitCartGrdPcenNdx(long array sizes, long array offsets, 
        double array deltas, double array origin, double array var, pointer var2)
 */

extern ContourInitCartGrdZcenNdx;
/* xxDOCUMENT  ContourInitCartGrdZcenNdx(sizes, offsets, deltas, origin, var, var2)
    Prepare to compute an iso-surface at the specified level on a regular grid.
    var has zone centered values.
    var2 is an auxiliary zone centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first real cell in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth arguments are the numbers of cells in 
    the i, j, and k directions in the global grid.
    xyz is the grid.
    The output triangles are specified by indices into vectors of coordinates, 
    gradients, and auxiliary variables.
*/
/* PROTOTYPE
   int ycInitCartGrdZcenNdx(long array sizes, long array offsets, 
        double array deltas, double array origin, double array var, pointer var2)
 */

extern ContourInitCrvGrdPcen;
/* xxDOCUMENT  ContourInitCrvGrdPcen(sizes, offsets, xyz, var, var2)
    Prepare to compute an iso-surface at the specified level on a curvi-linear grid.
    var has point centered values.
    var2 is an auxiliary point centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first point in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth arguments are the numbers of cells in 
    the i, j, and k directions in the global grid.
    xyz is the grid.
*/
/* PROTOTYPE
   int ycInitCrvGrdPcen(long array sizes, long array offsets, 
        double array xyz, double array var, pointer var2)
 */

extern ContourInitCrvGrdZcen;
/* xxDOCUMENT  ContourInitCrvGrdZcen(sizes, offsets, xyzn, var, var2)
    Prepare to compute an iso-surface at the specified level on a curvi-linear grid.
    var has zone centered values.
    var2 is an auxiliary zone centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first point in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth arguments are the numbers of cells in 
    the i, j, and k directions in the global grid.
    xyz is the grid.
*/
/* PROTOTYPE
   int ycInitCrvGrdZcen(long array sizes, long array offsets, 
        double array xyz, double array var, pointer var2)
 */

extern ContourInitCrvGrdPcenNdx;
/* xxDOCUMENT  ContourInitCrvGrdPcenNdx(sizes, offsets, xyz, var, var2)
    Prepare to compute an iso-surface at the specified level on a curvi-linear grid.
    var has point centered values.
    var2 is an auxiliary point centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first point in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth arguments are the numbers of cells in 
    the i, j, and k directions in the global grid.
    xyz is the grid.
    The output triangles are specified by indices into vectors of coordinates, 
    gradients, and auxiliary variables.
*/
/* PROTOTYPE
   int ycInitCrvGrdPcenNdx(long array sizes, long array offsets, 
        double array xyz, double array var, pointer var2)
 */

extern ContourInitCrvGrdZcenNdx;
/* xxDOCUMENT  ContourInitCrvGrdZcenNdx(sizes, offsets, xyz, var, var2)
    Prepare to compute an iso-surface at the specified level on a curvi-linear grid.
    var has zone centered values.
    var2 is an auxiliary zone centered variable that is the same size as var.
    sizes is the size of the chunk to be processed including the ghost
    points on both sides.
    offsets has 5 elements. The first three are the i, j, and k of
    the first point in this chunk of the grid relative to the global grid.
    They are 1-based, so they must be adjusted before use in C.
    The fourth, fifth, and sixth arguments are the numbers of cells in 
    the i, j, and k directions in the global grid.
    xyz is the grid.
    The output triangles are specified by indices into vectors of coordinates, 
    gradients, and auxiliary variables.
*/
/* PROTOTYPE
   int ycInitCrvGrdZcenNdx(long array sizes, long array offsets, 
        double array xyz, double array var, pointer var2)
 */


extern ContourTetHex;
/* xxDOCUMENT  ContourTetHex(level, ifirst, nzone, xyz, grad, hexndx, var, var2, triangles)
    compute an iso-surface at the specified level on an arbitrarily connected 
    grid of hexahedra.
    xyz, grad, and var are coordinates, gradients, and variable values 
    at the points making up the grid.
    var2 is an optional auxiliary variable that is the same size as var.
    hexndx has 8 indices into xyz etc. for each hex making up the grid.
    level is the contour level.
    triangles is a struct with space in which to return the result.
*/
/* PROTOTYPE
   int ycContourTetHex(double level, long ifirst, long nzone, double array xyz, 
        double array grad, long array hexndx, 
        double array var, pointer var2, pointer triangles)
 */

extern ContourTreeVarr;
/* xxDOCUMENT  ContourTreeVarr(deltas, origin, level, var, 
               triangles, tree)
    compute an iso-surface at the specified level on a regular grid.
    level is the contour level and var is the array to contour.
    var has cell centered values.
    triangles is a struct with space in which to return the result.
    tree contains a pre-computed octtree to be used in computing
    the iso-surface. The tree structure holds several pieces of 
    information about the grid etc.
    A triangle is specified by 3 indices into a list of 
    the unique vertices on the iso-surface which is stored in triangles.
*/
/* PROTOTYPE
   int ycContourTreeVarr(double array deltas, double array origin, 
         double level, double array var, pointer triangles, 
         pointer tree, long array edgndx)
 */

extern ContourTreeVarr2;
/* xxDOCUMENT  ContourTreeVarr2(deltas, origin, level, var, var2,
               triangles, tree)
    compute an iso-surface at the specified level on a regular grid.
    level is the contour level and var is the array to contour.
    var has cell centered values.
    triangles is a struct with space in which to return the result.
    tree contains a pre-computed octtree to be used in computing
    the iso-surface. The tree structure holds several pieces of 
    information about the grid etc.
    A triangle is specified by 3 indices into a list of 
    the unique vertices on the iso-surface which is stored in triangles.
    var2 is a second variable that will be interpolated to the vertices
    of the iso-surface (if var2 is non-zero).
*/
/* PROTOTYPE
   int ycContourTreeVarr2(double array deltas, double array origin, 
         double level, double array var, pointer var2,
         pointer triangles, pointer tree, long array edgndx)
 */

extern ContourTetArray;
/* xxDOCUMENT  ContourTetArray(make_strip, sizes, level, var, var2, 
               xyz, grid, flag, triangles)
    Compute an iso-surface at the specified level for a chunk within the
	full grid.
	A initialization routine suitable for the grid type must be called
    before this function for every chunk to establish the "mapping"
    between the full grid and the chunk.
    sizes is the number of points in each coordinate direction for the
    chunk passed in this call. There is a ghost point at both ends of each
    direction (included in sizes).
    level is the contour level and var is the array to contour.
    triangles is a struct with space in which to return the result.
    var2 is a second variable that will be interpolated to the vertices
    of the iso-surface (if var2 is non-zero).
*/
/* PROTOTYPE
   int ycContourTet_array(long make_strip, long array sizes, double level, 
              double array var, double array var2, double array xyz, 
              double array grd, char array flag, pointer triangles)
 */

extern ContourTetArrayNdx;
/* xxDOCUMENT  ContourTetArrayNdx(make_strip, sizes, level, var, var2, 
               xyz, grid, flag, ndx, triangles)
    Compute an iso-surface at the specified level for a chunk within the
    full grid.
    A initialization routine suitable for the grid type must be called
    before this function for every chunk to establish the "mapping"
    between the full grid and the chunk.
    sizes is the number of points in each coordinate direction for the
    chunk passed in this call. There is a ghost point at both ends of each
    direction (included in sizes).
    level is the contour level and var is the array to contour.
    triangles is a struct with space in which to return the result.
    A triangle is specified by 3 indices in ndx. They point into a list of 
    the unique vertices on the iso-surface which is stored in triangles.
    var2 is a second variable that will be interpolated to the vertices
    of the iso-surface (if var2 is non-zero).
*/
/* PROTOTYPE
   int ycContourTet_array_ndx(long make_strip, long array sizes, double level, 
              double array var, double array var2, double array xyz, 
              double array grd, char array flag, long array ndx, pointer triangles)
 */

extern ContourTetZone;
/* xxDOCUMENT  ContourTetZone(level, var, t_strips)
    Compute an iso-surface for a single hexahedral zone at the specified level.
    The result is the set of hex edges cut in the form of a (possibly
    more than one) tri-strip. The result is in the form of an array of
    IsoTriStrip.
    var is a 2-by-2-by-2 set of function values.
    level is the contour level.
*/
/* PROTOTYPE
   int ycTetIso_one_zone(double level, double array var, pointer t_strip)
 */

extern PrepIsoTet;
/* xxDOCUMENT  PrepIsoTet()
    Prepare the case table for iso-surfaces generated using a 6 tet 
    decompostion of a hexahedral zone.
*/
/* PROTOTYPE
   int ycPrepIsoTet(void)
 */
