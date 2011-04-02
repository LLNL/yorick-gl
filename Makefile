# Makefile for yorick-gl, the yorgl OpenGL package
# $Id: Makefile,v 1.1 2005-09-18 22:07:55 dhmunro Exp $

Y_MAKEDIR=
Y_EXE=
Y_EXE_PKGS=
Y_EXE_HOME=
Y_EXE_SITE=
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

# ------------------------------------------------ macros for this package

PKG_NAME=yorgl
PKG_EXENAME=yorgl

YISO_PKG_I=cntrfunc.i
YGL_PKG_I=glfunc.i

YISO_I=contour.i tristruct.i
YGL_I=dlist3d.i slicenew.i glcompat.i glprofile.i testgl.i testisotree.i

YISO_OBJS=ContourTets3D.o Gradient3D.o isotree.o slicetree.o
YGL_OBJS=glPolys.o glStrips.o gltexture.o glcode.o glfunc.o glMouse.o \
  glx11view.o glx11setup.o TriUtil.o dlist3d.o glTarray.o gltexsubs.o \
  glustub.o glGlyph.o $(OGLXW)

include ./Makegl

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# -------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package

GLCODE=glcode.h dlist3d.h playgl.h
GLFUNC=glfunc.h $(GLCODE)
GLWRAP=glWrappers.h TriStruct.h
TriUtil.o: $(GLWRAP) $(GLFUNC)
dlist3d.o: glBasic.h glStrips.h $(GLWRAP) $(GLFUNC)
glGlyph.o: TriStruct.h $(GLFUNC)
glInfo.o: glInfo.h playgl.h
glMouse.o: glMouse.h $(GLWRAP) $(GLFUNC)
glPolys.o: glPolys.h $(GLFUNC)
glStrips.o: glStrips.h $(GLWRAP) $(GLFUNC)
glTarray.o: glStrips.h $(GLWRAP) $(GLFUNC)
glcode.o: $(GLFUNC)
glfunc.o: glBasic.h glMouse.h glStrips.h $(GLWRAP) $(GLFUNC)
gltexpal.o: Contour3D.h gl3dtex.h glMouse.h glcubetex.h glviewpoint.h $(GLFUNC)
gltexsubs.o: gl3dtex.h glBasic.h glcubetex.h $(GLFUNC)
gltexture.o: Contour3D.h gl3dtex.h glBasic.h glMouse.h glcubetex.h $(GLFUNC)
glustub.o: $(GLCODE)
glx11setup.o: glBasic.h glMouse.h $(GLFUNC)
glx11view.o: glBasic.h glMouse.h $(GLFUNC)
glxWinUtils.o: glxWinUtils.h
ContourTets3D.o: Contour3D.h
Gradient3D.o: Contour3D.h
isotree.o: Contour3D.h isotree.h
slicetree.o: Contour3D.h slicetree.h
oglx.o: playgl.h oglx.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(D_MESA_PIXMAPS) -o $@ -c oglx.c
oglw.o: playgl.h

distclean::
	./configure --distclean

# -------------------------------------------------------- end of Makefile
