/*
 * $Id: glInfo.h,v 1.1 2005-09-18 22:07:46 dhmunro Exp $
 */
/* Copyright (c) 2005, The Regents of the University of California.
 * All rights reserved.
 * This file is part of yorick (http://yorick.sourceforge.net).
 * Read the accompanying LICENSE file for details.
 */
typedef enum
{
   Normal,
   Wide,
   Verbose
} InfoMode;


struct visual_attribs
{
   /* X visual attribs */
   int id;
   int klass;
   int depth;
   int redMask, greenMask, blueMask;
   int colormapSize;
   int bitsPerRGB;

   /* GL visual attribs */
   int supportsGL;
   int transparent;
   int bufferSize;
   int level;
   int rgba;
   int doubleBuffer;
   int stereo;
   int auxBuffers;
   int redSize, greenSize, blueSize, alphaSize;
   int depthSize;
   int stencilSize;
   int accumRedSize, accumGreenSize, accumBlueSize, accumAlphaSize;
   int numSamples, numMultisample;
   int visualCaveat;
};
