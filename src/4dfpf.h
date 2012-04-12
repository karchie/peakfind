/*
 * 4dfpf
 * Declarations for Fortran subroutines from Avi's 4dfp code
 * Copyright (c) 2012 Washington University
 */

#ifndef _4DFPF_H_
#define _4DFPF_H_

/**
 * Blur provided image using a hard sphere kernel.
 * From gauss_4dfp
 * Uses FFT, REALT from fftsun.f (in librms)
 *
 * @param image image data to be blurred in place
 * @param nx dimension 1 of image (must be multiple of 2)
 * @param ny dimension 2 of image
 * @param nz dimension 3 of image
 * @param mmppix voxel dimensions in mm
 * @param radius blurring sphere radius in mm
 
 */
extern void hsphere3d_(float *image,
                       const int *nx, const int *ny, const int *nz,
		       const float *mmppix, const float *radius);

/**
 * Return into v the interpolated value of imgt at location x.
 * From imgreg_4dfp
 * 
 * @param imgt image data
 * @param nx dimension 1 of imgt
 * @param ny dimension 2 of imgt
 * @param nz dimension 3 of imgt
 * @param center x,y,z location of image center
 * @param mmppix voxel dimensions in mm
 * @param x location
 * @param v value at x
 * @param lslice slice from which most of the data are taken
 */
extern void imgvalx_(const float *imgt,
                     const int *nx, const int *ny, const int *nz,
		     const float *center, const float *mmppix,
		     const float *x, float *v, int *lslice);
#endif
