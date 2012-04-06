/*
 * rms
 * Declarations for Avi's RMS Fortran subroutines
 * Copyright (c) 2012 Washington University
 */

#ifndef _RMS_H_
#define _RMS_H_

/**
 * How large should a padded vector be?
 * @param size data size
 * @param margin margin size
 * @return size of padded vector
 */
extern int npad_(int *size, int *margin);

/**
 * Center the 3D source data in a padded volume.
 * @param imag source data
 * @param nx dimension 1 of imag
 * @param ny dimension 2 of imag
 * @param nz dimension 3 of imag
 * @param imagp padded destination volume
 * @param dd0 dimension 1 of imagp
 * @param dd1 dimension 2 of imagp
 * @param dd2 dimension 3 of imagp
 */
extern void imgpad_(float *imag, int *nx, int *ny, int *nz,
		    float *imagp, int *dd0, int *dd1, int *dd2);

/**
 * Extract 3D data from a padded source volume, undoing imgpad_.
 * @param imag padded source data
 * @param nx dimension 1 of imag
 * @param ny dimension 2 of imag
 * @param nz dimension 3 of imag
 * @param imagp destination volume
 * @param dd0 dimension 1 of imagp
 * @param dd1 dimension 2 of imagp
 * @param dd2 dimension 3 of imagp
 */
extern void imgdap_(float *imag, int *nx, int *ny, int *nz,
		   float *imagp, int *dd0, int *dd1, int *dd2);

#endif
