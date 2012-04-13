/*
 * peakfind
 * Declarations for peak finder
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie <karchie@wustl.edu>
 *
 * emacs -*- indent-tabs-mode: nil; c-basic-offset: 4 -*-
 */

#ifndef _PEAKFIND_H_
#define _PEAKFIND_H_

#include <stdio.h>


typedef struct {
    float      x[3];            /**< location of peak */
    float      v;               /**<  value at peak */
    float      del2v;           /**<  Laplacian at peak */
    float      weight;
    int        nvox;
    int        killed;
} EXTREMUM;


extern void sphereblur(float *image, const int dim[3],
                       const float mmppixr[3], float radius);

extern void
find_peaks(float *image, const int dim[3],
           const float mmppixr[3], const float centerr[3],
           float vtneg, float vtpos, float ctneg, float ctpos,
           float dthresh,
           float *roi, float orad,
           int min_vox, const float *statmask,
	   int *npeaks, EXTREMUM **peaks);
#endif
