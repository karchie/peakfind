/** emacs -*- indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * sphereblur
 * Apply a spherical blur using Avi's FORTRAN routine hsphere3d.
 *
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie <karchie@wustl.edu>
 */

#include <stdlib.h>
#include "4dfpf.h"
#include "rms.h"
#include "pf_util.h"

/**
 * 
 * Convolves the provided image in place with a sphere of the provided
 * radius. Uses Avi's FORTRAN routine hsphere3d.
 *
 * @param image 3D image array
 * @param dimensions (x,y,z) of image
 * @param mmppixr mm/pixel for each dimension
 * @param radius radius of blurring sphere, in mm
 */
void sphereblur(float* image, const int dim[3],
                const float mmppixr[3], float radius) {
    int i;
    int pdim[3], psize = 1;
    int *vdim = (int*)dim;      /* really we don't change this */
    float *mmpvox = (float*)mmppixr;

    for (i = 0; i < 3; i++) {
        int margin = 2.0 * radius / mmppixr[i];
        pdim[i] = npad_((int*)dim + i, &margin);
        psize *= pdim[i];
    }
    pf_log(PF_LOG_IMAGE_PADDING,
           "image dimensions %d %d %d padded to %d %d %d\n",
           vdim[0], vdim[1], vdim[2], pdim[0], pdim[1], pdim[2]);
  
    float *padimage = calloc(psize, sizeof(float));
    if (!padimage) {
      pf_error(PF_ERR_ALLOCATION, "blur padding");
    }
  
    imgpad_(image, vdim, vdim+1, vdim+2, padimage, pdim, pdim+1, pdim+2);
    hsphere3d_(padimage, pdim, pdim+1, pdim+2, mmpvox, &radius);
    imgdap_(image, vdim, vdim+1, vdim+2, padimage, pdim, pdim+1, pdim+2);
}

