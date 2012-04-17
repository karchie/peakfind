/*
 * sphereblur
 * MEX file for calling sphereblur from MATLAB
 *
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie
 */

#include "matrix.h"
#include "mex.h"

#include "peakfind.h"

#define DIMS 3

#define IMAGE  prhs[0]
#define MMPVOX prhs[1]
#define RADIUS prhs[2]

#define BLURRED plhs[0]

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    float radius;

    if (3 != nrhs) {
        mexErrMsgIdAndTxt("MATLAB:sphereblur:invalidNumInputs",
                          "arguments: image, mmppixr, radius");
    }

    if (DIMS != mxGetNumberOfDimensions(IMAGE) || !mxIsSingle(IMAGE)) {
        mexErrMsgIdAndTxt("MATLAB:sphereblur:invalidNumDimensions",
                          "image must be single precision 3D matrix");
    }
    if (DIMS != mxGetN(MMPVOX) || !mxIsSingle(MMPVOX)) {
        mexErrMsgIdAndTxt("MATLAB:sphereblur:invalidSpaceParams",
                          "mmppixr must be single precision 3-vector");
    }
    if (1 != mxGetN(RADIUS) || 1 != mxGetM(RADIUS)) {
        mexErrMsgIdAndTxt("MATLAB:sphereblur:invalidParameter",
                          "radius must be a floating point number");
    }
    if (mxIsSingle(RADIUS)) {
      radius = *(float*)mxGetData(RADIUS);
    } else if (mxIsDouble(RADIUS)) {
      radius = (float)*(double*)mxGetData(RADIUS);
    } else {
        mexErrMsgIdAndTxt("MATLAB:sphereblur:invalidParameter",
                          "radius must be a floating point number");
    }

    BLURRED = mxDuplicateArray(IMAGE);
    sphereblur((float*)mxGetData(BLURRED), mxGetDimensions(BLURRED),
               (float*)mxGetData(MMPVOX), radius);
}
