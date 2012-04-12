/*
 * peakfind
 * MEX file for calling find_peaks and sphereblur from MATLAB
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
#define CENTER prhs[2]
#define ORAD   prhs[3]
#define MASK   prhs[4]

#define ROI    plhs[0]


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
  int dims[DIMS];
  float vtneg = 0, vtpos = 0, ctneg = 0, ctpos = 0, dthresh = 0;
  int minvox = 1;
  float mmppixr[DIMS], centerr[DIMS];

  /* REQUIRED: image, mmppixr, centerr, orad
   * OPTIONAL: mask, vtpos, vtneg, ctpos, ctneg, dthresh
   */
  if (nrhs < 5) {
    mexErrMsgIdAndTxt("MATLAB:peakfind:invalidNumInputs",
                      "arguments: image, mmppixr, centerr, orad, mask");
  }
      
  /* FIRST: peak cell array
   * SECOND: roi
   */
  if (nlhs > 2) {
    mexErrMsgIdAndTxt("MATLAB:peakfind:maxlhs",
                      "Too many output arguments.");
  }

  if (DIMS != mxGetNumberOfDimensions(IMAGE) || !mxIsSingle(IMAGE)) {
    mexErrMsgIdAndTxt("MATLAB:peakfind:invalidNumDimensions",
                      "PEAKFIND requires single precision 3D matrix");
  }

  mxCopyPtrToInteger4(mxGetDimensions(IMAGE), dims, DIMS);

  if (3 != mxGetN(MMPVOX) || !mxIsSingle(MMPVOX) ||
      3 != mxGetN(CENTER) || !mxIsSingle(CENTER)) {
    mexErrMsgIdAndTxt("MATLAB:peakfind:invalidSpaceParams",
                      "mmppixr, centerr must be single precision 3-vectors");
  }

  ROI = mxCreateNumericArray(DIMS, dims, mxSINGLE_CLASS, mxREAL);
  find_peaks((float*)mxGetData(IMAGE), dims,
             (float*)mxGetData(MMPVOX), (float*)mxGetData(CENTER),
             vtneg, vtpos, ctneg, ctpos, dthresh,
             (float*)mxGetData(ROI), *(float*)mxGetData(ORAD),
             minvox, (float*)mxGetData(MASK));
  return;
}
