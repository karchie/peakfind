/*
 * peakfind
 * MEX file for calling find_peaks from MATLAB
 *
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie
 */

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "mex.h"

#include "peakfind.h"
#include "pf_util.h"

#define DIMS 3

#define IMAGE  prhs[0]
#define MMPVOX prhs[1]
#define CENTER prhs[2]
#define ORAD   prhs[3]
#define MASK   prhs[4]

#define PEAKS  plhs[0]
#define ROI    plhs[1]

static const char *peak_fields[] = {
    "i", "ix", "x", "v", "del2v", "weight", "voxels"
};
static const int n_peak_fields = sizeof(peak_fields)/sizeof(*peak_fields);

#define INT_CLASS (4 == sizeof(int) ? mxINT32_CLASS : mxINT64_CLASS)

static const int KEYSIZE = 64;

float getFloatVal(const mxArray *a) {
    if (1 == mxGetM(a) && 1 == mxGetN(a)) {
        if (mxDOUBLE_CLASS == mxGetClassID(a)) {
            return (float)*(double*)mxGetData(a);
        } else if (mxSINGLE_CLASS == mxGetClassID(a)) {
            return *(float*)mxGetData(a);
        } else {
            mexErrMsgIdAndTxt("MATLAB:peakfind:invalidArgument",
                              mxGetClassName(a));
        }
    } else {
        mexErrMsgIdAndTxt("MATLAB:peakfind:invalidArgument",
                          "expected numeric value");
    }
}

int getIntVal(const mxArray *a) {
    return (int)floor(getFloatVal(a));
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    int i;
    const int *dims, *maskdims;
    float vtneg = 0, vtpos = 0, ctneg = 0, ctpos = 0, dthresh = 0, orad;
    int minvox = 1, ntop = 0;
    float *mmppixr, *centerr;
    int npeaks;
    EXTREMUM *peaks;
    FILE *logfp = 0;

    /* REQUIRED: image, mmppixr, centerr, orad, mask
     * OPTIONAL: vtpos, vtneg, ctpos, ctneg, dthresh
     */
    if (nrhs < 5) {
        mexErrMsgIdAndTxt("MATLAB:peakfind:invalidNumInputs",
                          "arguments: image, mmppixr, centerr, orad, mask");
    }
    
    for (i = 5; i < nrhs; i++) {
        char keybuf[KEYSIZE];
        mxArray *v;
        switch (mxGetString(prhs[i], keybuf, KEYSIZE-1)) {
        case 0:
            if (0 == strcmp("vtpos", keybuf)) {
                vtpos = getFloatVal(prhs[++i]);
            } else if (0 == strcmp("vtneg", keybuf)) {
                vtneg = getFloatVal(prhs[++i]);
            } else if (0 == strcmp("ctpos", keybuf)) {
                ctpos = getFloatVal(prhs[++i]);
            } else if (0 == strcmp("ctneg", keybuf)) {
                dthresh = getIntVal(prhs[++i]);
            } else if (0 == strcmp("minvox", keybuf)) {
                minvox = getIntVal(prhs[++i]);
            } else if (0 == strcmp("ntop", keybuf)) {
                ntop = getIntVal(prhs[++i]);
            } else if (0 == strcmp("log", keybuf)) {
                char logpath[PATH_MAX];
                int r = mxGetString(prhs[++i], logpath, PATH_MAX-1);
                if (0 == r) {
                    logfp = fopen(logpath, "w+");
                    if (logfp) {
                        pf_log_to(logfp);
                    } else {
                        mexErrMsgIdAndTxt("MATLAB:peakfind:invalidLogFile",
                                          strerror(errno));
                    }
                } else {
                    mexErrMsgIdAndTxt("MATLAB:peakfind:invalidLogFile",
                                      "unable to determine logfile path");
                }
            } else {
                mexErrMsgIdAndTxt("MATLAB:peakfind:invalidKey", keybuf);
            }
            break;
        case 1:
            mexErrMsgIdAndTxt("MATLAB:peakfind:invalidKey",
                              "key argument too long");
        default:
            mexErrMsgIdAndTxt("MATLAB:peakfind:invalidKey",
                              "valid keys: vtneg vtpos dthresh minvox");
        }
        
    }

    /* FIRST: peak cell array
     * SECOND: roi
     */
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:peakfind:maxlhs",
                          "[peaks roi] = peakfind(...)");
    }
    
    if (DIMS != mxGetNumberOfDimensions(IMAGE) || !mxIsSingle(IMAGE)) {
        mexErrMsgIdAndTxt("MATLAB:peakfind:invalidNumDimensions",
                          "image must be single precision 3D matrix");
    }
    dims = mxGetDimensions(IMAGE);
    
    if (DIMS != mxGetN(MMPVOX) || !mxIsSingle(MMPVOX) ||
        DIMS != mxGetN(CENTER) || !mxIsSingle(CENTER)) {
        mexErrMsgIdAndTxt("MATLAB:peakfind:invalidSpaceParams",
                          "mmppixr, centerr must be single precision 3-vectors");
    }
    
    if (DIMS != mxGetNumberOfDimensions(IMAGE) || !mxIsSingle(IMAGE)) {
        mexErrMsgIdAndTxt("MATLAB:peakfind:invalidNumDimensions",
                          "mask must be single precision 3D matrix");
    }
    
    maskdims = mxGetDimensions(MASK);
    for (i = 0; i < DIMS; i++) {
        if (dims[i] != maskdims[i]) {
            mexErrMsgIdAndTxt("MATLAB:peakfind:invalidMaskDimensions",
                              "mask dimensions must match image");
        }
    }

    if (mxIsSingle(ORAD)) {
        orad = *(float*)mxGetData(ORAD);
    } else if (mxIsDouble(ORAD)) {
        orad = (float)*(double*)mxGetData(ORAD);
    } else {
        mexErrMsgIdAndTxt("MATLAB:peakfind:invalidROIRadius",
                          "ROI radius must be floating point type");
    }

    ROI = mxCreateNumericArray(DIMS, dims, mxSINGLE_CLASS, mxREAL);
    find_peaks((float*)mxGetData(IMAGE), mxGetDimensions(IMAGE),
               mmppixr = (float*)mxGetData(MMPVOX),
               centerr = (float*)mxGetData(CENTER),
               vtneg, vtpos, ctneg, ctpos, dthresh,
               (float*)mxGetData(ROI), orad, ntop, minvox, 
               (float*)mxGetData(MASK), &npeaks, &peaks);

    PEAKS = mxCreateStructMatrix(1, npeaks, n_peak_fields, peak_fields);
    for (i = 0; i < npeaks; i++) {
        int j;
        EXTREMUM p = peaks[i];

        mxArray *ii = mxCreateNumericMatrix(1, DIMS, INT_CLASS, 0);
        int *ivals = (int *)mxGetData(ii);
        for (j = 0; j < DIMS; j++) {
            ivals[j] = p.i[j];
        }
        mxSetField(PEAKS, i, peak_fields[0], ii);

        mxArray *xi = mxCreateNumericMatrix(1, DIMS, mxSINGLE_CLASS, 0);
        float *ixvals = (float *)mxGetData(xi);
        for (j = 0; j < DIMS; j++) {
            ixvals[j] = (p.x[j] + centerr[j])/mmppixr[j];
        }
        mxSetField(PEAKS, i, peak_fields[1], xi);

        mxArray *x = mxCreateNumericMatrix(1, DIMS, mxSINGLE_CLASS, 0);
        float *xvals = (float *)mxGetData(x);
        for (j = 0; j < DIMS; j++) {
            xvals[j] = p.x[j];
        }
        mxSetField(PEAKS, i, peak_fields[2], x);

        mxArray *v = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, 0);
        *(float*)mxGetData(v) = p.v;
        mxSetField(PEAKS, i, peak_fields[3], v);

        mxArray *del2v = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, 0);
        *(float*)mxGetData(del2v) = p.del2v;
        mxSetField(PEAKS, i, peak_fields[4], del2v); 
       
        mxArray *w = mxCreateNumericMatrix(1, 1, INT_CLASS, 0);
        *(int*)mxGetData(w) = p.weight;
        mxSetField(PEAKS, i, peak_fields[5], w);

        mxArray *nvox = mxCreateNumericMatrix(1, 1, INT_CLASS, 0);
        *(int*)mxGetData(nvox) = p.nvox;
        mxSetField(PEAKS, i, peak_fields[6], nvox);
    }

    if (logfp) {
        pf_log_to_stdout();
        fclose(logfp);
    }
    free(peaks);
}
