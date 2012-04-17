/*
 * pf_set_log
 * MEX file for setting peakfind log output destination from MATLAB
 *
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie <karchie@wustl.edu>
 */

#include <stdio.h>

#include "matrix.h"
#include "mex.h"

#define PATH_MAX 1024

static FILE *log = 0;

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {
    switch (nrhs) {
    case 0: {
        if (log) {
            fclose(log);
            pf_log_to_stdout();
            log = 0;
        }
        return;
    }
            
    case 1: {
        char namebuf[PATH_MAX];
        int rval = mxGetString(prhs[0], namebuf, PATH_MAX-1);
        if (0 != rval) {
            mexErrMsgIdAndTxt("MATLAB:pf_set_log:invalidFilePath",
                              "file path argument could not be read");
        }
        if (log) {
            fclose(log);
        }
        log = fopen(namebuf, "w+");
        return;
    }

    default: mexErrMsgIdAndTxt("MATLAB:pf_set_log:invalidNumInputs",
                                   "at most one argument allowed");
    }
}
