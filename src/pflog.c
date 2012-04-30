/** emacs -*- indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * Logging for peakfind
 *
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie <karchie@wustl.edu>
 */

#include <stdarg.h>
#include <stdio.h>

#include "pf_util.h"

static int verbosity = ~0;
static FILE *logfp = 0;

void pf_log(int option, char *message, ...) {
    va_list ap;
  
    va_start(ap, message);
    if (!logfp) {
      logfp = stdout;
    }
    vfprintf(logfp, message, ap);
    va_end(ap);
    fflush(logfp);
}

void pf_log_to(FILE *fp) {
  logfp = fp;
}

void pf_log_to_stderr() {
  logfp = stderr;
}

void pf_log_to_stdout() {
  logfp = stdout;
}

void pf_log_image_padding(int should) {
    if (should) {
        verbosity |= PF_LOG_IMAGE_PADDING;
    } else {
        verbosity &= ~PF_LOG_IMAGE_PADDING;
    }
}

void pf_log_peak_params(int should) {
    if (should) {
        verbosity |= PF_LOG_PEAK_PARAMS;
    } else {
        verbosity &= ~PF_LOG_PEAK_PARAMS;
    }
}

void pf_log_peak_trace(int should) {
    if (should) {
        verbosity |= PF_LOG_PEAK_TRACE;
    } else {
        verbosity &= ~PF_LOG_PEAK_TRACE;
    }
}

void pf_log_peak_results(int should) {
    if (should) {
        verbosity |= PF_LOG_PEAK_RESULTS;
    } else {
        verbosity &= ~PF_LOG_PEAK_RESULTS;
    }
}

void pf_log_undef_point(int should) {
    if (should) {
        verbosity |= PF_LOG_UNDEF_POINT;
    } else {
        verbosity &= ~PF_LOG_UNDEF_POINT;
    }
}

void pf_log_all(int should) {
    verbosity = should ? ~0 : 0;
}

