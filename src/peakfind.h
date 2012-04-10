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

extern void peakf_set_log_image_padding(int should);
extern void peakf_set_log_peak_params(int should);
extern void peakf_set_log_peak_trace(int should);
extern void peakf_set_log_peak_results(int should);
extern void peakf_set_log_undef_point(int should);
extern void peakf_set_log(int should);

extern void peakf_set_error_handler(void (*)(int, const char*));

extern void sphereblur(float *image, const int dim[3],
                       const float mmppixr[3], float radius);

extern void find_peaks(float *image, const int dim[3],
                       const float mmppixr[3], const float centerr[3],
                       float vtneg, float vtpos,
                       float ctneg, float ctpos,
                       float dthresh,
                       float *roi, float orad,
                       int min_vox, int polarize_roi,
                       const float *statmask);
#endif
