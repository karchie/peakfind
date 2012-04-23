/** emacs -*- indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * Utility declarations for peakfind
 *
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie <karchie@wustl.edu>
 */

#ifndef _PFUTIL_H_
#define _PFUTIL_H_

#include <stdarg.h>
#include <stdio.h>

/*
 * error handling
 */

extern void (*pf_error)(int code, const char *message);

#define PF_ERR_ALLOCATION 0x0001

extern void pf_error_handler_exit(int, const char *);
extern void pf_error_handler_warn(int, const char *);

extern void pf_set_error_handler(void (*)(int, const char *));

/*
 * logging
 */

extern void pf_log(int feature, char *format, ...);

#define PF_LOG_IMAGE_PADDING 0x01
#define PF_LOG_PEAK_PARAMS   0x02
#define PF_LOG_PEAK_TRACE    0x04
#define PF_LOG_PEAK_RESULTS  0x08
#define PF_LOG_UNDEF_POINT   0x10

extern void pf_log_to(FILE *fp);
extern void pf_log_to_stderr();
extern void pf_log_to_stdout();
extern void pf_log_image_padding(int should);
extern void pf_log_peak_params(int should);
extern void pf_log_peak_trace(int should);
extern void pf_log_peak_results(int should);
extern void pf_log_undef_point(int should);
extern void pf_log_all(int should);

#endif
