/** emacs -*- indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * Error handling for peakfind
 *
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie <karchie@wustl.edu>
 */

#include <stdlib.h>

#include "pf_util.h"


static char *error_messages[] = {
    0,
    "unable to allocate memory" /* PF_ERR_ALLOCATION */
};


/*
 * Error handling
 */

void pf_error_handler_warn(int type, const char *message) {
    fputs(error_messages[type], stderr);
    if (message) {
        fprintf(stderr, ": %s", message);
    }
    fputs("\n", stderr);
}    

void pf_error_handler_exit(int type, const char *message) {
    pf_error_handler_warn(type, message);
    exit(type);
}

void (*pf_error)(int, const char *) = pf_error_handler_exit;

void pf_set_error_handler(void (*f)(int, const char *)) {
    pf_error = f;
}
