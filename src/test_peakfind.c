/** emacs -*- indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * Copyright (c) 2012 Washington University
 * Author: Kevin A. Archie <karchie@wustl.edu>
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "peakfind.h"

float gauss3d(float dx, float dy, float dz, float s) {
  return expf(-(dx*dx + dy*dy + dz*dz)/(s*s));
}

static const char display_chars[] = {' ','.','+','*','#'};
static const int num_display_levels = sizeof(display_chars)/sizeof(char);

void display(const float *x, const int dims[3], int slice, float range) {
  int iz, iy, ix, i;

  iz = slice;
  i = dims[0] * dims[1] * iz;
  for (iy = 0; iy < dims[1]; iy++) {
    putc('v', stdout);
  }
  putc('\n', stdout);
  for (iy = 0; iy < dims[1] + 2; iy++) {
    putc('>', stdout);
    for (ix = 0; ix < dims[0]; ix++, i++) {
      int dci = (x[i]/range * num_display_levels);
      if (dci < 0) {
        putc('-',stdout);
      } else if (dci > num_display_levels) {
        putc('!',stdout);
      } else if (dci == num_display_levels) {
        putc(display_chars[num_display_levels-1],stdout);
      } else {
        putc(display_chars[dci],stdout);
      }
    }
    fputs("<\n", stdout);
  }
  for (iy = 0; iy < dims[1] + 2; iy++) {
    putc('^', stdout);
  }
  putc('\n', stdout);
}

int main(int argc, const char *argv[]) {
  const int nx = 32, ny = 32, nz = 32, dim = nx*ny*nz;
  const int dims[3] = {nx, ny, nz};
  const float mmppixr[] = { 2, 2, 2 }, centerr[] = { 0, 0, 0};
  const float blur_radius = 4.0;
  int iz, iy, ix, i;
  float *image, *mask, *roi;

  pf_log_all(1);

  image = calloc(dim, sizeof(float));
  roi = calloc(dim, sizeof(float));
  mask = malloc(dim * sizeof(float));

  if (!image || !mask || !roi) {
    fprintf(stderr, "allocation failed");
    exit(1);
  }

  /* start with a trivial mask and a gaussian image*/
  float min = INT_MAX, max = INT_MIN;
  for (i = iz = 0; iz < nz; iz++) {
    for (iy = 0; iy < ny; iy++) {
      for (ix = 0; ix < nx; ix++, i++) {
        mask[i] = 1.0;
        image[i] = gauss3d(ix-(nx/2),iy-(ny/2),iz-(ny/2),nx/4);
        if (image[i] > max) {
          max = image[i];
        }
        if (image[i] < min) {
          min = image[i];
        }
      }
    }
  }
  fprintf(stdout, "image value range: [%g, %g]\n", min, max);
    

  display(image, dims, 17, max);

  sphereblur(image, dims, mmppixr, blur_radius);

  display(image, dims, 17, max);

  find_peaks(image, dims, mmppixr, centerr,
             0, 0, 0, 0, 1,
             roi, 4, 0, 8, mask, 0, 0);

  puts("slice at peak:\n");
  display(roi, dims, 17, max);
  puts("2 slices away from peak:\n");
  display(roi, dims, 15, max);
  puts("4 slices away from peak:\n");
  display(roi, dims, 13, max);
}
