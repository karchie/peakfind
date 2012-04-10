/**
 * Copyright (c) 2012 Washington University
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "peakfind.h"

void display(const float *x, const int dims[3], int slice, float thresh) {
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
      fputc(fabs(x[i]) < thresh ? ' ' : 'X', stdout);
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

  image = calloc(dim, sizeof(float));
  roi = calloc(dim, sizeof(float));
  mask = malloc(dim * sizeof(float));

  if (!image || !mask || !roi) {
    fprintf(stderr, "allocation failed");
    exit(1);
  }

  /* start with a trivial mask */
  for (i = iz = 0; iz < nz; iz++) {
    for (iy = 0; iy < ny; iy++) {
      for (ix = 0; ix < nx; ix++, i++) {
        mask[i] = 1.0;
      }
    }
  }

  /* make a cube */  
  for (iz = 20; iz < 26; iz++) {
    for (iy = 4; iy < 8; iy++) {
      for (ix = 4; ix < 8; ix++) {
        i = ix + nx * (iy + ny * iz);

        image[i] = 1.0;
      }
    }
  }

  display(image, dims, 22, 0.5);

  sphereblur(image, dims, mmppixr, blur_radius);

  display(image, dims, 22, 0.5);

  find_peaks(image, dims, mmppixr, centerr,
             0, 0, 0, 0, 1,
             roi, 4, 8, 0, mask);

  display(roi, dims, 22, 0.5);
}
