/**
 * peakfind.c
 * Find peaks in 3d slices of 4d images
 * Copyright (c) 2012 Washington University
 * Adapted from peak_4dfp by Avi Snyder and peak_nii by Tony Wilson
 * Author: Kevin A. Archie <karchie@wustl.edu>
 **/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <varargs.h>

/*
 * external fortran routines:
 * npad_(int *size, int *padded_size)
 * imgpad_(float *src, int*, int*, int*, float *dst, int*, int*, int*)
 * imgdap_(float *dst, int*, int*, int*, float *src, int*, int*, int*)
 */

typedef struct
{
    float	x[3];		/**< location of peak */
    float	v;
    float	del2v;
    float	weight;
    int	nvox;
    int	killed;
} EXTREMUM;

static const unsigned int MSIZE = 2048;

/**
 * logging options
 */
static const unsigned int LOG_IMAGE_PADDING = 0x01;
static const unsigned int LOG_PEAK_PARAMS = 0x02;
static const unsigned int LOG_PEAK_TRACE = 0x04;
static const unsigned int LOG_PEAK_RESULTS = 0x08;
static const unsigned int LOG_UNDEF_POINT = 0x10;

/*
 * error types
 */
static const int ERROR_ALLOCATION = 0;

static char *error_messages[] = {
  "FATAL: memory allocation error",
};

/**
 * logging flags
 */
static int verbosity = ~0;

static void log(int option, char *message, ...) {
  va_list ap;
  
  va_start(ap, message);
  vfprintf(stdout, message, ap);
  va_end(ap);
  fflush(stdout);
}

/**
 * Set whether image padding operations should be logged.
 *
 * @param should nonzero if image padding should be logged,
 *        zero if it should not
 */
void set_log_image_padding(int should) {
  if (should) {
    verbosity |= LOG_IMAGE_PADDING;
  } else {
    verbosity &= ~LOG_IMAGE_PADDING;
  }
}


static void default_error(int type, const char *message) {
  fputs(error_messages[type], stderr);
  if (message) {
    fprintf(stderr, ": %s", message);
  }
  fputs("\n", stderr);
  exit(type);
}

static void (*error_handler)(int, const char *) = default_error;

void set_error_handler(void (*f)(int, const char *)) {
  error_handler = f;
}


static void logpeaklist(EXTREMUM *plst, int nlst,
			float mmppixr[3], float centerr[3],
			int nvox_flag) {
  int i, k;
  float fndex[3];
  
  log(LOG_PEAK_RESULTS, "%-5s%10s%10s%10s%10s%10s%10s%10s%10s%",
      "ROI", "index_x", "index_y", "index_z",
      "atlas_x", "atlas_y", "atlas_z", "value", "curvature");
  if (nvox_flag) {
    log(LOG_PEAK_RESULTS, "%10s", "nvox");
  }
  log(LOG_PEAK_RESULTS,"\n");
  for (i = 0; i < ntop && i < nlst; i++) {
    if (plst[i].killed) continue;
    for (k = 0; k < 3; k++) {
      fndex[k] = (plst[i].x[k] + centerr[k]) / mmppixr[k];
    }
    log(LOG_PEAK_RESULTS,
	"%-5d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.6f",
	ntot++ + 1, fndex[0], fndex[1], fndex[2],
	plst[i].x[0], plst[i].x[1], plst[i].x[2], plst[i].v,
	-plst[i].del2v);
    if (nvox_flag) {
      log(LOG_PEAK_RESULTS, "%10d", plst[i].nvox);
    }
    log(LOG_PEAK_RESULTS, "\n");
  }
}

/**
 * Convolves the provided image in place with a sphere of the provided
 * radius. Uses Avi's FORTRAN routine hsphere3d.
 *
 * @param image 3D image array
 * @param dimensions (x,y,z) of image
 * @param mmppixr mm/pixel for each dimension
 * @param radius radius of blurring sphere, in mm
 */
void sphereblur(float* image, int dim[3], double mmppixr[3],
		float radius) {
  int i;
  int pdim[3], psize = 1;

  for (i = 0; i < 3; i++) {
    int margin = 2.0 * radius / mmppixr[i];
    pdim[i] = npad_(dim + i, &margin);
    psize *= pdim[i];
  }
  log(LOG_IMAGE_PADDING,
      "image dimensions %d %d %d padded to %d %d %d\n",
      dim[0], dim[1], dim[2], pdim[0], pdim[1], pdim[2]);
  }
  
  float *padimage = calloc(psize, sizeof(float));
  if (!padimage) {
    error_handler(ERROR_ALLOCATION, "blur padding");
  }
  
  imgpad_(image, dim, dim+1, dim+2, padimage, pdim, pdim+1, pdim+2);
  hsphere3d_(padimage, pdim, pdim+1, pdim+2, mmppixr, &radius);
  imgdap_(image, dim, dim+1, dim+2, padimage, pdim, pdim+1, pdim+2);
}


static int pcompare(const void *ptr1, const void *ptr2) {
  EXTREMUM *p1 = ptr1, *p2 = ptr2;
  if (p1->v == p2->v) {
    return 0;
  } else if (fabs(p1-v) > fabs(p2->v)) {
    return -1;
  } else {
    return 1;
  }
}

static double pdist(EXTREMUM *p1, EXTREMUM *p2) {
  double q;
  int k;
  
  for (k = 0; k < 3; k++) {
    q += (p2->x[k] - p1->x[k]) * (p2->x[k] - p1->x[k]);
  }
  return sqrt(q);
}

static void consolidate(EXTREMUM *plst, int nlst) {
  float dmin = 1.0e6;
  do {
    float x[3], t;
    int i, j, imin, jmin, k;
    int npair = 0;
  
    log(LOG_PEAK_RESULTS, "peak pairs closer than %.4f mm:\n", dthresh);
    for (i = 0; i < ntop && i < nlst; i++) {
      if (plst[i].killed) continue;
      for (j = i + 1; j < ntop && j < nlst; j++) {
	if (plst[j].killed) continue;
	t = (float) pdist (plst + i, plst + j);
	if (t < dthresh) {
	  log(LOG_PEAK_RESULTS, "%5d%5d%10.4f\n", i + 1, j + 1, t);
	  if (t < dmin) {
	    imin = i;
	    jmin = j;
	    dmin = t;
	  }
	  npair++;
	}
      }
    }
    log(LOG_PEAK_RESULTS, "npair = %d\n", npair);
    if (npair) {
      for (k = 0; k < 3; k++) {
	x[k] = plst[imin].weight*plst[imin].x[k]
	  + plst[jmin].weight*plst[jmin].x[k];
      }
      plst[imin].weight += plst[jmin].weight;
      for (k = 0; k < 3; k++) {
	plst[imin].x[k] = x[k] / plst[imin].weight;
      }
      plst[jmin].killed++;
      npair--;
    }
  } while (npair);
}


/*! \def IS_PEAK(img, d, i, dir, op)
 * \brief Builds the test for whether this is a local extremum.
 *
 * \param img image data
 * \param d image dimensions (int[3])
 * \param i index of point being checked
 * \param dir pos or neg
 * \param op > (for positive) or < (for negative)
 *
 * Assumes the containing scope includes variables:
 *
 * del2v
 * ct{dir}
 */
#define IS_PEAK(img, d, i, dir, op) \
  img[i] op 0.0 && del2v op##= -ct##dir \
    && img[i] op img[i-1]         && img[i] op img[i+1] \
    && img[i] op img[i-d[0]]      && img[i] op img[i+d[0]] \
    && img[i] op img[i-d[0]*d[1]] && img[i] op img[i+d[0]*d[1]]

#define IS_POS_PEAK(img, d, i) IS_PEAK(img, d, i, pos, >)
#define IS_NEG_PEAK(img, d, i) IS_PEAK(img, d, i, neg, <)

/*! \def ADD_PEAK(img, d, i, dir, op)
 * \brief Builds the code for adding a local extremum.
 *
 * \param img image data
 * \param d image dimensions (int[3])
 * \param i index of point being checked
 * \param dir pos or neg
 * \param op > (for positive) or < (for negative)
 *
 * Assumes the containing scope includes variables:
 *
 * ix, iy, iz (x-, y-, z- of i, in pixel index space)
 * centerr, mmppixr, x, t, l
 * vt{dir}, m{dir}, n{dir}, p{dir}
 */
#define ADD_PEAK(img, d, i, dir, op) \
  imgvalx_(img, d, d+1, d+2, centerr, mmppixr, x, &t, &l); \
  if (l < 1) { \
    log(LOG_UNDEF_POINT, "undefined imgvalx point %d %d %d", ix, iy, iz); \
    continue; \
  } \
  if (vt##dir op t) continue; \
  if (m##dir <= n##dir) { \
    m##dir += MSIZE; \
    p##dir = realloc(p##dir, m##dir * sizeof(EXTREMUM)); \
    if (!p##dir)
      error_handler(ERROR_ALLOCATION, "extrema array update"); \
    for (k = 0; k < 3; k++) p##dir[n##dir].x[k] = x[k]; \
    p##dir[n##dir].v = t; \
    p##dir[n##dir].del2v = del2v; \
    p##dir[n##dir].weight = 1.0; \
    p##dir[n##dir].nvox = p##dir[n##dir].killed = 0; \
    n##dir++; \
  }

#define ADD_POS_PEAK(img, d, i) ADD_PEAK(img, d, i, pos, >)
#define ADD_NEG_PEAK(img, d, i) ADD_PEAK(img, d, i, neg, <)


void find_peaks(float *image, int dim[3],
		float mmppixr[3], float centerr[3],
		float vtneg, float vtpos, float ctneg, float ctpos) {
  EXTREMUM *ppos, *pneg;
  int mpos = MSIZE, mneg = MSIZE;
  int iz, iy, ix;

  ppos = malloc(mpos * sizeof(EXTREMUM));
  if (!ppos) {
    error_handler(ERROR_ALLOCATION, "positive extremum array");
  }
  pneg = malloc(mpos * sizeof(EXTREMUM));
  if (!pneg) {
    error_handler(ERROR_ALLOCATION, "negative extremum array");
  }

  log(LOG_PEAK_PARAMS,
      "peak value     thresholds %10.4f to %10.4f\n", vtneg, vtpos);
  log(LOG_PEAK_PARAMS,
      "peak curvature thresholds %10.4f to %10.4f\n", ctneg, ctpos);
  log(LOG_PEAK_TRACE, "compiling extrema slice");

  int nx = dim[0], nxy = dim[0]*dim[1];
  for (iz = 1; iz < dim[2] - 1; iz++) {
    for (iy = 1; iy < dim[1] - 1; iy++) {
      for (iz = 1; iz < dim[0] - 1; iz++) {
	int i = ix + dim[0] * (iy + dim[1] * iz);
	float del2v = 0.0, dvdx[3], d2vdx2[3], w[3], fndex[3];
	
	dvdx[0] = (-image[i - 1] + image[i + 1]) / 2.0;
	dvdx[1] = (-image[i - nx] + image[i + nx]) / 2.0;
	dvdx[2] = (-image[i - nxy] + image[i + nxy]) / 2.0;

	d2vdx2[0] = -2.0 * image[i] + image[i - 1] + image[i + 1];
	d2vdx2[1] = -2.0 * image[i] + image[i - nx] + image[i + nx];
	d2vdx2[2] = -2.0 * image[i] + image[i - nxy] + image[i + nxy];

	for (k = 0; k < 3; k++) {
	  w[k] = -dvdx[k] / d2vdx2[k];
	  del2v += d2vdx2[k] / (mmppixr[k] * mmppixr[k]);
	}

	fndex[0] = w[0] + ix + 1;
	fndex[1] = w[1] + iy + 1;
	fndex[2] = w[2] + iz + 1;
	for (k = 0; k < 3; k++) {
	  x[k] = fndex[k] * mmppixr[k] - centerr[k];
	}

	if (IS_POS_PEAK(image, dim, i)) {
	  MAKE_POS_PEAK;
	} else if (IS_NEG_PEAK(image, dim, i)) {
	  MAKE_NEG_PEAK;
	}
      }
    }
  }

  log(LOG_PEAK_TRACE, "\nbefore sorting npos = %d nneg = %d\n",
      npos, nneg);

  /* Sort by value, then consolidate nearby extrema */
  qsort(ppos, npos, sizeof(EXTREMUM), pcompare);
  qsort(pneg, nneg, sizeof(EXTREMUM), pcompare);

  consolidate(ppos, npos);
  consolidate(pneg, nneg);
  int npos0 = npos, nneg0 = nneg;

  logpeaklist(ppos, npos, mmppixr, centerr, 1);
  logpeaklist(pneg, nneg, mmppixr, centerr, 1);

  /* TODO: combine positive and negative peak lists
	   create output ROI? */
  
}
