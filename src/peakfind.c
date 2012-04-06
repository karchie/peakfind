/**
 * peakfind.c
 * Find peaks in 3d slices of 4d images
 * Copyright (c) 2012 Washington University
 * Adapted from peak_4dfp by Avi Snyder and peak_nii by Tony Wilson
 * Author: Kevin A. Archie <karchie@wustl.edu>
 **/

#include <assert.h>
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

static const unsigned int NTOP = 1000;

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
			int ntop, int nvox_flag) {
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


/* value ordering for extrema */
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

/**
 * Returns the square of the distance between two peaks.
 *
 * @param p1 pointer to peak 1
 * @param p2 pointer to peak 2
 * @return square of distance between p1 and p2
 */
static float pdist2(EXTREMUM *p1, EXTREMUM *p2) {
  int k;
  float q = 0.0;
  
  for (k = 0; k < 3; k++) {
    float dk = p2->x[k] - p1->x[k];
    q += dk * dk;
  }
  return q;
}


/**
 * Combine the non-killed peaks from ppos and pneg into a single,
 * newly allocated array.
 *
 * @param ppall output array pointer
 * @param ppos positive peaks array
 * @param npos number of positive peaks
 * @param pneg negative peaks array
 * @param nneg number of negative peaks
 * @return number of combined peaks (i.e., size of *ppall)
 */
static int combine_extrema(float **ppall,
			   float *ppos, int npos,
			   float *pneg, int nneg) {
  int i, nall = 0;
  for (i = 0; i < npos; i++) {
    if (!ppos[i].killed) nall++;
  }
  for (i = 0; i < nneg; i++) {
    if (!pneg[i].killed) nall++;
  }

  if (nall > 0) {
    int iall = 0;

    *ppall = malloc(nall * sizeof(EXTREMUM));
    if (!*ppall) {
      error_handler(ERROR_ALLOCATION, "peak combination");
    }

    for (int i = 0; i < npos; i++) {
      if (!ppos[i].killed) {
	*ppall[iall++] = ppos[i];
      }
    }
    for (int i = 0; i < nneg; i++) {
      if (!pneg[i].killed) {
	*ppall[iall++] = ppos[i];
      }
    }
    assert(iall == nall);
  }
  return nall;
}
  

/**
 * Mark as killed any peaks that are too close to a strong peak, where
 * too close means within sqrt(d2thresh) mm. Move each surviving peak
 * to the center of mass of its cluster.
 *
 * @param plst array of peaks
 * @param nlst size of peak array
 * @param d2thresh square of threshold distance
 */
static void consolidate(EXTREMUM *plst, int nlst, float d2thresh) {
  int npair;
  do {
    float d2min = MAXFLOAT;
    int i, j, imin, jmin;

    /* On each pass, if any peaks are within threshold distance of
     * each other, consolidate just the two closest peaks.
     * Consolidation kills the larger-indexed peak, adds its weight to
     * the surviving peak, and moves the surviving peak to the
     * weighted average location.
     */

    /* Look for the closest within-threshold pair, if any */
    npair = 0;
    log(LOG_PEAK_RESULTS, "peak pairs closer than %.4f mm:\n",
	sqrt(d2thresh));
    for (i = 0; i < nlst; i++) {
      if (plst[i].killed) continue;
      for (j = i + 1; j < nlst; j++) {
	float t2;
	if (plst[j].killed) continue;
	t2 = pdist2(plst + i, plst + j);
	if (t2 < d2thresh) {
	  log(LOG_PEAK_RESULTS, "%5d%5d%10.4f\n", i + 1, j + 1, t2);
	  if (t2 < d2min) {
	    imin = i;
	    jmin = j;
	    d2min = t2;
	  }
	  npair++;
	}
      }
    }

    /* If any peak pairs were within threshold, consolidate the
     * closest pair.
     */
    log(LOG_PEAK_RESULTS, "npair = %d\n", npair);
    if (npair > 0) {
      int k;
      float x[3];

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
  } while (npair > 0);
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


/**
 * Fills roi with a peaks mask built from img.
 *
 * @param img original 3D image
 * @param roi output mask image: img where close to peak, 0 elsewhere
 * @param mask optional mask to apply to ROI (null if none)
 * @param dim dimensions of roi and (optional) mask
 * @param ppos array of local maxima
 * @param npos size of ppos
 * @param pneg array of local minima
 * @param nneg size of pneg
 * @param mmppixr mm/pixel for each dimension
 * @param centerr location offset, in mm
 * @param orad radius of peak spheres
 */
void build_mask(float *img, float *mask, float *roi, int dim[3],
		EXTREMUM pall[], int nall,
		float mmppixr[3], float centerr[3],
		float orad) {
  const float orad2 = orad * orad;
  EXTREMUM p;			/**< location (ix,iy,iz) */
  int i = 0, ix, iy, iz;
  
  // TODO: throw out local minima in positive territory,
  // local maxima in negative territory

  for (iz = 0; iz < dim[2]; iz++) {
    p.x[2] = (iz + 1) * mmppixr[2] - centerr[2];
    for (iy = 0; iy < dim[1]; iy++) {
      p.x[1] = (iy + 1) * mmppixr[1] - centerr[1];
      for (ix = 0; ix < dim[0]; ix++, i++) {
	float d2min = MAXFLOAT;
	int ip, ipmin, k;
	float tmin;

	p.x[0] = (ix + 1) * mmppixr[0] - centerr[0];

	/* find the closest peak to (ix,iy,iz) */
	for (ip = 0; ip < nall; ip++) {
	  float t2 = pdist2(&p, pall + ip);
	  if (t2 < d2min) {
	    ipmin = ip;
	    d2min = t2;
	  }
	}

	/* if the closest peak is within orad, mark this point. */
	/* TODO: check with Avi on mask values, 1e-37 comparison */
	if (d2min < orad2 && (!mask || fabs(mask[i]) > 1.0e-37))
	  roi[i] = img[i];
	  pall[ipmin].nvox++;
	} else {
	  roi[i] = 0.0;
	}
      }
    }
  }
}

/**
 * Find peaks in the provided 3d image.
 *
 * @param image 3D image data (as single-precision floats)
 * @param dim dimensions of image
 * @param mmppixr mm to pixel scaling factor
 * @param centerr location offset, in mm
 * @param vtneg
 * @param vtpos
 * @param ctneg
 * @param ctpos
 * @param dthresh peak distance threshold: if two peaks are
 *                within dthresh of each other, they are consolidated
 *                into a single, stronger peak
 */
void find_peaks(float *image, int dim[3],
		float mmppixr[3], float centerr[3],
		float vtneg, float vtpos, float ctneg, float ctpos,
		float dthresh) {
  EXTREMUM *ppos, *pneg, *pall;	/**< local maxima, minima, all extrema */
  int npos = 0, nneg = 0, nall;	/**< number of maxima, minima, extrema */
  int mpos = MSIZE, mneg = MSIZE; /**< array allocation sizes */
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

  int i = 0;
  for (iz = 0; iz < dim[2] - 1; iz++) {
    for (iy = 1; iy < dim[1] - 1; iy++) {
      for (ix = 1; ix < dim[0] - 1; ix++, i++) {
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
  consolidate(ppos, npos, d2thresh);
  qsort(pneg, nneg, sizeof(EXTREMUM), pcompare);
  consolidate(pneg, nneg, d2thresh);

  logpeaklist(ppos, npos, mmppixr, centerr, NTOP, 1);
  logpeaklist(pneg, nneg, mmppixr, centerr, NTOP, 1);

  nall = combine_extrema(&pall, ppos, npos, pneg, nneg);
  if (orad > 0) {
    build_mask(imgdata, null, roi, dim, pall, nall,
	       mmppixr, centerr, orad);
  }

  free(ppos);
  free(pneg);
  free(pall);
}
