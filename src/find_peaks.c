/** emacs -*- indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * peakfind.c
 * Find peaks in 3d volumes
 * Copyright (c) 2012 Washington University
 * Adapted from peak_4dfp by Avi Snyder and peak_nii by Tony Wilson
 * Author: Kevin A. Archie <karchie@wustl.edu>
 **/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "peakfind.h"
#include "pf_util.h"
#include "4dfpf.h"


static const unsigned int MSIZE = 2048;

static const unsigned int NTOP = 1000;

#define MIN(a,b) (a<b?a:b)

/**
 * qsort-compatible comparison function for sorting EXTREMUM structs
 * into descending value order
 * @param ptr1 one EXTREMUM
 * @param ptr2 another EXTREMUM
 * @return  0 if p1.v == p2.v,
 *         -1 if |p1.v| > |p2.v|,
 *          1 if |p1.v| < |p2.v|
 */
static int pcompare(const void *ptr1, const void *ptr2) {
    const EXTREMUM *p1 = ptr1, *p2 = ptr2;
    if (p1->v == p2->v) {
        return 0;
    } else if (fabs(p1->v) > fabs(p2->v)) {
        return -1;
    } else {
        return 1;
    }
}

struct peak_edge_t {
    int i1, i2;
    EXTREMUM *p1, *p2;
    float d2;
};

/**
 * qsort-compatible comparison function for sorting peak edges
 * into increasing distance order
 * @param pe1 one peak_edge_t
 * @param pe2 another peak_edge_t
 * @return  0 if pe1.d2 == pe2.d2
 *         -1 if pe1.d2 < pe2.d2
 *          1 if pe1.d2 > pe2.d2
 */
static int edge_compare(const void *pe1, const void *pe2) {
    const struct peak_edge_t *p1 = pe1, *p2 = pe2;
    if (p1->d2 == p2->d2) {
        return 0;
    } else {
        return p1->d2 < p2->d2 ? -1 : 1;
    }
}

/**
 * Returns the square of the distance between two peaks.
 *
 * @param p1 pointer to peak 1
 * @param p2 pointer to peak 2
 * @return square of distance between p1 and p2
 */
static float pdist2(const EXTREMUM *p1, const EXTREMUM *p2) {
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
 * @param ntop maximum number of peaks
 * @return number of combined peaks (i.e., size of *ppall)
 */
static int combine_extrema(EXTREMUM **ppall,
                           EXTREMUM *ppos, int npos,
                           EXTREMUM *pneg, int nneg,
                           int ntop) {
    int i, nall = 0;
    int np = MIN(npos,ntop), nn = MIN(nneg,ntop);

    pf_log(PF_LOG_PEAK_TRACE, "consolidating %d + %d peaks to 2*%d\n",
           npos, nneg, ntop);
    if (np > 0) {
        pf_log(PF_LOG_PEAK_TRACE, "first + peak value: %g\n", ppos[0].v);
    }
    if (nn > 0) {
        pf_log(PF_LOG_PEAK_TRACE, "first - peak value: %g\n", pneg[0].v);
    }

    for (i = 0; i < np; i++) {
        if (!ppos[i].killed) nall++;
    }
    for (i = 0; i < nn; i++) {
        if (!pneg[i].killed) nall++;
    }

    if (nall > 0) {
        int iall = 0;

        *ppall = malloc(nall * sizeof(EXTREMUM));
        if (!*ppall) {
            pf_error(PF_ERR_ALLOCATION, "peak combination");
        }

        for (i = 0; i < np; i++) {
            if (!ppos[i].killed) {
                (*ppall)[iall++] = ppos[i];
            }
        }
        for (i = 0; i < nn; i++) {
            if (!pneg[i].killed) {
                (*ppall)[iall++] = pneg[i];
            }
        }
        assert(iall == nall);
    } else {
        *ppall = 0;
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
 * @param ntop maximum number of peaks to consider
 * @param d2thresh square of threshold distance
 */
static void consolidate(EXTREMUM *plst, int nlst, int ntop, float d2thresh) {
    struct peak_edge_t *edges, *edge;
    int i, j, ip, nedges;
    float last_d2 = 0;
    int n = MIN(nlst,ntop);

    edges = malloc(n*n*sizeof(struct peak_edge_t));
    if (!edges) {
        pf_error(PF_ERR_ALLOCATION, "peak edges");
    }

    /* build the full set of edges */
    for (i = 0, nedges = 0; i < n; i++) {
        if (plst[i].killed) continue;
        for (j = i + 1; j < n; j++) {
            if (plst[j].killed) continue;
            int k = nedges++;
            edges[k].i1 = i;
            edges[k].i2 = j;
            edges[k].p1 = plst + i;
            edges[k].p2 = plst + j;
            edges[k].d2 = pdist2(plst + i, plst + j);
        }
    }

    pf_log(PF_LOG_PEAK_TRACE, "%d raw edges, sorting by distance...", nedges);

    qsort(edges, nedges, sizeof(struct peak_edge_t), edge_compare);

    pf_log(PF_LOG_PEAK_TRACE, "done.\n");

    for (i = 0; i < nedges && edges[i].d2 < d2thresh; i++) {
        int k;
        struct peak_edge_t *edge = edges + i;
        if (edge->p1 == edge->p2) {
            /* circular edges are former duplicates; skip */
            continue;
        }
        assert(!edge->p1->killed && !edge->p2->killed);

        float d2 = pdist2(edge->p1, edge->p2);
        if (d2 != edge->d2) {
            pf_log(PF_LOG_PEAK_TRACE, "!!! distance inconsistency: %g != %g\n",
                   d2, edge->d2);
        }
        if (last_d2 > edge->d2) {
            pf_log(PF_LOG_PEAK_RESULTS,
                   "!!! edge %i breaks distance ordering: %g > %g\n",
                   sqrt(last_d2), sqrt(edge->d2));
        }
        pf_log(PF_LOG_PEAK_TRACE,
               "merge %d <- %d (distance %g)\n",
               edge->i1+1, edge->i2+1, sqrt(edge->d2));
        
        if (edge->i1 > edge->i2) { /* take peaks in v order */
            int ti = edge->i2; edge->i2 = edge->i1; edge->i1 = ti;
            EXTREMUM *p = edge->p2; edge->p2 = edge->p1; edge->p1 = p;
        }
        edge->p2->killed = 1;
        float weight = edge->p1->weight + edge->p2->weight;
        for (k = 0; k < 3; k++) {
            edge->p1->x[k] = (edge->p1->weight * edge->p1->x[k]
                              + edge->p2->weight * edge->p2->x[k])
                / weight;
        }
        edge->p1->weight = weight;
        
        /* We're done with edge i now. p2 is killed and p1 has
           moved. Update following edges to fix distances for edges
           involving p1 or p2. This can create some duplicate
           edges. When the first of a set of duplicates is
           processed, the others will be converted to circular
           edges, bubbled to the front, and ignored. */
        for (j = i + 1; j < nedges; j++) {
            if (edges[j].p1 == edge->p2) {
                edges[j].i1 = edge->i1;
                edges[j].p1 = edge->p1;
                edges[j].d2 = pdist2(edges[j].p1,edges[j].p2);
            } else if (edges[j].p2 == edge->p2) {
                edges[j].i2 = edge->i1;
                edges[j].p2 = edge->p1;
                edges[j].d2 = pdist2(edges[j].p1,edges[j].p2);
            } else if (edges[j].p1 == edge->p1
                       || edges[j].p2 == edge->p1) {
                edges[j].d2 = pdist2(edges[j].p1,edges[j].p2);
            }
            
            /* edge j might now be closer than preceding edges; if
               it's within threshold, bubble it back as necessary to
               maintain distance ordering -- possibly all the way
               back to i+1. */
            int jbubb = j;
            for (k = j - 1; k > i && edges[jbubb].d2 < edges[k].d2; k--) {
                struct peak_edge_t temp = edges[jbubb];
                edges[jbubb] = edges[k];
                edges[k] = temp;
                jbubb = k;
            }
        }
    }

    free(edges);
}

#define UNMASKED(mask, i) (!mask || (fabs(mask[i]) > 1.0e-37))

/**
 * Fills roi with a peaks mask built from img.
 *
 * @param img original 3D image
 * @param mask statistical significance mask for constraining valid peaks
 * @param roi output mask image: img where close to peak, 0 elsewhere
 * @param dim dimensions of roi and (optional) mask
 * @param ppos array of local maxima
 * @param npos size of ppos
 * @param pneg array of local minima
 * @param nneg size of pneg
 * @param mmppixr mm/pixel for each dimension
 * @param centerr location offset, in mm
 * @param orad radius of peak spheres
 */
static void build_mask(const float *img, const float *mask,
                       float *roi, const int dim[3],
                       const float mmppixr[3], const float centerr[3],
                       const EXTREMUM *pall, int nall,
                       float orad) {
    const float orad2 = orad * orad;
    EXTREMUM loc;             /**< location (ix,iy,iz) in image space */
    int ix, iy, iz, i;

    if (0 == nall) {
        pf_log(PF_LOG_PEAK_RESULTS,
               "no peaks, skipping ROI mask build\n");
        return;
    } else {
        pf_log(PF_LOG_PEAK_TRACE, "starting ROI mask build\n");
    }

    i = 0;
    for (iz = 0; iz < dim[2]; iz++) {
        loc.x[2] = (iz + 1) * mmppixr[2] - centerr[2];
        for (iy = 0; iy < dim[1]; iy++) {
            loc.x[1] = (iy + 1) * mmppixr[1] - centerr[1];
            for (ix = 0; ix < dim[0]; ix++, i++) {
                float d2min = HUGE_VALF;
                int ip, ipmin = -1;

                loc.x[0] = (ix + 1) * mmppixr[0] - centerr[0];

                /* find the closest peak to (ix,iy,iz) */
                for (ip = 0; ip < nall; ip++) {
                    float d2 = pdist2(&loc, pall + ip);
                    if (d2 < d2min) {
                        ipmin = ip;
                        d2min = d2;
                    }
                }

                assert(ipmin >= 0);

                /* if the closest peak is within orad, mark this point. */
                if (d2min < orad2 && UNMASKED(mask, i)) {
                    roi[i] = ipmin + 2;
                } else {
                    roi[i] = 0.0;
                }
            }
        }
    }
}


/**
 * Counts pixels in prospective peak ROIs and cull ones that
 * are too small. Edits the provided peak list if peaks are culled.
 *
 * @param img source image
 * @param mask statistical significance mask
 * @param dim dimensions of source image (and statistical mask)
 * @param pallp pointer to array of peaks
 * @param nallp pointer to int size of pallp
 * @param mmppixr voxel dimensions in mm
 * @param centerr center coordinate location in mm
 * @param max_peaks largest number of peaks allowed
 * @param orad radius of ROI spheres
 * @param min_vox minimum voxel count for retained ROIs
 * @return number of culled peaks
 */
static int cull_small_extrema(const float *img, const float *mask,
                              const int dim[3],
                              const float mmppixr[3],
                              const float centerr[3],
                              EXTREMUM **pallp, int *nallp,
                              int max_peaks,
                              int orad, int min_vox) {
    const int orad2 = orad * orad;
    int iz, iy, ix, i, nkilled = 0, nall = MIN(*nallp,max_peaks);
    EXTREMUM loc;


    if (min_vox <= 0) {
        return 0;
    }
    if (0 == *nallp){
        pf_log(PF_LOG_PEAK_RESULTS,
               "no peaks found, skipping small peak cull\n");
        return 0;
    } else {
        pf_log(PF_LOG_PEAK_TRACE, "starting small extrema removal\n");
    }

    /* count voxels in each ROI */
    i = 0;
    for (iz = 0; iz < dim[2]; iz++) {
        loc.x[2] = (iz + 1)*mmppixr[2] - centerr[2];
        for (iy = 0; iy < dim[1]; iy++) {
            loc.x[1] = (iy + 1)*mmppixr[1] - centerr[1];
            for (ix = 0; ix < dim[0]; ix++, i++) {
                int ip, ipmin = -1;
                float d2min = HUGE_VALF;

                loc.x[0] = (ix + 1)*mmppixr[0] - centerr[0];
                for (ip = 0; ip < nall; ip++) {
                    float d2;
                    d2 = pdist2(&loc, *pallp + ip);
                    if (d2 < d2min) {
                        ipmin = ip;
                        d2min = d2;
                    }
                }
                assert(ipmin >= 0);

                if (d2min < orad2 && UNMASKED(mask, i)) {
                    (*pallp)[ipmin].nvox++;
                }
            }
        }
    }

    /* cull peaks that are too small by voxel count */
    for (i = 0; i < nall; i++) {
        if ((*pallp)[i].nvox < min_vox) {
            (*pallp)[i].killed = 1;
            nkilled++;
        }
    }
    /* also cull peaks that are past the peak count cutoff */
    for (; i < nall; i++) {
        (*pallp)[i].killed = 1;
        nkilled++;
    }

    /* if any peaks were killed, make a new, compacted peaks list */
    if (nkilled >= nall) {
        *pallp = 0;
        free(*pallp);
        *nallp = 0;
    } else if (nkilled > 0) {
        int j;
        EXTREMUM *pnew = malloc((nall - nkilled) * sizeof(EXTREMUM));
        if (!pnew) {
            pf_error(PF_ERR_ALLOCATION, "surviving peaks");
        }
        for (i = j = 0; i < nall; i++) {
            if (!(*pallp)[i].killed) {
                pnew[j++] = (*pallp)[i];
            }
        }
        assert(j == nall - nkilled);
        *nallp = j;
        *pallp = pnew;
    }
    return nkilled;
}


/*
 * Exported functions
 */



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
#define IS_PEAK(img, d, i, dir, op)                                     \
    img[i] op 0.0 && -ct ## dir op ## = del2v                           \
        && img[i] op img[i-1]         && img[i] op img[i+1]             \
        && img[i] op img[i-d[0]]      && img[i] op img[i+d[0]]          \
        && img[i] op img[i-d[0]*d[1]] && img[i] op img[i+d[0]*d[1]]

#define IS_POS_PEAK(img, d, i) IS_PEAK(img, d, i, pos, >)
#define IS_NEG_PEAK(img, d, i) IS_PEAK(img, d, i, neg, <)

/*! \def ADD_PEAK(img, d, i, dir, op)
 * \brief Builds the code for adding a local extremum.
 *
 * \param img image data
 * \param d image dimensions (int[3])
 * \param i index of point being checked
 * \param x
 * \param dir pos or neg
 * \param op > (for positive) or < (for negative)
 *
 * Assumes the containing scope includes variables:
 *
 * ix, iy, iz (x-, y-, z- of i, in pixel index space)
 * centerr, mmppixr, x, t
 * vt{dir}, m{dir}, n{dir}, p{dir}
 */
#define ADD_PEAK(img, d, i, dir, op)                                    \
    do {                                                                \
        float vx; /**< value at point x */                              \
        int slice; /**< slice most determining vx */                    \
        imgvalx_(img, (int*)d, (int*)d+1, (int*)d+2,                    \
                 (float*)centerr, (float*)mmppixr, x, &vx, &slice);     \
        if (slice < 1) {                                                \
            pf_log(PF_LOG_UNDEF_POINT,                                     \
                   "undefined imgvalx point %d %d %d\n", ix, iy, iz);     \
            continue;                                                   \
        }                                                               \
        if (vt ## dir op vx) continue;                                  \
        if (m ## dir <= n ## dir) {                                     \
            m ## dir += MSIZE;                                          \
            p ## dir = realloc(p ## dir, m ## dir * sizeof(EXTREMUM));  \
            if (!p ## dir)                                              \
                pf_error(PF_ERR_ALLOCATION, "extrema array update"); \
        }                                                               \
        for (k = 0; k < 3; k++) p ## dir[n ## dir].x[k] = x[k];         \
        p ## dir[n ## dir].i[0] = ix;                                   \
        p ## dir[n ## dir].i[1] = iy;                                   \
        p ## dir[n ## dir].i[2] = iz;                                   \
        p ## dir[n ## dir].v = vx;                                      \
        p ## dir[n ## dir].del2v = del2v;                               \
        p ## dir[n ## dir].weight = 1.0;                                \
        p ## dir[n ## dir].nvox = p ## dir[n ## dir].killed = 0;        \
        n ## dir++;                                                     \
    } while (0)

#define ADD_POS_PEAK(img, d, i) ADD_PEAK(img, d, i, pos, >)
#define ADD_NEG_PEAK(img, d, i) ADD_PEAK(img, d, i, neg, <)



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
 * @param roi 3D image space for output peaks mask
 * @param orad radius for peak spheres in roi
 * @param min_vox minimum voxel count for peak ROIs
 * @param statmask statistical significance mask
 * @param peaks peaks_t for returning peaks information
 */
void
find_peaks(float *image, const int dim[3],
           const float mmppixr[3], const float centerr[3],
           float vtneg, float vtpos, float ctneg, float ctpos,
           float dthresh,
           float *roi, float orad,
           int max_peaks, int min_vox, const float *statmask,
	   int *npeaksp, EXTREMUM **peaksp) {
    EXTREMUM *ppos, *pneg, *pall;        /**< local maxima, minima, all extrema */
    int npos = 0, nneg = 0, nall;        /**< number of maxima, minima, extrema */
    int mpos = MSIZE, mneg = MSIZE; /**< array allocation sizes */
    int i, iz, iy, ix;
    const int nx = dim[0], nxy = dim[0]*dim[1];
    const float d2thresh = dthresh * dthresh;

    if (max_peaks <= 0) {
        max_peaks = NTOP;
    }

    pf_log(PF_LOG_PEAK_TRACE,
           "starting find_peaks with orad=%g...\n",
           orad);

    ppos = malloc(mpos * sizeof(EXTREMUM));
    if (!ppos) {
        pf_error(PF_ERR_ALLOCATION, "positive extremum array");
    }

    pneg = malloc(mpos * sizeof(EXTREMUM));
    if (!pneg) {
        pf_error(PF_ERR_ALLOCATION, "negative extremum array");
    }

    pf_log(PF_LOG_PEAK_PARAMS,
           "peak value     thresholds %10.4f,%10.4f\n", vtneg, vtpos);
    pf_log(PF_LOG_PEAK_PARAMS,
           "peak curvature thresholds %10.4f,%10.4f\n", ctneg, ctpos);
    pf_log(PF_LOG_PEAK_TRACE, "compiling extrema slices\n");

    for (iz = 1; iz < dim[2] - 1; iz++) {
        for (iy = 1; iy < dim[1] - 1; iy++) {
            for (ix = 1; ix < dim[0] - 1; ix++) {
                int k;
                float del2v = 0.0, dvdx[3], d2vdx2[3], w[3], fndex[3], x[3];

                i = ix + dim[2]*(iy + dim[1]*iz);

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
                    ADD_POS_PEAK(image, dim, i);
                } else if (IS_NEG_PEAK(image, dim, i)) {
                    ADD_NEG_PEAK(image, dim, i);
                }
            }
        }
    }

    pf_log(PF_LOG_PEAK_TRACE, "\nbefore sorting npos = %d nneg = %d\n",
           npos, nneg);

    /* Sort by value, then consolidate nearby extrema */
    qsort(ppos, npos, sizeof(EXTREMUM), pcompare);
    consolidate(ppos, npos, max_peaks, d2thresh);
    qsort(pneg, nneg, sizeof(EXTREMUM), pcompare);
    consolidate(pneg, nneg, max_peaks, d2thresh);

    nall = combine_extrema(&pall, ppos, npos, pneg, nneg, max_peaks);

    pf_log(PF_LOG_PEAK_TRACE, "after consolidation nall = %d\n", nall);

    /* build a mask of spheres centered on the distinct extrema */
    if (orad > 0) {
        cull_small_extrema(image, statmask, dim, mmppixr, centerr,
                           &pall, &nall, max_peaks, orad, min_vox);
        build_mask(image, statmask, roi, dim, mmppixr, centerr,
                   pall, nall, orad);
    }

    free(ppos);
    free(pneg);

    if (npeaksp && peaksp) {
      *npeaksp = nall;
      *peaksp = pall;
    } else {
      free(pall);
    }

    pf_log(PF_LOG_PEAK_TRACE, "done.\n");
}
